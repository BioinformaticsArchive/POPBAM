#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "sam.h"

#ifdef _MSC_VER
#define inline __inline
#endif

typedef struct {
	int k;
	int x;
	int y;
	int end;
} cstate_t;

static cstate_t g_cstate_null = { -1, 0, 0, 0 };

typedef struct __linkbuf_t {
	bam1_t b;
	uint32_t beg;
	uint32_t end;
	cstate_t s;
	struct __linkbuf_t *next;
} lbnode_t;

/* --- BEGIN: Memory pool */

typedef struct {
	int cnt;
	int n;
	int max;
	lbnode_t **buf;
} mempool_t;

static mempool_t *mp_init(void)
{
	mempool_t *mp;
	mp = (mempool_t*)calloc(1, sizeof(mempool_t));
	return mp;
}

static void mp_destroy(mempool_t *mp)
{
	int k;

	for (k=0; k < mp->n; ++k)
	{
		free(mp->buf[k]->b.data);
		free(mp->buf[k]);
	}
	free(mp->buf);
	free(mp);
}

static inline lbnode_t *mp_alloc(mempool_t *mp)
{
	++mp->cnt;
	if (mp->n == 0)
		return (lbnode_t*)calloc(1, sizeof(lbnode_t));
	else
		return mp->buf[--mp->n];
}

static inline void mp_free(mempool_t *mp, lbnode_t *p)
{
	--mp->cnt;
	// clear lbnode_t::next here
	p->next = 0;
	if (mp->n == mp->max)
	{
		mp->max = mp->max ? mp->max<<1 : 256;
		mp->buf = (lbnode_t**)realloc(mp->buf, sizeof(lbnode_t*) * mp->max);
	}
	mp->buf[mp->n++] = p;
}

/* --- END: Memory pool */

/* --- BEGIN: Auxiliary functions */

/* s->k: the index of the CIGAR operator that has just been processed.
   s->x: the reference coordinate of the start of s->k
   s->y: the query coordiante of the start of s->k
 */
static inline int resolve_cigar2(bam_pileup1_t *p, uint32_t pos, cstate_t *s)
{
#define _cop(c) ((c)&BAM_CIGAR_MASK)
#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)

	bam1_t *b = p->b;
	bam1_core_t *c = &b->core;
	uint32_t *cigar = bam1_cigar(b);
	int k;
	int is_head = 0;

	// determine the current CIGAR operation
	// fprintf(stderr, "%s\tpos=%d\tend=%d\t(%d,%d,%d)\n", bam1_qname(b), pos, s->end, s->k, s->x, s->y);
	// never processed
	if (s->k == -1)
	{
		is_head = 1;
		
		// just one operation, save a loop
		if (c->n_cigar == 1)
		{
		  if ((_cop(cigar[0]) == BAM_CMATCH) || (_cop(cigar[0]) == BAM_CEQUAL) || (_cop(cigar[0]) == BAM_CDIFF))
			  s->k = 0, s->x = c->pos, s->y = 0;
		}
		// find the first match or deletion
		else
		{
			for (k=0, s->x=c->pos, s->y=0; k < c->n_cigar; ++k)
			{
				int op = _cop(cigar[k]);
				int l = _cln(cigar[k]);
				if ((op == BAM_CMATCH) || (op == BAM_CDEL) || (op == BAM_CEQUAL) || (op == BAM_CDIFF))
					break;
				else if (op == BAM_CREF_SKIP)
					s->x += l;
				else if ((op == BAM_CINS) || (op == BAM_CSOFT_CLIP))
					s->y += l;
			}
			assert(k < c->n_cigar);
			s->k = k;
		}
	}
	// the read has been processed before
	else
	{
		int op;
		int l = _cln(cigar[s->k]);

		// jump to the next operation
		if (pos - s->x >= l)
		{
			// otherwise a bug: this function should not be called in this case
			assert(s->k < c->n_cigar);
			op = _cop(cigar[s->k+1]);
			// jump to the next without a loop
			if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF)
			{
			  if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF)
				  s->y += l;
				s->x += l;
				++s->k;
			}
			// find the next M/D/N/=/X
			else
			{
			  if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF)
				  s->y += l;
				s->x += l;
				for (k=s->k+1; k < c->n_cigar; ++k)
				{
					op = _cop(cigar[k]), l = _cln(cigar[k]);
					if ((op == BAM_CMATCH) || (op == BAM_CDEL) || (op == BAM_CREF_SKIP) || (op == BAM_CEQUAL) || (op == BAM_CDIFF))
						break;
					else if ((op == BAM_CINS) || (op == BAM_CSOFT_CLIP))
						s->y += l;
				}
				s->k = k;
			}
			// otherwise a bug
			assert(s->k < c->n_cigar);
		} // else, do nothing
	}
	{
		// collect pileup information
		int op;
		int l;
		op = _cop(cigar[s->k]);
		l = _cln(cigar[s->k]);
		p->is_del = p->indel = p->is_refskip = 0;
		// peek the next operation
		if ((s->x + l - 1 == pos) && (s->k + 1 < c->n_cigar))
		{
			int op2 = _cop(cigar[s->k+1]);
			int l2 = _cln(cigar[s->k+1]);
			if (op2 == BAM_CDEL)
				p->indel = -(int)l2;
			else if (op2 == BAM_CINS)
				p->indel = l2;
			// no working for adjacent padding
			else if ((op2 == BAM_CPAD) && (s->k + 2 < c->n_cigar))
			{
				int l3 = 0;
				for (k=s->k+2; k < c->n_cigar; ++k)
				{
					op2 = _cop(cigar[k]);
					l2 = _cln(cigar[k]);
					if (op2 == BAM_CINS)
						l3 += l2;
					else if (op2 == BAM_CDEL || op2 == BAM_CMATCH || op2 == BAM_CREF_SKIP || op2 == BAM_CEQUAL || op2 == BAM_CDIFF)
						break;
				}
				if (l3 > 0)
					p->indel = l3;
			}
		}
		if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF))
		{
			p->qpos = s->y + (pos - s->x);
		}
		else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP))
		{
			p->is_del = 1; 
			// FIXME: distinguish D and N!!!!!
			p->qpos = s->y;
			p->is_refskip = (op == BAM_CREF_SKIP);
		}
		// cannot be other operations; otherwise a bug
		p->is_head = (pos == c->pos);
		p->is_tail = (pos == s->end);
	}
	return 1;
}

/* --- END: Auxiliary functions */

/*******************
 * pileup iterator *
 *******************/

struct __bam_plp_t {
	mempool_t *mp;
	lbnode_t *head;
	lbnode_t *tail;
	lbnode_t *dummy;
	int32_t tid;
	int32_t pos;
	int32_t max_tid;
	int32_t max_pos;
	int is_eof;
	int flag_mask;
	int max_plp;
	int error;
	int maxcnt;
	bam_pileup1_t *plp;
	// for the "auto" interface only
	bam1_t *b;
	bam_plp_auto_f func;
	void *data;
};

bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data)
{
	bam_plp_t iter;
	iter = (bam_plp_t)calloc(1, sizeof(struct __bam_plp_t));
	iter->mp = mp_init();
	iter->head = iter->tail = mp_alloc(iter->mp);
	iter->dummy = mp_alloc(iter->mp);
	iter->max_tid = iter->max_pos = -1;
	iter->flag_mask = BAM_DEF_MASK;
	iter->maxcnt = 8000;
	if (func)
	{
		iter->func = func;
		iter->data = data;
		iter->b = bam_init1();
	}
	return iter;
}

void bam_plp_destroy(bam_plp_t iter)
{
	mp_free(iter->mp, iter->dummy);
	mp_free(iter->mp, iter->head);
	if (iter->mp->cnt != 0)
		fprintf(stderr, "[bam_plp_destroy] memory leak: %d. Continue anyway.\n", iter->mp->cnt);
	mp_destroy(iter->mp);
	if (iter->b)
		bam_destroy1(iter->b);
	free(iter->plp);
	free(iter);
}

const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
	if (iter->error)
	{
		*_n_plp = -1;
		return 0;
	}
	*_n_plp = 0;
	if (iter->is_eof && (iter->head->next == 0))
		return 0;
	while (iter->is_eof || (iter->max_tid > iter->tid) || ((iter->max_tid == iter->tid) && (iter->max_pos > iter->pos)))
	{
		int n_plp = 0;
		lbnode_t *p;
		lbnode_t *q;
		// write iter->plp at iter->pos
		iter->dummy->next = iter->head;
		for (p=iter->head, q=iter->dummy; p->next; q = p, p = p->next)
		{
			// then remove
			if ((p->b.core.tid < iter->tid) || ((p->b.core.tid == iter->tid) && (p->end <= iter->pos)))
			{
				q->next = p->next; mp_free(iter->mp, p); p = q;
			}
			// here: p->end > pos; then add to pileup
			else if ((p->b.core.tid == iter->tid) && (p->beg <= iter->pos))
			{
				// then double the capacity
				if (n_plp == iter->max_plp)
				{
					iter->max_plp = iter->max_plp ? iter->max_plp<<1 : 256;
					iter->plp = (bam_pileup1_t*)realloc(iter->plp, sizeof(bam_pileup1_t) * iter->max_plp);
				}
				iter->plp[n_plp].b = &p->b;
				// actually always true...
				if (resolve_cigar2(iter->plp + n_plp, iter->pos, &p->s))
					++n_plp;
			}
		}
		// dummy->next may be changed
		iter->head = iter->dummy->next;
		*_n_plp = n_plp;
		*_tid = iter->tid;
		*_pos = iter->pos;
		// update iter->tid and iter->pos
		if (iter->head->next)
		{
			if (iter->tid > iter->head->b.core.tid)
			{
				fprintf(stderr, "[%s] unsorted input. Pileup aborts.\n", __FUNCTION__);
				iter->error = 1;
				*_n_plp = -1;
				return 0;
			}
		}
		// come to a new reference sequence
		if (iter->tid < iter->head->b.core.tid)
		{
			iter->tid = iter->head->b.core.tid;
			// jump to the next reference
			iter->pos = iter->head->beg;
		}
		// here: tid == head->b.core.tid
		else if (iter->pos < iter->head->beg)
		{
			// jump to the next position
			iter->pos = iter->head->beg;
		}
		// scan contiguously
		else
			++iter->pos;
		// return
		if (n_plp)
			return iter->plp;
		if (iter->is_eof && (iter->head->next == 0))
			break;
	}
	return 0;
}

int bam_plp_push(bam_plp_t iter, const bam1_t *b)
{
	if (iter->error)
		return -1;
	if (b)
	{
		if (b->core.tid < 0)
			return 0;
		if (b->core.flag & iter->flag_mask)
			return 0;
		if (iter->tid == b->core.tid && iter->pos == b->core.pos && iter->mp->cnt > iter->maxcnt)
			return 0;
		bam_copy1(&iter->tail->b, b);
		iter->tail->beg = b->core.pos;
		iter->tail->end = bam_calend(&b->core, bam1_cigar(b));
		iter->tail->s = g_cstate_null;

		// initialize cstate_t
		iter->tail->s.end = iter->tail->end - 1;
		if (b->core.tid < iter->max_tid)
		{
			fprintf(stderr, "[bam_pileup_core] the input is not sorted (chromosomes out of order)\n");
			iter->error = 1;
			return -1;
		}
		if ((b->core.tid == iter->max_tid) && (iter->tail->beg < iter->max_pos))
		{
			fprintf(stderr, "[bam_pileup_core] the input is not sorted (reads out of order)\n");
			iter->error = 1;
			return -1;
		}
		iter->max_tid = b->core.tid;
		iter->max_pos = iter->tail->beg;
		if ((iter->tail->end > iter->pos) || (iter->tail->b.core.tid > iter->tid))
		{
			iter->tail->next = mp_alloc(iter->mp);
			iter->tail = iter->tail->next;
		}
	}
	else
		iter->is_eof = 1;
	return 0;
}

const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
	const bam_pileup1_t *plp;
	
	if ((iter->func == 0) || iter->error)
	{
		*_n_plp = -1;
		return 0;
	}
	if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0)
		return plp;
	// no pileup line can be obtained; read alignments
	else
	{
		*_n_plp = 0;
		if (iter->is_eof)
			return 0;
		while (iter->func(iter->data, iter->b) >= 0)
		{
			if (bam_plp_push(iter, iter->b) < 0)
			{
				*_n_plp = -1;
				return 0;
			}
			if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0)
				return plp;
			// otherwise no pileup line can be returned; read the next alignment.
		}
		bam_plp_push(iter, 0);
		if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0)
			return plp;
		return 0;
	}
}

void bam_plp_reset(bam_plp_t iter)
{
	lbnode_t *p;
	lbnode_t *q;

	iter->max_tid = iter->max_pos = -1;
	iter->tid = iter->pos = 0;
	iter->is_eof = 0;
	for (p=iter->head; p->next;)
	{
		q = p->next;
		mp_free(iter->mp, p);
		p = q;
	}
	iter->head = iter->tail;
}

void bam_plp_set_mask(bam_plp_t iter, int mask)
{
	iter->flag_mask = mask < 0 ? BAM_DEF_MASK : (BAM_FUNMAP | mask);
}

void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt)
{
	iter->maxcnt = maxcnt;
}

/*****************
 * callback APIs *
 *****************/

int bam_pileup_file(bamFile fp, int mask, bam_pileup_f func, void *func_data)
{
	bam_plbuf_t *buf;
	int ret;
	bam1_t *b;

	b = bam_init1();
	buf = bam_plbuf_init(func, func_data);
	bam_plbuf_set_mask(buf, mask);

	while ((ret = bam_read1(fp, b)) >= 0)
		bam_plbuf_push(b, buf);

	bam_plbuf_push(0, buf);
	bam_plbuf_destroy(buf);
	bam_destroy1(b);

	return 0;
}

void bam_plbuf_set_mask(bam_plbuf_t *buf, int mask)
{
	bam_plp_set_mask(buf->iter, mask);
}

void bam_plbuf_reset(bam_plbuf_t *buf)
{
	bam_plp_reset(buf->iter);
}

bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data)
{
	bam_plbuf_t *buf;

	buf = (bam_plbuf_t*)calloc(1, sizeof(bam_plbuf_t));
	buf->iter = bam_plp_init(0, 0);
	buf->func = func;
	buf->data = data;

	return buf;
}

void bam_plbuf_destroy(bam_plbuf_t *buf)
{
	bam_plp_destroy(buf->iter);
	free(buf);
}

int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf)
{
	int ret;
	int n_plp;
	int tid;
	int pos;
	const bam_pileup1_t *plp;
	
	ret = bam_plp_push(buf->iter, b);

	if (ret < 0)
		return ret;

	while ((plp = bam_plp_next(buf->iter, &tid, &pos, &n_plp)) != 0)
		buf->func(tid, pos, n_plp, plp, buf->data);

	return 0;
}

/***********
 * mpileup *
 ***********/

struct __bam_mplp_t {
	int n;
	uint64_t min, *pos;
	bam_plp_t *iter;
	int *n_plp;
	const bam_pileup1_t **plp;
};

bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data)
{
	int i;
	bam_mplp_t iter;

	iter = (bam_mplp_t)calloc(1, sizeof(struct __bam_mplp_t));
	iter->pos = (uint64_t*)calloc((size_t)n, 8);
	iter->n_plp = (int*)calloc((size_t)n, sizeof(int));
	iter->plp = (const bam_pileup1_t**)calloc((size_t)n, sizeof(const bam_pileup1_t*));
	iter->iter = (bam_plp_t*)calloc((size_t)n, sizeof(bam_plp_t));
	iter->n = n;
	iter->min = (uint64_t)-1;

	for (i=0; i < n; ++i)
	{
		iter->iter[i] = bam_plp_init(func, data[i]);
		iter->pos[i] = iter->min;
	}

	return iter;
}

void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt)
{
	int i;

	for (i=0; i < iter->n; ++i)
		iter->iter[i]->maxcnt = maxcnt;
}

void bam_mplp_destroy(bam_mplp_t iter)
{
	int i;

	for (i=0; i < iter->n; ++i)
		bam_plp_destroy(iter->iter[i]);
	free(iter->iter);
	free(iter->pos);
	free(iter->n_plp);
	free(iter->plp);
	free(iter);
}

int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp)
{
	int i;
	int ret = 0;
	uint64_t new_min = (uint64_t)-1;

	for (i=0; i < iter->n; ++i)
	{
		if (iter->pos[i] == iter->min)
		{
			int tid;
			int pos;
			iter->plp[i] = bam_plp_auto(iter->iter[i], &tid, &pos, &iter->n_plp[i]);
			iter->pos[i] = (uint64_t)tid<<32 | pos;
		}
		if (iter->plp[i] && (iter->pos[i] < new_min))
			new_min = iter->pos[i];
	}
	iter->min = new_min;
	if (new_min == (uint64_t)-1)
		return 0;
	*_tid = new_min>>32;
	*_pos = (uint32_t)new_min;
	for (i=0; i < iter->n; ++i)
	{
		if (iter->pos[i] == iter->min)
		{
			// FIXME: valgrind reports "uninitialised value(s) at this line"
			n_plp[i] = iter->n_plp[i], plp[i] = iter->plp[i];
			++ret;
		}
		else
		{
			n_plp[i] = 0;
			plp[i] = 0;
		}
	}

	return ret;
}