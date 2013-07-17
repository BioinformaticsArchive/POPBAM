/** \file pop_sample.cpp
 *  \brief Functions for parsing samples in a BAM file
 *  \author Daniel Garrigan
 *  \version 0.3
*/

#include "popbam.h"
#include "khash.h"

KHASH_MAP_INIT_STR(sm, int);

static void add_sample_pair(bam_sample_t*, khash_t(sm)*, const char*, const char*);
static void add_pop_pair(bam_sample_t*, khash_t(sm)*, const char*, const char*);

int popbamData::bam_smpl_add(void)
{
	int n;
	const char *p;
	const char *q;
	const char *r;
	const char *s;
	kstring_t buf;
	kstring_t bug;
	khash_t(sm) *sm2id = (khash_t(sm)*)sm->sm2id;
	khash_t(sm) *pop2sm = (khash_t(sm)*)sm->pop2sm;

	p = h->text;
	n = 0;

	memset(&buf, 0, sizeof(kstring_t));
	memset(&bug, 0, sizeof(kstring_t));

	while ((q = strstr(p, "@RG")) != 0)
	{
		p = q+3;
		r = q = s = 0;
		if ((q = strstr(p, "\tID:")) != 0)
			q += 4;
		if ((r = strstr(p, "\tSM:")) != 0)
			r += 4;
		if ((s = strstr(p, "\tPO:")) != 0)
			s += 4;
		// if no PO tag is found in the header
		if (r && q && !s)
		{
			char *u;
			char *v;
			int oq;
			int or1;
			for (u=(char*)q; *u && *u != '\t' && *u != '\n'; ++u);
			for (v=(char*)r; *v && *v != '\t' && *v != '\n'; ++v);
			oq = *u;
			or1 = *v;
			*u = *v = '\0';
			buf.l = 0;
			kputs(bamfile.c_str(), &buf);
			kputc('/', &buf);
			kputs(q, &buf);
			add_sample_pair(sm, sm2id, buf.s, r);
			*u = oq;
			*v = or1;
		}
		// if PO tag is found
		else if (r && q && s)
		{
			char *u, *v, *w;
			int oq, or1, os;
			for (u=(char *)q; *u && *u != '\t' && *u != '\n'; ++u);
			for (v=(char *)r; *v && *v != '\t' && *v != '\n'; ++v);
			for (w=(char *)s; *w && *w != '\t' && *w != '\n'; ++w);
			oq = *u;
			or1 = *v;
			os = *w;
			*u = *v = *w = '\0';
			buf.l = 0;
			bug.l = 0;
			kputs(bamfile.c_str(), &buf);
			kputs(bamfile.c_str(), &bug);
			kputc('/', &buf);
			kputc('/', &bug);
			kputs(q, &buf);
			kputs(r, &bug);
			add_sample_pair(sm, sm2id, buf.s, r);
			add_pop_pair(sm, pop2sm, bug.s, s);
			*u = oq;
			*v = or1;
			*w = os;
		}
		// no read group header line is found
		else
			break;
		p = q > r ? q : r;
		p = p > s ? p : s;
		++n;
	}

	if (n == 0)
	{
		add_sample_pair(sm, sm2id, bamfile.c_str(), bamfile.c_str());
		add_pop_pair(sm, sm2id, bamfile.c_str(), bamfile.c_str());
	}

	free(buf.s);
	free(bug.s);

	return 0;
}

void popbamData::bam_smpl_init(void)
{
	sm = (bam_sample_t*)calloc(1, sizeof(bam_sample_t));
	sm->sm2popid = kh_init(sm);
	sm->rg2smid = kh_init(sm);
	sm->sm2id = kh_init(sm);
	sm->pop2sm = kh_init(sm);
}

void popbamData::bam_smpl_destroy(void)
{
	int i;
	khint_t k;
	khash_t(sm) *rg2smid = (khash_t(sm)*)sm->rg2smid;
	khash_t(sm) *sm2popid = (khash_t(sm)*)sm->sm2popid;

	if (sm == 0)
		return;

	for (i=0; i < sm->n; i++)
		free(sm->smpl[i]);
	free(sm->smpl);

	for (i=0; i < sm->npops; i++)
		free(sm->popul[i]);
	free(sm->popul);

	//dealloc strdups
	for (k=kh_begin(rg2smid); k != kh_end(rg2smid); ++k)
		if (kh_exist(rg2smid, k))
			free((char*)kh_key(rg2smid, k));

	for (k=kh_begin(sm2popid); k != kh_end(sm2popid); ++k)
		if (kh_exist(sm2popid, k))
			free((char*)kh_key(sm2popid, k));

	kh_destroy(sm, static_cast<kh_sm_t*>(sm->sm2popid));
	kh_destroy(sm, static_cast<kh_sm_t*>(sm->rg2smid));
	kh_destroy(sm, static_cast<kh_sm_t*>(sm->sm2id));
	kh_destroy(sm, static_cast<kh_sm_t*>(sm->pop2sm));
	free(sm);
}

static void add_sample_pair(bam_sample_t *sm, khash_t(sm) *sm2id, const char *key, const char *val)
{
	int ret;
	khint_t k_rg, k_sm;
	khash_t(sm) *rg2smid = (khash_t(sm)*)sm->rg2smid;

	k_rg = kh_get(sm, rg2smid, key);
	if (k_rg != kh_end(rg2smid))
		return;
	else
		k_rg = kh_put(sm, rg2smid, strdup(key), &ret);

	k_sm = kh_get(sm, sm2id, val);
	if (k_sm == kh_end(sm2id))
	{
		//need to allocate for more samples?
		if (sm->n == sm->m)
		{
			sm->m = sm->m ? sm->m << 1 : 1;
			sm->smpl = static_cast<char**>(realloc(sm->smpl, sizeof(void*)*sm->m));
		}

		//add sample name
		sm->smpl[sm->n] = strdup(val);
		
		//add new entry in sm2id hash table
		k_sm = kh_put(sm, sm2id, sm->smpl[sm->n], &ret);
		if (ret > 0)
			kh_val(sm2id, k_sm) = sm->n++;
	}

	//update rg2smid hash table
	kh_val(rg2smid, k_rg) = kh_val(sm2id, k_sm);
}

static void add_pop_pair(bam_sample_t *sm, khash_t(sm) *pop2sm, const char *key, const char *val)
{
	int ret;
	khint_t k_rg, k_sm;
	khash_t(sm) *sm2popid = (khash_t(sm)*)sm->sm2popid;

	//is the key already in the sm2popid hash table?
	k_rg = kh_get(sm, sm2popid, key);

	//duplicated @RG-ID
	if (k_rg != kh_end(sm2popid))
		return;
	else
		k_rg = kh_put(sm, sm2popid, strdup(key), &ret);

	//does val appear as a key in the pop2sm hash table?
	k_sm = kh_get(sm, pop2sm, val);

	//if val does not appears as a key
	if (k_sm == kh_end(pop2sm))
	{
		//need to allocate more samples?
		if (sm->b == sm->npops)
		{
			sm->b = sm->b ? sm->b << 1 : 1;
			sm->popul = static_cast<char**>(realloc(sm->popul, sizeof(void*)*sm->b));
		}

		//add sample name
		sm->popul[sm->npops] = strdup(val);

		//add new entry in pop2sm hash table
		k_sm = kh_put(sm, pop2sm, sm->popul[sm->npops], &ret);
		if (ret > 0)
			kh_val(pop2sm, k_sm) = sm->npops++;
	}

	//update sm2popid hash table
	kh_val(sm2popid, k_rg) = kh_val(pop2sm, k_sm);
}

int bam_smpl_rg2smid(const bam_sample_t *sm, const char *fn, const char *rg, kstring_t *str)
{
	khint_t k;
	khash_t(sm) *rg2smid = (khash_t(sm)*)sm->rg2smid;

	if (rg)
	{
		str->l = 0;
		kputs(fn, str);
		kputc('/', str);
		kputs(rg, str);
		k = kh_get(sm, rg2smid, str->s);
	}
	else
		k = kh_get(sm, rg2smid, fn);

	return k == kh_end(rg2smid) ? -1 : kh_val(rg2smid, k);
}

int bam_smpl_sm2popid(const bam_sample_t *sm, const char *fn, const char *smpl, kstring_t *str)
{
	khint_t k;
	khash_t(sm) *sm2popid = (khash_t(sm)*)sm->sm2popid;

	if (smpl)
	{
		str->l = 0;
		kputs(fn, str);
		kputc('/', str);
		kputs(smpl, str);
		k = kh_get(sm, sm2popid, str->s);
	}
	else
		k = kh_get(sm, sm2popid, fn);

	return k == kh_end(sm2popid) ? -1 : kh_val(sm2popid, k);
}
