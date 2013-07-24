#include <ctype.h>
#include "bam.h"
#include "khash.h"

typedef char *str_p;

KHASH_MAP_INIT_STR(s, int)
KHASH_MAP_INIT_STR(r2l, str_p)

void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data)
{
    int ori_len = b->data_len;

    b->data_len += 3+len;
    b->l_aux += 3+len;

    if (b->m_data < b->data_len)
    {
        b->m_data = b->data_len;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
    }

    b->data[ori_len] = tag[0];
    b->data[ori_len+1] = tag[1];
    b->data[ori_len+2] = type;
    memcpy(b->data+ori_len+3, data, len);
}

uint8_t *bam_aux_get_core(bam1_t *b, const char tag[2])
{
    return bam_aux_get(b, tag);
}

#define __skip_tag(s) do { \
		int type = toupper(*(s)); \
		++(s); \
		if ((type == 'Z') || (type == 'H')) \
		{ \
			while (*(s)) \
				++(s); \
			++(s); \
		} \
		else if (type == 'B') \
			(s) += 5 + bam_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
		else \
			(s) += bam_aux_type2size(type); \
	} while(0)

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
    uint8_t *s;
    int y = tag[0]<<8 | tag[1];

    s = bam1_aux(b);
    while ((s < b->data) + b->data_len)
    {
        int x = (int)s[0]<<8 | s[1];
        s += 2;
        if (x == y)
            return s;
        __skip_tag(s);
    }

    return 0;
}

// s MUST BE returned by bam_aux_get()
int bam_aux_del(bam1_t *b, uint8_t *s)
{
    uint8_t *p;
    uint8_t *aux;

    aux = bam1_aux(b);
    p = s-2;
    __skip_tag(s);
    memmove(p, s, b->l_aux - (s-aux));
    b->data_len -= s-p;
    b->l_aux -= s-p;

    return 0;
}

int bam_aux_drop_other(bam1_t *b, uint8_t *s)
{
    if (s)
    {
        uint8_t *p;
        uint8_t *aux;

        aux = bam1_aux(b);
        p = s - 2;
        __skip_tag(s);
        memmove(aux, p, s-p);
        b->data_len -= b->l_aux - (s-p);
        b->l_aux = s-p;
    }
    else
    {
        b->data_len -= b->l_aux;
        b->l_aux = 0;
    }

    return 0;
}

void bam_destroy_header_hash(bam_header_t *header)
{
    if (header->hash)
        kh_destroy(s, (khash_t(s)*)header->hash);
}

int32_t bam_get_tid(const bam_header_t *header, const char *seq_name)
{
    khint_t k;
    khash_t(s) *h = (khash_t(s)*)header->hash;

    k = kh_get(s, h, seq_name);

    return k == kh_end(h) ? -1 : kh_value(h, k);
}

int32_t bam_aux2i(const uint8_t *s)
{
    int type;

    if (s == 0)
        return 0;
    type = *s++;
    if (type == 'c')
        return (int32_t)*(int8_t*)s;
    else if (type == 'C')
        return (int32_t)*(uint8_t*)s;
    else if (type == 's')
        return (int32_t)*(int16_t*)s;
    else if (type == 'S')
        return (int32_t)*(uint16_t*)s;
    else if ((type == 'i') || (type == 'I'))
        return *(int32_t*)s;
    else
        return 0;
}

float bam_aux2f(const uint8_t *s)
{
    int type = *s++;

    if (s == 0)
        return 0.0;
    if (type == 'f')
        return *(float*)s;
    else
        return 0.0;
}

double bam_aux2d(const uint8_t *s)
{
    int type = *s++;

    if (s == 0)
        return 0.0;
    if (type == 'd')
        return *(double*)s;
    else
        return 0.0;
}

char bam_aux2A(const uint8_t *s)
{
    int type = *s++;

    if (s == 0)
        return 0;
    if (type == 'A')
        return *(char*)s;
    else
        return 0;
}

char *bam_aux2Z(const uint8_t *s)
{
    int type = *s++;

    if (s == 0)
        return 0;
    if ((type == 'Z') || (type == 'H'))
        return (char*)s;
    else
        return 0;
}

