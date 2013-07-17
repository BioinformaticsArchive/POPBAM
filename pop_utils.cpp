/** \file pop_utils.cpp
 *  \brief Utility functions evolutionary analysis of BAM files
 *  \author Daniel Garrigan
 *  \version 0.3
 *
 * Notes:
 * Byte/bit ordering of unsigned long long consensus base data type ("cb"):
 * Byte 1:     Boolean flags (cb[i]&0xff)
 *    bit 1:      Does this site pass the quality filters? (cb[i]&0x1)
 *    bit 2:      Is there a variant present at this site? (cb[i]&0x2)
 *    bit 3:      Not implemented
 *    bit 4:      Not implemented
 *    bit 5:      Not implemented
 *    bit 6:      Not implemented
 *    bit 7:      Not implemented
 *    bit 8:      Not implemented
 * Byte 2:     unsigned char-- the IUPAC consensus genotype (cb[i]>>8)&0xff
 *    bits 1-2:   unsigned char-- the IUPAC base call for allele 1 (cb[i]>>8)&0x3
 *    bits 3-4:   unsigned char-- the IUPAC base call for allele 2 (cb[i]>>10)&0x3
 *    bits 5-8:   Not implemented
 * Bytes 3-4:  unsigned short-- the number of reads mapped to that site (cb[i]>>16)&0xffff
 * Bytes 5-6:  unsigned short-- the SNP quality score (snpQ) (cb[i] >> 32)&0xffff
 * Bytes 7-8:  unsigned short-- the root-mean quality score (rmsQ) (cb[i]>>48)&0xffff
 *
**/
#include "popbam.h"
#include "ksort.h"
#include "khash.h"
#include "gamma.h"

#define M_LN2 0.69314718055994530942
#define M_LN10 2.30258509299404568402

typedef char *str_p;

KHASH_MAP_INIT_STR(s, int)
KHASH_MAP_INIT_STR(r2l, str_p)
KSORT_INIT_GENERIC(uint16_t)

void bam_init_header_hash(bam_header_t *header);

unsigned char iupac_rev[256] = {
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14, 0,14, 1, 14,14,14, 2, 14,14,14,14, 14,14,14,14,
	14,14,14,14,  3,14,14,14, 14,14,14,14, 14,14,14,14,
	14, 0,14, 1, 14,14,14, 2, 14,14,14,14, 14,14,14,14,
	14,14,14,14,  3,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14
};

unsigned short popcount64(unsigned long long x)
{
	x = (x&0x5555555555555555ULL)+((x>>1)&0x5555555555555555ULL);
	x = (x&0x3333333333333333ULL)+((x>>2)&0x3333333333333333ULL);
	x = (x&0x0F0F0F0F0F0F0F0FULL)+((x>>4)&0x0F0F0F0F0F0F0F0FULL);
	return (x*0x0101010101010101ULL)>>56;
}

unsigned int hamming_distance(unsigned long long x, unsigned long long y)
{
	unsigned int dist = 0;
	unsigned long long val = x^y;

	// count the number of set bits
	while (val)
	{
		++dist;
		val &= val-1;
	}

	return dist;
}

unsigned long long gl2cns(float q[16], unsigned short k)
{
	unsigned char i, j;
	unsigned short min_ij = 0;
	unsigned long long snp_quality;
	unsigned long long num_reads;
	unsigned long long genotype;
	float min = FLT_MAX;
	float min_next = FLT_MAX;
	float likelihood;
	std::string msg;

	for (i=0; i < NBASES; ++i)
	{
		for (j=i; j < NBASES; ++j)
		{
			likelihood = q[i<<2|j];
			if (likelihood < min)
			{
				min_ij = i<<2|j;
				min_next = min;
				min = likelihood;
			}
			else
				if (likelihood < min_next)
					min_next = likelihood;
		}
	}

	// return consensus base
	snp_quality = (unsigned long long)((min_next-min)+0.499)<<(CHAR_BIT*4);
	num_reads = (unsigned long long)k<<(CHAR_BIT*2);
	genotype = (unsigned long long)min_ij<<CHAR_BIT;

	return snp_quality+num_reads+genotype;
}

unsigned long long qfilter(int num_samples, unsigned long long *cb, int min_rmsQ, int min_depth, int max_depth)
{
	unsigned short rms;
	unsigned short num_reads;
	unsigned long long coverage = 0;

	for (int i=0; i < num_samples; ++i)
	{
		rms = (cb[i]>>(CHAR_BIT*6))&0xffff;
		num_reads = (cb[i]>>(CHAR_BIT*2))&0xffff;
		if ((rms >= min_rmsQ) && (num_reads >= min_depth) && (num_reads <= max_depth))
		{
			cb[i] |= 0x1ULL;
			coverage |= 0x1ULL<<i;
		}
	}

	return coverage;
}

int segbase(int num_samples, unsigned long long *cb, char ref, int min_snpq)
{
	int i, j, k;
	unsigned char genotype;
	unsigned char allele1;
	unsigned char allele2;
	unsigned short snp_quality;
	int baseCount[NBASES] = {0, 0, 0, 0};

	for (i=0; i < num_samples; ++i)
	{
		genotype = (cb[i]>>CHAR_BIT)&0xff;
		allele1 = (genotype>>2)&0x3;
		allele2 = genotype&0x3;
		snp_quality = (cb[i]>>(CHAR_BIT*4))&0xffff;

		// if homozygous and different from reference with high SNP quality
		if ((allele1 == allele2) && (iupac[genotype] != ref) && (snp_quality >= min_snpq))
		{
			cb[i] |= 0x2ULL;
			++baseCount[allele1];
		}
		// if SNP quality is low, revert both alleles to the reference allele
		else if ((allele1 == allele2) && (iupac[genotype] != ref) && (snp_quality < min_snpq))
		{
			cb[i] -= (genotype-iupac_rev[ref])<<CHAR_BIT;
			cb[i] -= (genotype-iupac_rev[ref])<<(CHAR_BIT+2);
		}
		else
			continue;
	}

	// check for infinite sites model
	for (i=0, j=0, k=0; i < NBASES; ++i)
	{
		if (baseCount[i] > 0)
		{
			++j;
			k = i;
		}
	}

	if (j > 1)
		return -1;
	else
		return baseCount[k];
}

void clean_heterozygotes(int num_samples, unsigned long long *cb, int ref, int min_snpq)
{
	unsigned short snp_quality;
	unsigned char genotype;
	unsigned char allele1;
	unsigned char allele2;

	for (int i=0; i < num_samples; ++i)
	{
		genotype = (cb[i]>>CHAR_BIT)&0xff;
		allele1 = (genotype>>2)&0x3;
		allele2 = genotype&0x3;
		snp_quality = (cb[i]>>(CHAR_BIT*4))&0xffff;

		// if heterozygous and high quality SNP--make homozygous derived
		if ((allele1 != allele2) && (snp_quality >= min_snpq))
		{
			if (allele1 == iupac_rev[ref])
				cb[i] += (allele2-allele1)<<(CHAR_BIT+2);
			if (allele2 == iupac_rev[ref])
				cb[i] -= (allele2-allele1)<<CHAR_BIT;
		}
		// if heterozygous but poor quality--make homozygous ancestral
		if ((allele1 != allele2) && (snp_quality < min_snpq))
		{
			if (allele1 != iupac_rev[ref])
				cb[i] += (allele2-allele1)<<(CHAR_BIT+2);
			if (allele2 != iupac_rev[ref])
				cb[i] -= (allele2-allele1)<<CHAR_BIT;
		}
	}
}

static errmod_coef_t *cal_coef(double depcorr, double eta)
{
	int k;
	int n;
	int q;
	long double sum;
	long double sum1;
	double *lC;
	errmod_coef_t *ec;

	ec = (errmod_coef_t*)calloc(1, sizeof(errmod_coef_t));
	
	// initialize ->fk
	ec->fk = (double*)calloc(256, sizeof(double));
	ec->fk[0] = 1.0;
	for (n=1; n != 256; ++n)
		ec->fk[n] = pow(1.0-depcorr, n)*(1.0-eta)+eta;
	
	// initialize ->coef
	ec->beta = (double*)calloc(256*256*64, sizeof(double));
	lC = (double*)calloc(256*256, sizeof(double));
	for (n=1; n != 256; ++n)
	{
		double lgn = LogGamma(n+1);
		for (k=1; k <= n; ++k)
			lC[n<<8|k] = lgn-LogGamma(k+1)-LogGamma(n-k+1);
	}
	for (q=1; q != 64; ++q)
	{
		double e = pow(10.0, -q/10.0);
		double le = log(e);
		double le1 = log(1.0-e);
		for (n=1; n <= 255; ++n)
		{
			double *beta = ec->beta+(q << 16 | n << 8);
			sum1 = sum = 0.0;
			for (k=n; k >= 0; --k, sum1 = sum)
			{
				sum = sum1 + expl(lC[n<<8|k] + k*le + (n-k)*le1);
				beta[k] = -10.0 / M_LN10 * logl(sum1 / sum);
			}
		}
	}

	// initialize ->lhet
	ec->lhet = (double*)calloc(256*256, sizeof(double));
	for (n=0; n < 256; ++n)
		for (k=0; k < 256; ++k)
			ec->lhet[n<<8|k] = lC[n<<8|k] - M_LN2 * n;
	free(lC);

	return ec;
}

errmod_t *errmod_init(float depcorr)
{
	errmod_t *em;

	em = (errmod_t*)calloc(1, sizeof(errmod_t));
	em->depcorr = depcorr;
	em->coef = cal_coef(depcorr, 0.03);

	return em;
}

void errmod_destroy(errmod_t *em)
{
	if (em == 0)
		return;
	free(em->coef->lhet);
	free(em->coef->fk);
	free(em->coef->beta);
	free(em->coef);
	free(em);
}

// qual:6, strand:1, base:4
int errmod_cal(const errmod_t *em, unsigned short n, int m, unsigned short *bases, float *q)
{
	call_aux_t aux;
	int i, j, k, w[32];

	if (m > m)
		return -1;
	memset(q, 0, m*m*sizeof(float));
	if (n == 0)
		return 0;

	// calculate aux.esum and aux.fsum
	// then sample 255 bases
	if (n > 255)
	{
		ks_shuffle(uint16_t, n, bases);
		n = 255;
	}
	ks_introsort(uint16_t, n, bases);
	memset(w, 0, 32*sizeof(int));
	memset(&aux, 0, sizeof(call_aux_t));

	// calculate esum and fsum
	for (j=n-1; j >= 0; --j)
	{
		unsigned short b = bases[j];
		int q = b >> 5 < NBASES ? NBASES : b >> 5;
		if (q > 63)
			q = 63;
		k = b&0x1f;
		aux.fsum[k&0xf] += em->coef->fk[w[k]];
		aux.bsum[k&0xf] += em->coef->fk[w[k]] * em->coef->beta[q<<16|n<<8|aux.c[k&0xf]];
		++aux.c[k&0xf];
		++w[k];
	}
	// generate likelihood
	for (j=0; j != m; ++j)
	{
		float tmp1, tmp3;
		int tmp2, bar_e;
		// homozygous
		for (k=0, tmp1=tmp3=0.0, tmp2=0; k != m; ++k)
		{
			if (k == j)
				continue;
			tmp1 += aux.bsum[k];
			tmp2 += aux.c[k];
			tmp3 += aux.fsum[k];
		}
		if (tmp2)
		{
			bar_e = (int)(tmp1/tmp3+0.499);
			if (bar_e > 63)
				bar_e = 63;
			q[j*m+j] = tmp1;
		}
		// heterozygous
		for (k=j+1; k < m; ++k)
		{
			int cjk = aux.c[j]+aux.c[k];
			for (i=0, tmp2=0, tmp1=tmp3=0.0; i < m; ++i)
			{
				if ((i == j) || (i == k))
					continue;
				tmp1 += aux.bsum[i];
				tmp2 += aux.c[i];
				tmp3 += aux.fsum[i];
			}
			if (tmp2)
			{
				bar_e = (int)(tmp1/tmp3+0.499);
				if (bar_e > 63)
					bar_e = 63;
				q[j*m+k] = q[k*m+j] = -4.343*em->coef->lhet[cjk<<8|aux.c[k]]+tmp1;
			}
			// all the bases are either j or k
			else
				q[j*m+k] = q[k*m+j] = -4.343*em->coef->lhet[cjk<<8|aux.c[k]];
		}
		for (k=0; k != m; ++k)
			if (q[j*m+k] < 0.0)
				q[j*m+k] = 0.0;
	}

	return 0;
}

void bam_init_header_hash(bam_header_t *header)
{
	int i;

	if (header->hash == 0)
	{
		int ret = 0;
		khiter_t iter;
		khash_t(s) *h;

		header->hash = h = kh_init(s);
		for (i=0; i < header->n_targets; ++i)
		{
			iter = kh_put(s, h, header->target_name[i], &ret);
			kh_value(h, iter) = i;
		}
	}
}

int bam_parse_region(bam_header_t *header, std::string region, int *ref_id, int *beg, int *end)
{
	std::size_t l;
	std::size_t name_end;
	khiter_t iter;
	khash_t(s) *h;

	bam_init_header_hash(header);
	h = (khash_t(s)*)header->hash;
	*ref_id = *beg = *end = -1;
	name_end = l = region.length();

	// remove spaces and commas
	std::string::iterator end_pos = std::remove(region.begin(), region.end(), ' ');
	region.erase(end_pos, region.end());
	end_pos = std::remove(region.begin(), region.end(), ',');
	region.erase(end_pos, region.end());
	l = region.length();

	// determine the sequence name
	// look for colon from the end
	name_end = region.find(":");
	if (name_end == std::string::npos)
		name_end = l;

	// check if this is really the end
	if (name_end < l)
	{
		std::string coords = region.substr(name_end+1);
		std::size_t n_hyphen = std::count(coords.begin(), coords.end(), '-');
		std::size_t n_nondigits = coords.find_first_not_of("0123456789,-");

		// malformated region string; then take str as the name
		if ((n_nondigits != std::string::npos) || (n_hyphen > 1))
			name_end = l;
		std::string scaffold_name = region.substr(0, name_end);
		iter = kh_get(s, h, scaffold_name.c_str());

		// cannot find the sequence name
		if (iter == kh_end(h))
		{
			// try str as the name
			iter = kh_get(s, h, region.c_str());
			if (iter == kh_end(h))
			{
				std::cerr << "Cannot find sequence name " << region << " in header" << std::endl;
				return -1;
			}
		}
	}
	else
		iter = kh_get(s, h, region.c_str());
	if (iter == kh_end(h))
		return -1;
	*ref_id = kh_val(h, iter);

	// parse the interval
	if (name_end < l)
	{
		std::string coords = region.substr(name_end+1);
		std::size_t parse = coords.find("-");
		std::string first = coords.substr(0, parse);
		*beg = atoi(first.c_str());
		if (*beg > 0)
			--*beg;
		std::string last = coords.substr(parse+1);
		*end = atoi(last.c_str());
	}
	else
	{
		*beg = 0;
		*end = header->target_len[*ref_id];
	}

	return *beg <= *end ? 0 : -1;
}

char *get_refid(char *htext)
{
	char *u;
	char *v;
	char *w;
	const int idblock = 200;
	int z;
	char *refid = NULL;

	u = htext;
	v = strstr(htext, "AS:");

	if (!v)
	{
		std::string msg("Unable to parse reference sequence name\nBe sure the AS tag is defined in the sequence dictionary");
		fatal_error(msg, __FILE__, __LINE__, 0);
	}

	u = v+3;

	for (z=0, w=(char*)u; *w && *w != '\t' && *w != '\n'; ++w, ++z);

	try
	{
		refid = new char [idblock];
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
	}
	refid[0] = '\0';
	strncpy(refid, u, z);
	refid[z] = '\0';

	return refid;
}

int fetch_func(const bam1_t *b, void *data)
{
	bam_plbuf_t *buf;

	buf = (bam_plbuf_t*)data;
	bam_plbuf_push(b, buf);

	return 0;
}

void fatal_error(std::string msg, const char* file, int line, void (*err_func)(void))
{
	std::cerr << "popbam runtime error:" << std::endl;
	std::cerr << msg << std::endl;
	std::cerr << "In " << file << " on line " << line << std::endl;
	if (err_func)
		err_func();
	std::cerr << "Exiting program" << std::endl;
	exit(EXIT_FAILURE);
}
