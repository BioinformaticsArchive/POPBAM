/** \file popbam.cpp
 *  \brief Main entry point for evolutionary analysis of BAM files
 *  \author Daniel Garrigan
 *  \version 0.3
*/
#include "popbam.h"
#include "tables.h"

int bam_nt16_nt4_table[16] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

const char iupac[16] = {'A','M','R','W','N','C','S','Y','N','N','G','K','N','N','N','T'};

unsigned char bam_nt16_table[256] =
{
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

unsigned char iupac_rev[256] =
{
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

int main(int argc, char *argv[])
{
	if (argc < 2)
		return popbam_usage();
	if (!strcmp(argv[1], "snp"))
		return main_snp(argc-1, argv+1);
	else if (!strcmp(argv[1], "haplo"))
		return main_haplo(argc-1, argv+1);
	else if (!strcmp(argv[1], "diverge"))
		return main_diverge(argc-1, argv+1);
	else if (!strcmp(argv[1], "tree"))
		return main_tree(argc-1, argv+1);
	else if (!strcmp(argv[1], "nucdiv"))
		return main_nucdiv(argc-1, argv+1);
	else if (!strcmp(argv[1], "ld"))
		return main_ld(argc-1, argv+1);
	else if (!strcmp(argv[1], "sfs"))
		return main_sfs(argc-1, argv+1);
	else
	{
		std::cerr << "Error: unrecognized command: " << argv[1] << std::endl;
		return 1;
	}
	return 0;
}

popbamData::popbamData(void)
{
	flag = 0x0;
	num_sites = 0;
	tid = -1;
	beg = 0;
	end = 0x7fffffff;
	min_rmsQ = 25;
	min_snpQ = 25;
	min_depth = 3;
	max_depth = 255;
	min_mapQ = 13;
	min_baseQ = 13;
	het_prior = 0.0001;
}

void popbamData::checkBAM(void)
{
	std::string msg;

	bam_in = samopen(bamfile.c_str(), "rb", 0);

	// check if BAM file is readable
	if (!bam_in)
	{
		msg = "Cannot read BAM file " + bamfile;
		fatal_error(msg, __FILE__, __LINE__, 0);
	}

	// check if BAM header is returned
	if (!bam_in->header)
	{
		msg = "Cannot read BAM header from file " + bamfile;
		fatal_error(msg, __FILE__, __LINE__, 0);
	}
	else
		h = bam_in->header;

	//read in new header text
	if (flag & BAM_HEADERIN)
	{
		std::ifstream headin(headfile);
		headin.seekg(0, std::ios::end);
		h->l_text = headin.tellg();
		headin.seekg(0, std::ios::beg);
		h->text = (char*)realloc(h->text, (size_t)h->l_text);
		headin.read(h->text, h->l_text);
		headin.close();
	}

	// check for bam index file
	if (!(idx = bam_index_load(bamfile.c_str())))
	{
		msg = "Index file not available for BAM file " + bamfile;
		fatal_error(msg, __FILE__, __LINE__, 0);
	}

	// check if fastA reference index is available
	fai_file = fai_load(reffile.c_str());
	if (!fai_file)
	{
		msg = "Failed to load index for fastA reference file: " + reffile;
		fatal_error(msg, __FILE__, __LINE__, 0);
	}
}

void popbamData::assign_pops(void)
{
	int si = -1;
	kstring_t buf;

	memset(&buf, 0, sizeof(kstring_t));

	for (int i=0; i < sm->n; i++)
	{
		if (sm->smpl[i])
			si = bam_smpl_sm2popid(sm, bamfile.c_str(), sm->smpl[i], &buf);

		if (si < 0)
			si = bam_smpl_sm2popid(sm, bamfile.c_str(), 0, &buf);

		if (si < 0)
		{
			std::string msg;
			std::string missing_sample(sm->smpl[i]);
			msg = "Sample " + missing_sample + " not assigned to a population.\nPlease check BAM header file definitions";
			fatal_error (msg, __FILE__, __LINE__, 0);
		}

		pop_mask[si] |= 0x1ULL<<i;
		pop_nsmpl[si]++;
	}
}

unsigned long long popbamData::cal_site_type(unsigned long long *cb)
{
	unsigned long long site_type = 0;

	for (int i=0; i < sm->n; i++)
	{
		if ((cb[i]&0x3ULL) == 0x3ULL)
			site_type |= 0x1ULL<<i;
	}

	return site_type;
}

void popbamData::call_base(int n, const bam_pileup1_t *pl, unsigned long long *cb)
{
	int i, j;
	int qq;
	int baseQ;
	int tmp_baseQ;
	int b;
	int si = -1;
	int rmsq;
	const int n_smpl = sm->n;
	unsigned short k;
	unsigned short *bases = NULL;
	unsigned long long rms = 0;
	std::vector<int> depth(sm->n, 0);
	unsigned char *s = NULL;
	float q[16];
	std::string msg;
	bam_pileup1_t ***p = NULL;
	kstring_t buf;

	try
	{
		p = new bam_pileup1_t** [n_smpl];
		for (i=0; i < n_smpl; i++)
			p[i] = new bam_pileup1_t* [max_depth];
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
	}

	memset(&buf, 0, sizeof(kstring_t));

	// partition pileup according to sample
	for (i=0; i < n; i++)
	{
		if ((pl+i)->is_del || (pl+i)->is_refskip || ((pl+i)->b->core.flag & BAM_FUNMAP))
			continue;
		s = bam_aux_get((pl+i)->b, "RG");

		// skip reads with no read group tag
		if (!s)
			continue;
		else
			si = bam_smpl_rg2smid(sm, bamfile.c_str(), (char*)(s+1), &buf);
		if (si < 0)
			si = bam_smpl_rg2smid(sm, bamfile.c_str(), 0, &buf);

		if (si < 0)
		{
			free(buf.s);
			std::string rogue_rg(bam_aux2Z(s));
			msg = "Problem assigning read group " + rogue_rg + " to a sample.\nPlease check BAM header for correct SM and PO tags";
			fatal_error(msg, __FILE__, __LINE__, 0);
		}

		if (depth[si] < max_depth)
		{
			p[si][depth[si]] = const_cast<bam_pileup1_t*>(pl+i);
			depth[si]++;
		}
		else
			continue;
	}

	// fill in the base array
	for (j=0; j < n_smpl; ++j)
	{
		rmsq = 0;
		if (depth[j] > 0)
		{
			try
			{
				bases = new unsigned short [depth[j]]();
			}
			catch (std::bad_alloc& ba)
			{
				std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
			}

			for (i=k=0; i < depth[j]; ++i)
			{
				tmp_baseQ = bam1_qual(p[j][i]->b)[p[j][i]->qpos];
				if (flag & BAM_ILLUMINA)
					baseQ = tmp_baseQ > 31 ? tmp_baseQ - 31 : 0;
				else
					baseQ = tmp_baseQ;
				assert(baseQ >= 0);
				if ((baseQ < min_baseQ) || (p[j][i]->b->core.qual < min_mapQ))
					continue;
				b = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p[j][i]->b), p[j][i]->qpos)];
				if (b > 3)
					continue;
				qq = baseQ < p[j][i]->b->core.qual ? baseQ : p[j][i]->b->core.qual;
				if (qq < 4)
					qq = 4;
				if (qq > 63)
					qq = 63;
				bases[k++] = qq << 5 | (unsigned short)bam1_strand(p[j][i]->b) << 4 | b;
				rmsq += SQ(p[j][i]->b->core.qual);
			}

			// calculate genotype likelihoods
			errmod_cal(em, k, NBASES, bases, q);

			// finalize root mean quality score
			rms = (unsigned long long)(sqrt((float)(rmsq)/k)+0.499);

			// get consensus base call
			cb[j] = gl2cns(q, k);

			// add root-mean map quality score to cb array
			cb[j] |= rms<<(CHAR_BIT*6);

			// take out some garbage
			delete [] bases;
			bases = NULL;
		}
		else
			continue;
	}

	// take out garbage
	for (i=0; i < n_smpl; i++)
		delete [] p[i];
	delete [] p;
	free(buf.s);
}

int popbam_usage(void)
{
	std::cerr << std::endl;
	std::cerr << "Program: popbam " << std::endl;
	std::cerr << "(Tools to perform evolutionary analysis from BAM files)" << std::endl;
	std::cerr << "Version: " << POPBAM_RELEASE << std::endl; //"  (Revision: " << svn_version() << ")" << std::endl << std::endl;
	std::cerr << "Usage: popbam <command> [options] <in.bam> [region]"  << std::endl << std::endl;
	std::cerr << "Commands:  snp       output consensus base calls" << std::endl;
	std::cerr << "           haplo     output haplotype-based analyses" << std::endl;
	std::cerr << "           diverge   output divergence from reference" << std::endl;
	std::cerr << "           tree      output neighbor-joining trees" << std::endl;
	std::cerr << "           nucdiv    output nucleotide diversity statistics" << std::endl;
	std::cerr << "           ld        output linkage disequilibrium analysis" << std::endl;
	std::cerr << "           sfs       output site frequency spectrum analysis" << std::endl;
	std::cerr << std::endl;
	return 1;
}
