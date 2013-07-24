/** \file pop_haplo.cpp
 *  \brief Functions for calculating haplotype-based statistics
 *  \author Daniel Garrigan
 *  \version 0.3
*/

#include "pop_haplo.h"
#include "tables.h"

int main_haplo(int argc, char *argv[])
{
	int chr;                  //! chromosome identifier
	int beg;                  //! beginning coordinate for analysis
	int end;                  //! end coordinate for analysis
	int ref;                  //! ref
	long num_windows;         //! number of windows
	std::string msg;          //! string for error message
	bam_plbuf_t *buf;         //! pileup buffer
	haploData t;

	// parse the command line options
	std::string region = t.parseCommandLine(argc, argv);

	// check input BAM file for errors
	t.checkBAM();

	// initialize the sample data structure
	t.bam_smpl_init();

	// add samples
	t.bam_smpl_add();

	// initialize error model
	t.em = errmod_init(1.-0.83);

	// parse genomic region
	int k = bam_parse_region(t.h, region, &chr, &beg, &end);
	if (k < 0)
	{
		msg = "Bad genome coordinates: " + region;
		fatal_error(msg, __FILE__, __LINE__, 0);
	}

	// fetch reference sequence
	t.ref_base = faidx_fetch_seq(t.fai_file, t.h->target_name[chr], 0, 0x7fffffff, &(t.len));

	// calculate the number of windows
	if (t.flag & BAM_WINDOW)
		num_windows = ((end-beg)-1)/t.win_size;
	else
	{
		t.win_size = (end-beg);
		num_windows = 1;
	}

	// iterate through all windows along specified genomic region
	for (long cw=0; cw < num_windows; cw++)
	{
		// construct genome coordinate string
		std::string scaffold_name(t.h->target_name[chr]);
		std::ostringstream winc(scaffold_name);
		winc.seekp(0, std::ios::end);
		winc << ":" << beg+(cw*t.win_size)+1 << "-" << ((cw+1)*t.win_size)+(beg-1);
		std::string winCoord = winc.str();

		// initialize number of sites to zero
		t.num_sites = 0;

		// parse the BAM file and check if region is retrieved from the reference
		if (t.flag & BAM_WINDOW)
		{
			k = bam_parse_region(t.h, winCoord, &ref, &(t.beg), &(t.end));
			if (k < 0)
			{
				msg = "Bad window coordinates " + winCoord;
				fatal_error(msg, __FILE__, __LINE__, 0);
			}
		}
		else
		{
			ref = chr;
			t.beg = beg;
			t.end = end;
			if (ref < 0)
			{
				msg = "Bad scaffold name: " + region;
				fatal_error(msg, __FILE__, __LINE__, 0);
			}
		}

		// initialize diverge specific variables
		t.init_haplo();

		// create population assignments
		t.assign_pops();

		// initialize pileup
		buf = bam_plbuf_init(make_haplo, &t);

		// fetch region from bam file
		if ((bam_fetch(t.bam_in->x.bam, t.idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
		{
			msg = "Failed to retrieve region " + region + " due to corrupted BAM index file";
			fatal_error(msg, __FILE__, __LINE__, 0);
		}

		// finalize pileup
		bam_plbuf_push(0, buf);

		// calculate haplotype-based statistics
		t.calc_haplo();

		// print results to stdout
		t.print_haplo(chr);

		// take out the garbage
		t.destroy_haplo();
		bam_plbuf_destroy(buf);
	}
	// end of window interation

	errmod_destroy(t.em);
	samclose(t.bam_in);
	bam_index_destroy(t.idx);
	t.bam_smpl_destroy();
	free(t.ref_base);

	return 0;
}

int make_haplo(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
{
	int i;
	int fq;
	unsigned long long sample_cov;
	unsigned long long *cb = NULL;
	haploData *t = NULL;

	// get control data structure
	t = (haploData*)data;

	// only consider sites located in designated region
	if ((t->beg <= (int)pos) && (t->end > (int)pos))
	{
		// allocate memory pileup data
		try
		{
			cb = new unsigned long long [t->sm->n]();
		}
		catch (std::bad_alloc& ba)
		{
			std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
		}

		// call bases
		t->call_base(n, pl, cb);

		// resolve heterozygous sites
		if (!(t->flag & BAM_HETEROZYGOTE))
			clean_heterozygotes(t->sm->n, cb, (int)t->ref_base[pos], t->min_snpQ);

		// determine if site is segregating
		fq = segbase(t->sm->n, cb, t->ref_base[pos], t->min_snpQ);

		// determine how many samples pass the quality filters
		sample_cov = qfilter(t->sm->n, cb, t->min_rmsQ, t->min_depth, t->max_depth);

		for (i=0; i < t->sm->npops; i++)
			t->pop_sample_mask[i] = sample_cov & t->pop_mask[i];

		if (popcount64(sample_cov) == t->sm->n)
		{
			// calculate the site type
			t->types[t->num_sites] = t->cal_site_type(cb);

			if (fq > 0)
			{
				t->hap.pos[t->segsites] = pos;
				t->hap.ref[t->segsites] = (unsigned char)bam_nt16_table[(int)t->ref_base[pos]];
				for (i=0; i < t->sm->n; i++)
				{
					t->hap.rms[i][t->segsites] = (cb[i]>>(CHAR_BIT*6))&0xffff;
					t->hap.snpq[i][t->segsites] = (cb[i]>>(CHAR_BIT*4))&0xffff;
					t->hap.num_reads[i][t->segsites] = (cb[i]>>(CHAR_BIT*2))&0xffff;
					t->hap.base[i][t->segsites] = bam_nt16_table[(int)iupac[(cb[i]>>CHAR_BIT)&0xff]];
					if (cb[i]&0x2ULL)
						t->hap.seq[i][t->segsites/64] |= 0x1ULL << t->segsites%64;
				}
				t->hap.idx[t->segsites] = t->num_sites;
				t->segsites++;
			}
			t->num_sites++;
		}

		// take out the garbage
		delete [] cb;
	}
	return 0;
}

void haploData::calc_haplo(void)
{
	calc_diff_matrix(*(this));
	haplo_func do_haplo[3] = {&haploData::calc_nhaps, &haploData::calc_ehhs, &haploData::calc_minDxy};
	(this->*do_haplo[output])();
}

void haploData::calc_nhaps(void)
{
	int i, j, k, f;
	int nelem;

	for (i=0; i < sm->npops; i++)
	{
		nelem = pop_nsmpl[i];
		if (nelem > 1)
		{
			std::vector<int> b;
			std::vector<int>::iterator it1;
			std::vector<int>::iterator it2;

			// initialize the haplotype identity vector
			for (j=0; j < sm->n; j++)
				if (CHECK_BIT(pop_mask[i], j))
					b.push_back(j);

			// assign haplotype identifiers
			for (j=0, it1=b.begin(); j < nelem-1; j++, it1++)
			{
				it2 = b.begin();
				for (k=j+1, std::advance(it2, j+1); k < nelem; k++, it2++)
					if ((diff_matrix[j][k] == 0) && (*it2 > *it1))
						b.at(k)=j;
			}

			// count number of haplotypes and calculate haplotype diversity
			int ff = 0;
			double sh;
			for (j=0; j < (int)b.size(); j++)
			{
				if ((f = count(b.begin(), b.end(), j)) > 0)
					++nhaps[i];
				ff += f*f;
			}
			sh = (double)(ff)/(double)SQ(nelem);
			hdiv[i] = 1.0-((1.0-sh)*(double)(nelem/(nelem-1)));
		}
		else
		{
			nhaps[i] = 1;
			hdiv[i] = 1.0;
		}
	}
}

void haploData::calc_ehhs(void)
{
	int i, j;

	calc_nhaps();

	for (i=0; i < sm->npops; i++)
	{
		if (pop_nsmpl[i] < 4)
			ehhs[i] = std::numeric_limits<double>::quiet_NaN();
		else
		{
			unsigned short popf;
			unsigned long long pop_type;
			std::list<unsigned long long> pop_site;

			// make list container of all non-singleton partitions present in population i
			for (j=0; j < segsites; j++)
			{
				pop_type = types[hap.idx[j]] & pop_mask[i];
				popf = popcount64(pop_type);
				if ((popf > 1) && (popf < pop_nsmpl[i]-1))
					pop_site.push_back(pop_type);
			}

			// count unique partitions
			int before;
			int after;
			int part_count = 0;
			int part_max_count = 0;
			unsigned long long part_type = 0;
			unsigned long long part_type_comp = 0;
			unsigned long long max_site = 0;
			double sh;
			std::list<unsigned long long> uniq(pop_site);
			std::list<unsigned long long>::iterator it;
			uniq.sort();
			uniq.unique();

			for (it=uniq.begin(); it != uniq.end(); it++)
			{
				part_type = *it;

				// find the complement of part_type
				for (j=0; j < sm->n; j++)
					if (~CHECK_BIT(part_type, j) && CHECK_BIT(pop_mask[i], j))
						part_type_comp |= 0x1ULL << j;
				before = static_cast<int>(pop_site.size());
				pop_site.remove(part_type);
				pop_site.remove(part_type_comp);
				after = static_cast<int>(pop_site.size());
				part_count = (before-after)+1;
				if (part_count > part_max_count)
				{
					part_max_count = part_count;
					max_site = part_type;
				}
			}

			// calculate site heterozygosity
			popf = popcount64(max_site);
			sh = (1.0-((double)(SQ(popf)+((pop_nsmpl[i]-popf)*(pop_nsmpl[i]-popf)))/SQ(pop_nsmpl[i])))*(double)(pop_nsmpl[i]/(pop_nsmpl[i]-1));

			// calculate site-specific extended haplotype homozygosity
			ehhs[i] = hdiv[i]/(1.0-sh);
		}
	}
}

void haploData::calc_minDxy(void)
{
	int i, j, v, w;
	int npops = sm->npops;
	int n = sm->n;

	for (i=0; i < npops; i++)
	{
		for (j=i; j < npops; j++)
		{
			if (i != j)
				minDxy[i*npops+(j-(i+1))] = std::numeric_limits<unsigned int>::max();
			for (v=0; v < n-1; v++)
			{
				for (w=v+1; w < n; w++)
				{
					if (CHECK_BIT(pop_mask[i], v) && CHECK_BIT(pop_mask[j], w))
					{
						if (i == j)
							piw[i] += (double)diff_matrix[v][w];
						else
						{
							pib[i*npops+(j-(i+1))] += (double)diff_matrix[v][w];
							minDxy[i*npops+(j-(i+1))] = minDxy[i*npops+(j-(i+1))] < diff_matrix[v][w] ? minDxy[i*npops+(j-(i+1))] : diff_matrix[v][w];
						}
					}
				}
			}
			if (i != j)
				pib[i*npops+(j-(i+1))] *= 1.0/(double)(pop_nsmpl[i]*pop_nsmpl[j]);
			else
			{
				piw[i] *= 2.0/(double)(pop_nsmpl[i]*(pop_nsmpl[i]-1));
				if (isnan(piw[i]))
					piw[i] = 0.0;
			}
		}
	}
}

void haploData::print_haplo(int chr)
{
	int i, j;

	//print coordinate information and number of aligned sites
	std::cout << h->target_name[chr] << "\t" << beg+1 << "\t" << end+1 << "\t" << num_sites;

	switch(output)
	{
		case 0:
			for (i=0; i < sm->npops; i++)
			{
				if (num_sites >= min_sites)
				{
					std::cout << "\tK[" << sm->popul[i] << "]:\t" << nhaps[i];
					std::cout << "\tKdiv[" << sm->popul[i] << "]:";
					std::cout << "\t" << std::fixed << std::setprecision(5) << 1.0 - hdiv[i];
				}
				else
				{
					std::cout << "\tK[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
					std::cout << "\tKdiv[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
				}
			}
			break;
		case 1:
			for (i=0; i < sm->npops; i++)
			{
				if (num_sites >= min_sites)
				{
					if (isnan(ehhs[i]))
						std::cout << "\tEHHS[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
					else
					{
						std::cout << "\tEHHS[" << sm->popul[i] << "]:";
						std::cout << "\t" << std::fixed << std::setprecision(5) << ehhs[i];
					}
				}
				else
					std::cout << "\tEHHS[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
			}
			break;
		case 2:
			for (i=0; i < sm->npops; i++)
			{
				if (num_sites >= min_sites)
				{
					std::cout << "\tpi[" << sm->popul[i] << "]:";
					std::cout << "\t" << std::fixed << std::setprecision(5) << piw[i];
				}
				else
					std::cout << "\tpi[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
			}

			for (i=0; i < sm->npops-1; i++)
			{
				for (j=i+1; j < sm->npops; j++)
				{
					if (num_sites >= min_sites)
					{
						std::cout << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:";
						std::cout << "\t" << std::fixed << std::setprecision(5) << pib[i*sm->npops+(j-(i+1))];
						std::cout << "\tmin[" << sm->popul[i] << "-" << sm->popul[j] << "]:";
						std::cout << "\t" << minDxy[i*sm->npops+(j-(i+1))];
					}
					else
					{
						std::cout << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:\t" << std::setw(7) << "NA";
						std::cout << "\tmin[" << sm->popul[i] << "-" << sm->popul[j] << "]:\t" << std::setw(7) << "NA";
					}
				}
			}
			break;
		default:
			break;
	}
	std::cout << std::endl;
}

void calc_diff_matrix(haploData &h)
{
	int i, j, k;
	int n = h.sm->n;
	int segs = h.segsites;

	// calculate number of pairwise differences
	for (i=0; i < n-1; i++)
		for (j=i+1; j < n; j++)
		{
			for (k=0; k <= SEG_IDX(segs); k++)
				h.diff_matrix[j][i] += hamming_distance(h.hap.seq[i][k], h.hap.seq[j][k]);
			h.diff_matrix[i][j] = h.diff_matrix[j][i];
		}
}

std::string haploData::parseCommandLine(int argc, char *argv[])
{
#ifdef _MSC_VER
	struct _stat finfo;
#else
	struct stat finfo;
#endif
	std::vector<std::string> glob_opts;
	std::string msg;

	GetOpt::GetOpt_pp args(argc, argv);
	args >> GetOpt::Option('f', reffile);
	args >> GetOpt::Option('h', headfile);
	args >> GetOpt::Option('o', output);
	args >> GetOpt::Option('m', min_depth);
	args >> GetOpt::Option('x', max_depth);
	args >> GetOpt::Option('q', min_rmsQ);
	args >> GetOpt::Option('s', min_snpQ);
	args >> GetOpt::Option('a', min_mapQ);
	args >> GetOpt::Option('b', min_baseQ);
	args >> GetOpt::Option('k', min_sites);
	args >> GetOpt::Option('w', win_size);
	if (args >> GetOpt::OptionPresent('w'))
	{
		win_size *= 1000;
		flag |= BAM_WINDOW;
	}
	if (args >> GetOpt::OptionPresent('h'))
		flag |= BAM_HEADERIN;
	if (args >> GetOpt::OptionPresent('i'))
		flag |= BAM_ILLUMINA;
	args >> GetOpt::GlobalOption(glob_opts);

	// run some checks on the command line

	// check if output option is valid
	if ((output < 0) || (output > 2))
	{
		msg = "Not a valid output option";
		fatal_error(msg, __FILE__, __LINE__, &haploUsage);
	}

	// if no input BAM file is specified -- print usage and exit
	if (glob_opts.size() < 2)
	{
		msg = "Need to specify BAM file name";
		fatal_error(msg, __FILE__, __LINE__, &haploUsage);
	}
	else
		bamfile = glob_opts[0];

	// check if specified BAM file exists on disk
	if ((stat(bamfile.c_str(), &finfo)) != 0)
	{
		msg = "Specified input file: " + bamfile + " does not exist";
		switch(errno)
		{
			case ENOENT:
				std::cerr << "File not found" << std::endl;
				break;
			case EINVAL:
				std::cerr << "Invalid parameter to stat" << std::endl;
				break;
			default:
				std::cerr << "Unexpected error in stat" << std::endl;
				break;
		}
		fatal_error(msg, __FILE__, __LINE__, 0);
	}

	// check if fastA reference file is specified
	if (reffile.empty())
	{
		msg = "Need to specify fastA reference file";
		fatal_error(msg, __FILE__, __LINE__, &haploUsage);
	}

	// check is fastA reference file exists on disk
	if ((stat(reffile.c_str(), &finfo)) != 0)
	{
		switch(errno)
		{
			case ENOENT:
				std::cerr << "File not found" << std::endl;
				break;
			case EINVAL:
				std::cerr << "Invalid parameter to stat" << std::endl;
				break;
			default:
				std::cerr << "Unexpected error in stat" << std::endl;
				break;
		}
		msg = "Specified reference file: " + reffile + " does not exist";
		fatal_error(msg, __FILE__, __LINE__, 0);
	}

	//check if BAM header input file exists on disk
	if (flag & BAM_HEADERIN)
	{
		if ((stat(headfile.c_str(), &finfo)) != 0)
		{
			switch(errno)
			{
				case ENOENT:
					std::cerr << "File not found" << std::endl;
					break;
				case EINVAL:
					std::cerr << "Invalid parameter to stat" << std::endl;
					break;
				default:
					std::cerr << "Unexpected error in stat" << std::endl;
					break;
			}
			msg = "Specified header file: " + headfile + " does not exist";
			fatal_error(msg, __FILE__, __LINE__, 0);
		}
	}

	// return the index of first non-optioned argument
	return glob_opts[1];
}

haploData::haploData(void)
{
	derived_type = HAPLO;
	output = 0;
	min_sites = 10;
	win_size = 0;
}

void haploData::init_haplo(void)
{
	int i;
	int length = end-beg;
	int n = sm->n;
	int npops = sm->npops;

	segsites = 0;

	try
	{
		types = new unsigned long long [length]();
		pop_mask = new unsigned long long [npops]();
		pop_nsmpl = new unsigned char [npops]();
		pop_sample_mask = new unsigned long long [npops]();
		nhaps = new int [npops]();
		hdiv = new double [npops]();
		piw = new double [npops]();
		pib = new double [npops*(npops-1)]();
		ehhs = new double [npops]();
		minDxy = new unsigned short [npops*(npops-1)]();
		hap.pos = new unsigned int [length]();
		hap.idx = new unsigned int [length]();
		hap.ref = new unsigned char [length]();
		hap.seq = new unsigned long long* [n];
		hap.base = new unsigned char* [n];
		hap.rms = new unsigned short* [n];
		hap.snpq = new unsigned short* [n];
		hap.num_reads = new unsigned short* [n];
		diff_matrix = new unsigned short* [n];
		for (i=0; i < n; i++)
		{
			hap.seq[i] = new unsigned long long [length]();
			hap.base[i] = new unsigned char [length]();
			hap.rms[i] = new unsigned short [length]();
			hap.snpq[i] = new unsigned short [length]();
			hap.num_reads[i] = new unsigned short [length]();
			diff_matrix[i] = new unsigned short [n]();
		}
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
	}
}

void haploData::destroy_haplo(void)
{
	int i;

	delete [] pop_mask;
	delete [] types;
	delete [] pop_nsmpl;
	delete [] pop_sample_mask;
	delete [] hdiv;
	delete [] nhaps;
	delete [] piw;
	delete [] pib;
	delete [] ehhs;
	delete [] minDxy;
	delete [] hap.pos;
	delete [] hap.idx;
	delete [] hap.ref;
	for (i=0; i < sm->n; i++)
	{
		delete [] hap.seq[i];
		delete [] hap.base[i];
		delete [] hap.num_reads[i];
		delete [] hap.snpq[i];
		delete [] hap.rms[i];
		delete [] diff_matrix[i];
	}
	delete [] hap.seq;
	delete [] hap.base;
	delete [] hap.snpq;
	delete [] hap.rms;
	delete [] hap.num_reads;
	delete [] diff_matrix;
}

void haploData::haploUsage(void)
{
	std::cerr << std::endl;
	std::cerr << "Usage:   popbam haplo [options] <in.bam> [region]" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Options: -i          base qualities are Illumina 1.3+     [ default: Sanger ]" << std::endl;
	std::cerr << "         -h  FILE    Input header file                    [ default: none ]" << std::endl;
	std::cerr << "         -w  INT     use sliding window of size (kb)" << std::endl;
	std::cerr << "         -k  INT     minimum number of sites in window    [ default: 10 ]" << std::endl;
	std::cerr << "         -o  INT     analysis to output                   [ default: 0 ]" << std::endl;
	std::cerr << "                     0 : number of haplotypes" << std::endl;
	std::cerr << "                     1 : extended haplotype homozygosity statistic" << std::endl;
	std::cerr << "                     2 : minimum Dxy statistic" << std::endl;
	std::cerr << "         -f  FILE    Reference fastA file" << std::endl;
	std::cerr << "         -m  INT     minimum read coverage                [ default: 3 ]" << std::endl;
	std::cerr << "         -x  INT     maximum read coverage                [ default: 255 ]" << std::endl;
	std::cerr << "         -q  INT     minimum rms mapping quality          [ default: 25 ]" << std::endl;
	std::cerr << "         -s  INT     minimum snp quality                  [ default: 25 ]" << std::endl;
	std::cerr << "         -a  INT     minimum map quality                  [ default: 13 ]" << std::endl;
	std::cerr << "         -b  INT     minimum base quality                 [ default: 13 ]" << std::endl;
	std::cerr << std::endl;
	exit(EXIT_FAILURE);
}
