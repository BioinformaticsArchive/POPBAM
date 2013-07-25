/** \file pop_nucdiv.cpp
 *  \brief Functions for calculating nucleotide diversity statistics
 *  \author Daniel Garrigan
 *  \version 0.3
*/

#include "pop_nucdiv.h"
#include "tables.h"

int main_nucdiv(int argc, char *argv[])
{
    int chr;                  //! chromosome identifier
    int beg;                  //! beginning coordinate for analysis
    int end;                  //! end coordinate for analysis
    int ref;                  //! ref
    long num_windows;         //! number of windows
    std::string msg;          //! string for error message
    bam_plbuf_t *buf;         //! pileup buffer
    nucdivData t;

    // parse the command line options
    std::string region = t.parseCommandLine(argc, argv);

    // check input BAM file for errors
    t.checkBAM();

    // initialize the sample data structure
    t.bam_smpl_init();

    // add samples
    t.bam_smpl_add();

    // initialize error model
    t.em = errmod_init(1.0-0.83);

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

        // initialize nucdiv variables
        t.init_nucdiv();

        // create population assignments
        t.assign_pops();

        // set default minimum sample size as
        // the number of samples in the population
        t.set_min_pop_n();

        // initialize pileup
        buf = bam_plbuf_init(make_nucdiv, &t);

        // fetch region from bam file
        if ((bam_fetch(t.bam_in->x.bam, t.idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
        {
            msg = "Failed to retrieve region " + region + " due to corrupted BAM index file";
            fatal_error(msg, __FILE__, __LINE__, 0);
        }

        // finalize pileup
        bam_plbuf_push(0, buf);

        // calculate nucleotide diversity in window
        calc_diff_matrix(t);
        t.calc_nucdiv();

        // print results to stdout
        t.print_nucdiv(chr);

        // take out the garbage
        t.destroy_nucdiv();
        bam_plbuf_destroy(buf);
    }
    // end of window iteration

    errmod_destroy(t.em);
    samclose(t.bam_in);
    bam_index_destroy(t.idx);
    t.bam_smpl_destroy();
    free(t.ref_base);

    return 0;
}

int make_nucdiv(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
{
    int i;
    int fq;
    unsigned long long sample_cov;
    unsigned long long *cb = NULL;
    nucdivData *t = NULL;

    // get control data structure
    t = (nucdivData*)data;

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

        if (popcount64(sample_cov) == t->sm->n)
        {
            // calculate the site type
            t->types[t->num_sites] = t->cal_site_type(cb);

            if (fq > 0)
            {
                t->hap.pos[t->segsites] = pos;
                t->hap.ref[t->segsites] = bam_nt16_table[(int)t->ref_base[pos]];
                for (i=0; i < t->sm->n; i++)
                {
                    t->hap.rms[i][t->segsites] = (cb[i]>>(CHAR_BIT*6))&0xffff;
                    t->hap.snpq[i][t->segsites] = (cb[i]>>(CHAR_BIT*4))&0xffff;
                    t->hap.num_reads[i][t->segsites] = (cb[i]>>(CHAR_BIT*2))&0xffff;
                    t->hap.base[i][t->segsites] = bam_nt16_table[(int)iupac[(cb[i]>>CHAR_BIT)&0xff]];
                    if (cb[i]&0x2ULL)
                        t->hap.seq[i][t->segsites/64] |= 0x1ULL << t->segsites%64;
                }
                for (i=0; i < t->sm->npops; i++)
                    t->pop_sample_mask[i][t->segsites] = sample_cov & t->pop_mask[i];
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

void nucdivData::calc_nucdiv(void)
{
    int i, j, v, w;
    const int npops = sm->npops;
    const int n = sm->n;

    for (i=0; i < npops; i++)
    {
        for (j=i; j < npops; j++)
        {
            // calculate pairwise differences
            for (v=0; v < n-1; v++)
                for (w=v+1; w < n; w++)
                {
                    if (CHECK_BIT(pop_mask[i], v) && CHECK_BIT(pop_mask[j], w))
                    {
                        if (i == j)
                            piw[i] += (double)diff_matrix[v][w];
                        else
                            pib[i*npops+(j-(i+1))] += (double)diff_matrix[v][w];
                    }
                }

            if (i != j)
                pib[i*npops+(j-(i+1))] *= 1.0/(double)(pop_nsmpl[i]*pop_nsmpl[j]);
            else
            {
                piw[i] *= 2.0/(double)(pop_nsmpl[i]*(pop_nsmpl[i]-1));
                if (isnan(piw[i]))
                    piw[i] = 0.;
            }
        }
    }
}

// overloaded calc_diff_matrix functions
void calc_diff_matrix(nucdivData &h)
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

void nucdivData::print_nucdiv(int chr)
{
    int i, j;

    std::cout << h->target_name[chr] << "\t" << beg+1 << "\t" << end+1 << "\t" << num_sites;

    for (i=0; i < sm->npops; i++)
    {
        if (num_sites >= min_sites)
        {
            std::cout << "\tpi[" << sm->popul[i] << "]:";
            std::cout << "\t" << std::fixed << std::setprecision(5) << piw[i]/num_sites;
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
                std::cout << "\t" << std::fixed << std::setprecision(5) << pib[i*sm->npops+(j-(i+1))]/num_sites;
            }
            else
                std::cout << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:\t" << std::setw(7) << "NA";
        }
    }
    std::cout << std::endl;
}

void nucdivData::set_min_pop_n(void)
{
    for (int j=0; j < sm->npops; j++)
        min_pop_n[j] = (unsigned short)pop_nsmpl[j];
}

std::string nucdivData::parseCommandLine(int argc, char *argv[])
{
    std::vector<std::string> glob_opts;
    std::string msg;
#ifdef _MSC_VER
    struct _stat finfo;
#else
    struct stat finfo;
#endif

    //read command line options
    GetOpt::GetOpt_pp args(argc, argv);
    args >> GetOpt::Option('f', reffile);
    args >> GetOpt::Option('h', headfile);
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
    if (args >> GetOpt::OptionPresent('p'))
        flag |= BAM_OUTGROUP;
    if (args >> GetOpt::OptionPresent('i'))
        flag |= BAM_ILLUMINA;
    if (args >> GetOpt::OptionPresent('n'))
        flag |= BAM_MINPOPSAMPLE;
    args >> GetOpt::GlobalOption(glob_opts);

    // run some checks on the command line

    // if no input BAM file is specified -- print usage and exit
    if (glob_opts.size() < 2)
    {
        msg = "Need to specify BAM file name";
        fatal_error(msg, __FILE__, __LINE__, &nucdivUsage);
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
        msg = "Need to specify fasta reference file";
        fatal_error(msg, __FILE__, __LINE__, &nucdivUsage);
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

nucdivData::nucdivData(void)
{
    derived_type = NUCDIV;
    min_sites = 10;
    win_size = 0;
}

void nucdivData::init_nucdiv(void)
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
        pop_sample_mask = new unsigned long long* [npops];
        piw = new double [npops]();
        pib = new double [npops*(npops-1)]();
        min_pop_n = new unsigned short [npops]();
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
        for (i=0; i < npops; i++)
            pop_sample_mask[i] = new unsigned long long [length]();
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
    }
}

void nucdivData::destroy_nucdiv(void)
{
    int i;
    int npops = sm->npops;

    delete [] pop_mask;
    delete [] types;
    delete [] pop_nsmpl;
    delete [] min_pop_n;
    delete [] piw;
    delete [] pib;
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
    for (i=0; i < npops; i++)
        delete [] pop_sample_mask[i];
    delete [] pop_sample_mask;
    delete [] hap.seq;
    delete [] hap.base;
    delete [] hap.snpq;
    delete [] hap.rms;
    delete [] hap.num_reads;
    delete [] diff_matrix;
}

void nucdivData::nucdivUsage(void)
{
    std::cerr << std::endl;
    std::cerr << "Usage:   popbam nucdiv [options] <in.bam> [region]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options: -i          base qualities are Illumina 1.3+     [ default: Sanger ]" << std::endl;
    std::cerr << "         -h  FILE    Input header file                    [ default: none ]" << std::endl;
    std::cerr << "         -w  INT     use sliding window of size (kb)" << std::endl;
    std::cerr << "         -k  INT     minimum number of sites in window    [ default: 10 ]" << std::endl;
    std::cerr << "         -n  INT     minimum sample size per population   [ default: all samples present ]" << std::endl;
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
