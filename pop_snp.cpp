/** \file pop_snp.cpp
 *  \brief Functions for extracting SNP calls from BAM files
 *  \author Daniel Garrigan
 *  \version 0.3
*/
#include "pop_snp.h"
#include "tables.h"

int main_snp(int argc, char *argv[])
{
    bool found = false;       //! is the outgroup found?
    int chr;                  //! chromosome identifier
    int beg;                  //! beginning coordinate for analysis
    int end;                  //! end coordinate for analysis
    int ref;                  //! reference allele
    long num_windows;         //! total number of windows
    std::string msg;          //! string for error message
    bam_plbuf_t *buf;         //! pileup buffer
    snpData t;                //! data object for the snp function

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

    // if outgroup option is used check to make sure it exists
    if (t.flag & BAM_OUTGROUP)
    {
        for (int i=0; i < t.sm->n; i++)
            if (strcmp(t.sm->smpl[i], t.outgroup.c_str()) == 0)
            {
                t.outidx = i;
                found = true;
            }
        if (!found)
        {
            msg = "Specified outgroup " + t.outgroup + " not found";
            fatal_error(msg, __FILE__, __LINE__, 0);
        }
    }

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
        t.init_snp();

        // create population assignments
        t.assign_pops();

        // print ms header if first window iteration
        if ((t.output == 2) && (cw == 0))
            t.print_ms_header(num_windows);

        // initialize pileup
        buf = bam_plbuf_init(make_snp, &t);

        // fetch region from bam file
        if ((bam_fetch(t.bam_in->x.bam, t.idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
        {
            msg = "Failed to retrieve region " + region + " due to corrupted BAM index file";
            fatal_error(msg, __FILE__, __LINE__, 0);
        }

        // finalize pileup
        bam_plbuf_push(0, buf);

        // print results to stdout
        t.print_snp(chr);

        // take out the garbage
        t.destroy_snp();
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

int make_snp(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
{
    int i;
    int fq;
    unsigned long long sample_cov;
    unsigned long long *cb = NULL;
    snpData *t = NULL;

    // get control data structure
    t = (snpData*)data;

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

void snpData::print_snp(int chr)
{
    snp_func fp[3] = {&snpData::print_popbam_snp, &snpData::print_sweep, &snpData::print_ms};
    (this->*fp[output])(chr);
}

void snpData::print_popbam_snp(int chr)
{
    int i, j;

    for (i=0; i < segsites; i++)
    {
        std::cout << h->target_name[chr] << "\t" << hap.pos[i]+1 << "\t";
        std::cout << bam_nt16_rev_table[hap.ref[i]];
        for (j=0; j < sm->n; j++)
        {
            std::cout << "\t" << bam_nt16_rev_table[hap.base[j][i]];
            std::cout << "\t" << hap.snpq[j][i];
            std::cout << "\t" << hap.rms[j][i];
            std::cout << "\t" << hap.num_reads[j][i];
        }
        std::cout << std::endl;
    }
}

void snpData::print_sweep(int chr)
{
    int i, j;
    unsigned short freq;
    unsigned short pop_n;
    unsigned long long pop_type;

    for (i=0; i < segsites; i++)
    {
        std::cout << h->target_name[chr] << "\t" << hap.pos[i]+1;
        for (j=0; j < sm->npops; j++)
        {
            // population-specific site type
            pop_type = types[hap.idx[i]] & pop_sample_mask[j][i];

            // assign derived allele counts and sample sizes
            pop_n = popcount64(pop_sample_mask[j][i]);
            if ((flag & BAM_OUTGROUP) && CHECK_BIT(types[hap.idx[i]], outidx))
                freq = pop_n - popcount64(pop_type);
            else
                freq = popcount64(pop_type);
            std::cout << "\t" << freq << "\t" << pop_n;
        }
        std::cout << std::endl;
    }
}

void snpData::print_ms(int chr)
{
    int i, j;

    std::cout << "//" << std::endl;
    std::cout << "segsites: " << segsites << std::endl;
    std::cout << "positions: ";
    for (i=0; i < segsites; i++)
        std::cout << std::setprecision(8) << (double)(hap.pos[i]-beg)/(end-beg) << " ";
    std::cout << std::endl;
    for (i=0; i < sm->n; i++)
    {
        for (j=0; j < segsites; j++)
        {
            if ((flag & BAM_OUTGROUP) && CHECK_BIT(types[hap.idx[j]], outidx))
            {
                if (CHECK_BIT(hap.seq[i][j/64], j%64))
                    std::cout << "0";
                else
                    std::cout << "1";
            }
            else
            {
                if (CHECK_BIT(hap.seq[i][j/64], j%64))
                    std::cout << "1";
                else
                    std::cout << "0";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

}

void snpData::print_ms_header(long nwindows)
{
    if (sm->npops > 1)
    {
        std::cout << "ms " << sm->n << " " << nwindows << " -t 5.0 -I " << sm->npops << " ";
        for (int i=0; i < sm->npops; i++)
            std::cout << (int)pop_nsmpl[i] << " ";
    }
    else
        std::cout << "ms " << sm->n << " " << nwindows << " -t 5.0 ";

    std::cout << std::endl << "1350154902" << std::endl << std::endl;
}

std::string snpData::parseCommandLine(int argc, char *argv[])
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
    args >> GetOpt::Option('m', min_depth);
    args >> GetOpt::Option('x', max_depth);
    args >> GetOpt::Option('q', min_rmsQ);
    args >> GetOpt::Option('s', min_snpQ);
    args >> GetOpt::Option('a', min_mapQ);
    args >> GetOpt::Option('b', min_baseQ);
    args >> GetOpt::Option('o', output);
    args >> GetOpt::Option('z', het_prior);
    args >> GetOpt::Option('p', outgroup);
    args >> GetOpt::Option('w', win_size);
    if (args >> GetOpt::OptionPresent('w'))
    {
        win_size *= 1000;
        flag |= BAM_WINDOW;
    }
    if (args >> GetOpt::OptionPresent('h'))
        flag |= BAM_HEADERIN;
    if (args >> GetOpt::OptionPresent('v'))
        flag |= BAM_VARIANT;
    if (args >> GetOpt::OptionPresent('i'))
        flag |= BAM_ILLUMINA;
    if (args >> GetOpt::OptionPresent('z'))
        flag |= BAM_HETEROZYGOTE;
    if (args >> GetOpt::OptionPresent('p'))
        flag |= BAM_OUTGROUP;
    args >> GetOpt::GlobalOption(glob_opts);

    // run some checks on the command line

    // check if output option is valid
    if ((output < 0) || (output > 2))
    {
        msg = "Not a valid output option";
        fatal_error(msg, __FILE__, __LINE__, &snpUsage);
    }

    // if no input BAM file is specified -- print usage and exit
    if (glob_opts.size() < 2)
    {
        msg = "Need to specify input BAM file name";
        fatal_error(msg, __FILE__, __LINE__, &snpUsage);
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
        msg = "Need to specify a fasta reference file";
        fatal_error(msg, __FILE__, __LINE__, &snpUsage);
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

snpData::snpData(void)
{
    derived_type = SNP;
    output = 0;
    outidx = 0;
    win_size = 0;
}

void snpData::init_snp(void)
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
        hap.pos = new unsigned int [length]();
        hap.idx = new unsigned int [length]();
        hap.ref = new unsigned char [length]();
        hap.seq = new unsigned long long* [n];
        hap.base = new unsigned char* [n];
        hap.rms = new unsigned short* [n];
        hap.snpq = new unsigned short* [n];
        hap.num_reads = new unsigned short* [n];
        for (i=0; i < n; i++)
        {
            hap.seq[i] = new unsigned long long [length]();
            hap.base[i] = new unsigned char [length]();
            hap.rms[i] = new unsigned short [length]();
            hap.snpq[i] = new unsigned short [length]();
            hap.num_reads[i] = new unsigned short [length]();
        }
        for (i=0; i < npops; i++)
            pop_sample_mask[i] = new unsigned long long [length]();
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
    }
}

void snpData::destroy_snp(void)
{
    int i;
    int npops = sm->npops;

    delete [] pop_mask;
    delete [] pop_nsmpl;
    delete [] types;
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
    }
    for (i=0; i < npops; i++)
        delete [] pop_sample_mask[i];
    delete [] pop_sample_mask;
    delete [] hap.seq;
    delete [] hap.base;
    delete [] hap.snpq;
    delete [] hap.rms;
    delete [] hap.num_reads;
}

void snpData::snpUsage(void)
{
    std::cerr << std::endl;
    std::cerr << "Usage:   popbam snp [options] <in.bam> [region]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options: -i          base qualities are Illumina 1.3+     [ default: Sanger ]" << std::endl;
    std::cerr << "         -h  FILE    Input header file                    [ default: none ]" << std::endl;
    std::cerr << "         -v          output variant sites only            [ default: All sites ]" << std::endl;
    std::cerr << "         -z  FLT     output heterozygous base calls       [ default: Consensus ]" << std::endl;
    std::cerr << "         -w  INT     use sliding window of size (kb)" << std::endl;
    std::cerr << "         -p  STR     sample name of outgroup              [ default: reference ]" << std::endl;
    std::cerr << "         -o  INT     output format                        [ default: 0 ]" << std::endl;
    std::cerr << "                     0 : popbam snp format" << std::endl;
    std::cerr << "                     1 : SweepFinder snp format" << std::endl;
    std::cerr << "                     2 : MS format" << std::endl;
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
