/** \file pop_tree.cpp
 *  \brief Functions for constructing phylogenetic trees file from BAM files
 *  \author Daniel Garrigan
 *  \version 0.3
*/

#include "pop_tree.h"
#include "tables.h"

int main_tree(int argc, char *argv[])
{
    int chr;                  //! chromosome identifier
    int beg;                  //! beginning coordinate for analysis
    int end;                  //! end coordinate for analysis
    int ref;                  //! ref
    long num_windows;         //! number of windows
    std::string msg;          //! string for error message
    bam_plbuf_t *buf;         //! pileup buffer
    treeData t;

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

    // extract name of reference sequence
    t.refid = get_refid(t.h->text);

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

        // initialize tree-specific variables
        t.init_tree();

        // create population assignments
        t.assign_pops();

        // initialize pileup
        buf = bam_plbuf_init(make_tree, &t);

        // fetch region from bam file
        if ((bam_fetch(t.bam_in->x.bam, t.idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
        {
            msg = "Failed to retrieve region " + region + " due to corrupted BAM index file";
            fatal_error(msg, __FILE__, __LINE__, 0);
        }

        // finalize pileup
        bam_plbuf_push(0, buf);

        // count pairwise differences
        calc_diff_matrix(t);

        // construct distance matrix
        t.calc_dist_matrix();

        // construct nj tree
        t.make_nj(chr);

        // take out the garbage
        t.destroy_tree();
        bam_plbuf_destroy(buf);
    }
    // end of window interation

    errmod_destroy(t.em);
    samclose(t.bam_in);
    bam_index_destroy(t.idx);
    t.bam_smpl_destroy();
    delete [] t.refid;
    free(t.ref_base);

    return 0;
}

int make_tree(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
{
    int i;
    int fq;
    unsigned long long sample_cov;
    unsigned long long *cb = NULL;
    treeData *t = NULL;

    // get control data structure
    t = (treeData*)data;

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

void treeData::make_nj(int chr)
{
    int i;
    tree curtree;
    node **cluster = NULL;

    if ((num_sites < min_sites) || (segsites < 1))
    {
        std::cout << h->target_name[chr] << "\t" << beg+1 << "\t" << end+1 << "\t" << num_sites;
        std::cout << "\tNA" << std::endl;
        return;
    }

    try
    {
        cluster = new node* [ntaxa];
        enterorder = new int [ntaxa];
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
    }

    for (i=0; i < ntaxa; i++)
        enterorder[i] = i+1;

    tree_init(&curtree.nodep);
    node *p;
    p = curtree.nodep[2*ntaxa-2]->next;
    curtree.nodep[2*ntaxa-2]->next = curtree.nodep[2*ntaxa-2];
    free(p->next);
    free(p);
    setup_tree(&curtree);

    for (i=0; i < ntaxa; i++)
        cluster[i] = curtree.nodep[i];

    join_tree(curtree, cluster);
    curtree.start = curtree.nodep[0]->back;
    std::cout << h->target_name[chr] << "\t" << beg+1 << "\t" << end+1 << "\t" << num_sites << "\t";
    print_tree(curtree.start, curtree.start);
    free_tree(&curtree.nodep);
    delete [] cluster;
    delete [] enterorder;
}

void treeData::join_tree(tree curtree, node **cluster)
{
    int nc, nextnode, mini = 0, minj = 0, i, j, jj, ii, ia, ja, nude, iter;
    int el[3];
    int *oc = NULL;
    double fotu2, total = 0, tmin, dio, djo, bi, bj, bk, dmin = 0, da;
    double *av = NULL;
    double **x = NULL;
    double *R = NULL;

    try
    {
        av = new double [ntaxa];
        oc = new int [ntaxa];
        R = new double [ntaxa];
        x = new double* [ntaxa];
        for (i=0; i < ntaxa; i++)
        {
            x[i] = new double [ntaxa];
            memcpy(x[i], dist_matrix[i], ntaxa*sizeof(double));
        }
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
    }

    for (i=0; i < ntaxa-1; i++)
    {
        for (j=i+1; j < ntaxa; j++)
        {
            da = (x[i][j]+x[j][i])/2.0;
            x[i][j] = da;
            x[j][i] = da;
        }
    }

    // first initialization
    fotu2 = ntaxa-2.0;
    nextnode = ntaxa+1;
    for (i=0; i < ntaxa; i++)
    {
        av[i] = 0.0;
        oc[i] = 1;
    }

    // enter main cycle
    iter = ntaxa-3;
    for (nc=1; nc <= iter; nc++)
    {
        for (j=2; j <= ntaxa; j++)
            for (i=0; i <= j-2; i++)
                x[j-1][i] = x[i][j-1];

        tmin = DBL_MAX;

        // compute sij and minimize
        for (i=0; i < ntaxa; i++)
            R[i] = 0.0;
        for (ja=2; ja <= ntaxa; ja++)
        {
            jj = enterorder[ja-1];
            if (cluster[jj-1] != 0)
            {
                for (ia=0; ia <= ja-2; ia++)
                {
                    ii = enterorder[ia];
                    if (cluster[ii-1] != 0)
                    {
                        R[ii-1] += x[ii-1][jj-1];
                        R[jj-1] += x[ii-1][jj-1];
                    }
                }
            }
        }

        for (ja=2; ja <= ntaxa; ja++)
        {
            jj = enterorder[ja-1];
            if (cluster[jj-1] != 0)
            {
                for (ia=0; ia <= ja-2; ia++)
                {
                    ii = enterorder[ia];
                    if (cluster[ii-1] != 0)
                        total = fotu2*x[ii-1][jj-1]-R[ii-1]-R[jj-1];
                    if (total < tmin)
                    {
                        tmin = total;
                        mini = ii;
                        minj = jj;
                    }
                }
            }
        }
        dio = 0.0;
        djo = 0.0;
        for (i=0; i < ntaxa; i++)
        {
            dio += x[i][mini-1];
            djo += x[i][minj-1];
        }
        dmin = x[mini-1][minj-1];
        dio = (dio-dmin)/fotu2;
        djo = (djo-dmin)/fotu2;
        bi = (dmin+dio-djo)*0.5;
        bj = dmin-bi;
        bi -= av[mini-1];
        bj -= av[minj-1];
        hookup(curtree.nodep[nextnode-1]->next, cluster[mini-1]);
        hookup(curtree.nodep[nextnode-1]->next->next, cluster[minj-1]);
        cluster[mini-1]->v = bi;
        cluster[minj-1]->v = bj;
        cluster[mini-1]->back->v = bi;
        cluster[minj-1]->back->v = bj;
        cluster[mini-1] = curtree.nodep[nextnode-1];
        cluster[minj-1] = 0;
        nextnode++;
        av[mini-1] = dmin*0.5;

        // re-initialization
        fotu2 -= 1.0;
        for (j=0; j < ntaxa; j++)
        {
            if (cluster[j] != 0)
            {
                da = (x[mini-1][j]+x[minj-1][j])*0.5;
                if (mini-j-1 < 0)
                    x[mini-1][j] = da;
                if (mini-j-1 > 0)
                    x[j][mini-1] = da;
            }
        }
        for (j=0; j < ntaxa; j++)
        {
            x[minj-1][j] = 0.0;
            x[j][minj-1] = 0.0;
        }
        oc[mini-1] += oc[minj-1];
    }

    // the last cycle
    nude = 1;
    for (i=1; i <= ntaxa; i++)
    {
        if (cluster[i-1] != 0)
        {
            el[nude-1] = i;
            nude++;
        }
    }
    bi = (x[el[0]-1][el[1]-1]+x[el[0]-1][el[2]-1]-x[el[1]-1][el[2]-1])*0.5;
    bj = x[el[0]-1][el[1]-1]-bi;
    bk = x[el[0]-1][el[2]-1]-bi;
    bi -= av[el[0]-1];
    bj -= av[el[1]-1];
    bk -= av[el[2]-1];
    hookup(curtree.nodep[nextnode-1], cluster[el[0]-1]);
    hookup(curtree.nodep[nextnode-1]->next, cluster[el[1]-1]);
    hookup(curtree.nodep[nextnode-1]->next->next, cluster[el[2]-1]);
    cluster[el[0]-1]->v = bi;
    cluster[el[1]-1]->v = bj;
    cluster[el[2]-1]->v = bk;
    cluster[el[0]-1]->back->v = bi;
    cluster[el[1]-1]->back->v = bj;
    cluster[el[2]-1]->back->v = bk;
    curtree.start = cluster[el[0]-1]->back;

    // take out the garbage
    delete [] av;
    delete [] oc;
    delete [] R;
    for (i=0; i < ntaxa; i++)
        delete [] x[i];
    delete [] x;
}

void treeData::hookup(node *p, node *q)
{
    assert(p != 0);
    assert(q != 0);
    p->back = q;
    q->back = p;
}

void treeData::print_tree(node *p, node *start)
{
    if (p->tip)
    {
        if (p->index == 1)
            std::cout << refid;
        else
            std::cout << sm->smpl[p->index-2];
    }
    else
    {
        std::cout << "(";
        print_tree(p->next->back, start);
        std::cout << ",";
        print_tree(p->next->next->back, start);
        if (p == start)
        {
            std::cout << ",";
            print_tree(p->back, start);
        }
        std::cout << ")";
    }
    if (p == start)
        std::cout << ";" << std::endl;
    else
    {
        if (p->v < 0)
            std::cout << ":0.00000";
        else
            std::cout << ":" << std::fixed << std::setprecision(5) << p->v;
    }
}

void calc_diff_matrix(treeData &h)
{
    int i, j, k;
    int n = h.sm->n;
    int segs = h.segsites;

    // calculate number of differences with reference sequence
    for (i=0; i < n; i++)
    {
        for (k=0; k <= SEG_IDX(segs); k++)
            h.diff_matrix[i+1][0] += popcount64(h.hap.seq[i][k]);
        h.diff_matrix[0][i+1] = h.diff_matrix[i+1][0];
    }

    // calculate number of pairwise differences
    for (i=0; i < n-1; i++)
        for (j=i+1; j < n; j++)
        {
            for (k=0; k <= SEG_IDX(segs); k++)
                h.diff_matrix[j+1][i+1] += hamming_distance(h.hap.seq[i][k], h.hap.seq[j][k]);
            h.diff_matrix[i+1][j+1] = h.diff_matrix[j+1][i+1];
        }
}

void treeData::calc_dist_matrix(void)
{
    int i, j;

    for (i=0; i < ntaxa-1; i++)
    {
        for (j=i+1; j < ntaxa; j++)
        {
            // p-distance
            dist_matrix[i][j] = (double)diff_matrix[i][j]/num_sites;
            dist_matrix[j][i] = dist_matrix[i][j];
            if (dist == "jc")
            {
                // Jukes-Cantor distance
                dist_matrix[i][j] = -0.75*log(1.0-(4.0*dist_matrix[i][j]/3.0));
                dist_matrix[j][i] = dist_matrix[i][j];
            }
        }
    }
}

void treeData::tree_init(ptarray *treenode)
{
    int i, j;
    int nnodes = 2*ntaxa-1;
    node *p, *q;

    *treenode = (ptarray)malloc(nnodes*sizeof(node*));

    for (i=0; i < ntaxa; i++)
        (*treenode)[i] = (node*)malloc(sizeof(node));
    for (i=ntaxa; i < nnodes; i++)
    {
        q = 0;
        for (j=1; j <= 3; j++)
        {
            p = (node*)malloc(sizeof(node));
            p->next = q;
            q = p;
        }
        p->next->next->next = p;
        (*treenode)[i] = p;
    }
}

void treeData::setup_tree(tree *curtree)
{
    int i;
    int nnodes = 2*ntaxa-1;
    node *p;

    for (i=1; i <= nnodes; i++)
    {
        curtree->nodep[i-1]->back = 0;
        curtree->nodep[i-1]->tip = (i <= ntaxa);
        curtree->nodep[i-1]->index = i;
        curtree->nodep[i-1]->v = 0.0;
        if (i > ntaxa)
        {
            p = curtree->nodep[i-1]->next;
            while (p != curtree->nodep[i-1])
            {
                p->back = 0;
                p->tip = 0;
                p->index = i;
                p = p->next;
            }
        }
    }
    curtree->start = curtree->nodep[0];
}

void treeData::free_tree(ptarray *treenode)
{
    int i;
    node *p, *q;

    for (i=0; i < ntaxa; i++)
        free((*treenode)[i]);
    for (i=ntaxa; i < 2*ntaxa-1; i++)
    {
        p = (*treenode)[i];
        q = p->next;
        while (q != p)
        {
            node *r = q;
            q = q->next;
            free(r);
        }
        free(p);
    }
    free(*treenode);
}

std::string treeData::parseCommandLine(int argc, char *argv[])
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
    args >> GetOpt::Option('k', min_sites);
    args >> GetOpt::Option('w', win_size);
    args >> GetOpt::Option('d', dist);
    if (args >> GetOpt::OptionPresent('w'))
    {
        win_size *= 1000;
        flag |= BAM_WINDOW;
    }
    if (args >> GetOpt::OptionPresent('h'))
        flag |= BAM_HEADERIN;
    if (args >> GetOpt::OptionPresent('i'))
        flag |= BAM_ILLUMINA;
    if ((dist != "pdist") && (dist != "jc"))
    {
        msg = dist + " is not a valid distance option";
        fatal_error(msg, __FILE__, __LINE__, &treeUsage);
    }
    args >> GetOpt::GlobalOption(glob_opts);

    // run some checks on the command line

    // check to see if distance was specified, else use p-distance
    if (dist.empty())
        dist = "pdist";

    // if no input BAM file is specified -- print usage and exit
    if (glob_opts.size() < 2)
    {
        msg = "Need to specify BAM file name";
        fatal_error(msg, __FILE__, __LINE__, &treeUsage);
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
        fatal_error(msg, __FILE__, __LINE__, &treeUsage);
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

treeData::treeData(void)
{
    derived_type = TREE;
    refid = NULL;
    dist = "pdist";
    min_sites = 10;
    win_size = 0;
}

void treeData::init_tree(void)
{
    int i;
    int length = end-beg;
    int n = sm->n;
    int npops = sm->npops;

    ntaxa = n+1;
    segsites = 0;

    try
    {
        types = new unsigned long long [length]();
        pop_mask = new unsigned long long [npops]();
        pop_nsmpl = new unsigned char [npops]();
        pop_sample_mask = new unsigned long long [npops]();
        hap.pos = new unsigned int [length]();
        hap.idx = new unsigned int [length]();
        hap.ref = new unsigned char [length]();
        hap.seq = new unsigned long long* [n];
        hap.base = new unsigned char* [n];
        hap.rms = new unsigned short* [n];
        hap.snpq = new unsigned short* [n];
        hap.num_reads = new unsigned short* [n];
        diff_matrix = new unsigned short* [ntaxa];
        dist_matrix = new double* [ntaxa];
        for (i=0; i < n; i++)
        {
            hap.seq[i] = new unsigned long long [length]();
            hap.base[i] = new unsigned char [length]();
            hap.rms[i] = new unsigned short [length]();
            hap.snpq[i] = new unsigned short [length]();
            hap.num_reads[i] = new unsigned short [length]();
        }
        for (i=0; i < ntaxa; i++)
        {
            diff_matrix[i] = new unsigned short [ntaxa]();
            dist_matrix[i] = new double [ntaxa]();
        }
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
    }
}

void treeData::destroy_tree(void)
{
    int i;

    delete [] pop_mask;
    delete [] types;
    delete [] pop_nsmpl;
    delete [] pop_sample_mask;
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
    for (i=0; i < ntaxa; i++)
    {
        delete [] diff_matrix[i];
        delete [] dist_matrix[i];
    }
    delete [] hap.seq;
    delete [] hap.base;
    delete [] hap.snpq;
    delete [] hap.rms;
    delete [] hap.num_reads;
    delete [] diff_matrix;
    delete [] dist_matrix;
}

void treeData::treeUsage(void)
{
    std::cerr << std::endl;
    std::cerr << "Usage:   popbam tree [options] <in.bam> [region]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options: -i          base qualities are Illumina 1.3+     [ default: Sanger ]" << std::endl;
    std::cerr << "         -h  FILE    Input header file                    [ default: none ]" << std::endl;
    std::cerr << "         -d  STR     distance (pdist or jc)               [ default: pdist ]" << std::endl;
    std::cerr << "         -w  INT     use sliding window of size (kb)" << std::endl;
    std::cerr << "         -k  INT     minimum number of sites in window    [ default: 10 ]" << std::endl;
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
