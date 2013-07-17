/** \file pop_sfs.h
 *  \brief Header for the pop_sfs.cpp file
 *  \author Daniel Garrigan
 *  \version 0.3
*/

#include "popbam.h"

///
/// Additional include headers
///
#include <limits>

///
/// Definitions
///

//
// Define data structures
//
/*!
 * \class sfsData
 * \brief A derived class for passing parameters and data to the sfs function
 */
class sfsData: public popbamData
{
	public:
		// constructor
		sfsData();

		// destructor
		~sfsData() {}

		// member variables
		unsigned int win_size;                  //!< Size of sliding window in kilobases
		unsigned long long *pop_sample_mask;    //!< Bit mask for samples covered from a specific population
		int min_sites;                          //!< User-specified minimum number of aligned sites to perform analysis
		int *num_snps;                          //!< Number of SNPs in a given window
		hData_t hap;                            //!< Structure to hold haplotype data
		std::string outgroup;                   //!< Sample name of outgroup to use
		int outidx;                             //!< Index of outgroup sequence
		double *a1;                             //!< Constants for Tajima's D calculation
		double *a2;                             //!< Constants for Tajima's D calculation
		double *e1;                             //!< Constants for Tajima's D calculation
		double *e2;                             //!< Constants for Tajima's D calculation
		double *td;                             //!< Pointer to the array of Tajima's D statistics
		double *fwh;                            //!< Pointer to the array of standardized Fay and Wu's H statistics

		// member functions
		std::string parseCommandLine(int, char**);
		void init_sfs(void);
		void destroy_sfs(void);
		void print_sfs(int);
		static void sfsUsage(void);
		void calc_sfs(void);
		void calc_a1(void);
		void calc_a2(void);
		void calc_e1(void);
		void calc_e2(void);
};

///
/// Function prototypes
///

/*!
 * \fn int make_sfs(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Runs the site frequency spectrum analysis
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int make_sfs(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);
