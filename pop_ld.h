/** \file pop_ld.h
 *  \brief Header for the pop_ld.cpp file
 *  \author Daniel Garrigan
 *  \version 0.3
*/

#include "popbam.h"

///
/// Additional include headers
///
#include <vector>
#include <algorithm>

///
/// Definitions
///

//
// Define data structures
//

/*!
 * \class ldData
 * \brief A derived class for passing parameters and data to the ld function
 */
class ldData: public popbamData
{
	public:
		// constructor
		ldData();

		// destructor
		~ldData() {}

		// member variables
		int output;                             //!< Analysis output option
		unsigned int win_size;                  //!< Size of sliding window in kilobases
		unsigned long long *pop_sample_mask;    //!< Bit mask for samples covered from a specific population
		hData_t hap;                            //!< Structure to hold haplotype data
		int min_snps;                           //!< Minimum number of snps for a window to be considered
		int min_sites;                          //!< Minimum number of sites for a window to be considered
		unsigned short min_freq;                //!< Minimum allele count in LD calculation
		int *num_snps;                          //!< Number of SNPs in a given window
		double *omegamax;                       //!< Pointer to array of omega_max values
		double *wallb;                          //!< Pointer to array of Wall's B statistic
		double *wallq;                          //!< Pointer to array of Wall's Q statistic
		double *zns;                            //!< Pointer to array of ZnS values

		// member functions
		void calc_zns(void);
		void calc_omegamax(void);
		void calc_wall(void );
		std::string parseCommandLine(int, char**);
		void init_ld(void);
		void destroy_ld(void);
		void print_ld(int);
		static void ldUsage(void);
};

///
/// Function prototypes
///

/*!
 * \fn int make_ld(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Runs the linkage disequilibrium analysis
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int make_ld(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);

typedef void(ldData::*ld_func)(void);
