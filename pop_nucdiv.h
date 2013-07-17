/** \file pop_nucdiv.h
 *  \brief Header for the pop_nucdiv.cpp file
 *  \author Daniel Garrigan
 *  \version 0.3
*/

#include "popbam.h"

///
/// Definitions
///

/*! \def EPSILON
 *  \brief Value of epsilon for testing whether double is zero
 */
#define EPSILON 1e-08

//
// Define data structures
//

/*!
 * \class nucdivData
 * \brief A derived class for passing parameters and data to the nucdiv function
 */
class nucdivData: public popbamData
{
	public:
		// constructor
		nucdivData();

		// destructor
		~nucdivData() {}

		friend void calc_diff_matrix(nucdivData&);

		// member public variables
		hData_t hap;                            //!< Structure to hold haplotype data
		unsigned int win_size;                  //!< Size of sliding window in kilobases
		unsigned long long **pop_sample_mask;   //!< Bit mask for samples covered from a specific population
		unsigned short *min_pop_n;              //!< Minimum sample size per population

		// member public functions
		void calc_nucdiv(void);
		std::string parseCommandLine(int, char**);
		void set_min_pop_n(void);
		void init_nucdiv(void);
		void destroy_nucdiv(void);
		void print_nucdiv(int);

	private:
		// member private variables
		int min_sites;                          //!< User-specified minimum number of aligned sites to perform analysis
		unsigned short **diff_matrix;           //!< Array of pairwise sequence differences
		double *piw;                            //!< Array of within-population nucleotide diversity
		double *pib;                            //!< Array of between-population Dxy values

		// member private functions
		static void nucdivUsage(void);
};

///
/// Function prototypes
///

/*!
 * \fn int make_nucdiv(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Runs the nucleotide diversity calculations
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int make_nucdiv(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);
