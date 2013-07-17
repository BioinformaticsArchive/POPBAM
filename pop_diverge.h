/** \file pop_diverge.h
 *  \brief Header for the pop_diverge.cpp file
 *  \author Daniel Garrigan
 *  \version 0.3
*/

#include "popbam.h"

//
// Define data structures
//

/*!
 * \class divergeData
 * \brief A derived class for passing parameters and data to the diverge function
 */
class divergeData: public popbamData
{
	public:
		// constructor
		divergeData();

		// destructor
		~divergeData() {}

		// member public variables
		int output;                             //!< Analysis output option
		std::string outgroup;                   //!< Sample name of outgroup to use
		int outidx;                             //!< Index of outgroup sequence
		hData_t hap;                            //!< Structure to hold haplotype data
		unsigned int win_size;                  //!< Size of sliding window in kilobases
		unsigned long long *pop_sample_mask;    //!< Bit mask for samples covered from a specific population
		std::string dist;                       //!< Pointer to the name of the desired distance metric	(-d switch)
		unsigned short *min_pop_n;              //!< Minimum sample size per population
		int *num_snps;                          //!< Number of SNPs in a given window

		// member public functions
		std::string parseCommandLine(int, char**);
		void calc_diverge(void);
		void init_diverge(void);
		void set_min_pop_n(void);
		void destroy_diverge(void);
		void print_diverge(int);

	private:
		// member private variables
		int min_sites;                          //!< User-specified minimum number of aligned sites to perform analysis
		unsigned short *pop_div;                //!< Array of mean population divergence calculations
		unsigned short *ind_div;                //!< Array of individual divergence calculations

		// member private functions
		static void divergeUsage(void);
};

///
/// Function prototypes
///

/*!
 * \fn int make_diverge(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Calculate divergence with reference genome sequence
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int make_diverge(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);
