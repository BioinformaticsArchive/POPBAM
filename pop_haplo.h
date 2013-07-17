/** \file pop_haplo.h
 *  \brief Header for the pop_haplo.cpp file
 *  \author Daniel Garrigan
 *  \version 0.3
*/

#include "popbam.h"

///
/// Additional include headers
///
#include <limits>
#include <algorithm>
#include <vector>
#include <list>

///
/// Definitions
///

//
// Define data structures
//

/*!
 * \class haploData
 * \brief A derived class for passing parameters and data to the haplo function
 */
class haploData: public popbamData
{
	public:
		// constructor
		haploData();

		// destructor
		~haploData() {}

		friend void calc_diff_matrix(haploData&);

		// member public variables
		hData_t hap;                            //!< Structure to hold haplotype data
		unsigned int win_size;                  //!< User-specified window size in kilobases
		unsigned long long *pop_sample_mask;    //!< Bit mask for samples covered from a specific population

		// member public functions
		std::string parseCommandLine(int, char**);
		void init_haplo(void);
		void calc_haplo(void);
		void print_haplo(int);
		void destroy_haplo(void);

	private:
		// member private variables
		int output;                             //!< Analysis output option
		int min_sites;                          //!< User-specified minimum number of aligned sites to perform analysis
		int *nhaps;                             //!< Array of number of haplotypes within populations
		double *hdiv;                           //!< Array of haplotype diversities within populations
		unsigned short **diff_matrix;           //!< Array of pairwise sequence differences
		unsigned short *minDxy;                 //!< Array of minimum between-population Dxy
		double *piw;                            //!< Array of within-population heterozygosities
		double *pib;                            //!< Array of between-population heterozygosities
		double *ehhs;                           //!< Array of site-specfic haplotype homozygosity statistics within populations

		// member private function
		void calc_nhaps(void);
		void calc_ehhs(void);
		void calc_minDxy(void);
		static void haploUsage(void);
};

///
/// Function prototypes
///

/*!
 * \fn int make_haplo(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Calculate haplotype-based statistics
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int make_haplo(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);

typedef void(haploData::*haplo_func)(void);
