/** \file pop_snp.h
 *  \brief Header for the pop_snp.cpp file
 *  \author Daniel Garrigan
 *  \version 0.3
*/

#include "popbam.h"

char bam_nt16_rev_table[16] = {'=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};

/*!
 * \class snpData
 * \brief A derived class for passing parameters and data to the snp function
 */
class snpData: public popbamData
{
	public:
		// constructor
		snpData();

		// destructor
		~snpData() {}

		// member public variables
		hData_t hap;                            //!< Structure to hold haplotype data
		unsigned int win_size;                  //!< User-specified window size in kilobases
		unsigned long long **pop_sample_mask;   //!< Bit mask for samples covered from a specific population
		int output;                             //!< User-specified output mode
		std::string outgroup;                   //!< Sample name of outgroup to use
		int outidx;                             //!< Index of outgroup sequence

		// member public functions
		std::string parseCommandLine(int, char**);
		void init_snp(void);
		void print_ms_header(long);
		void print_snp(int);
		void destroy_snp(void);

	private:
		// member private functions
		void print_popbam_snp(int);
		void print_sweep(int);
		void print_ms(int);
		static void snpUsage(void);
};

///
/// Function prototypes
///

/*!
 * \fn int make_snp(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Runs the SNP analysis
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int make_snp(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);

typedef void(snpData::*snp_func)(int);
