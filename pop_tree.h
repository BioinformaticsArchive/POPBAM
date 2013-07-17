/** \file pop_tree.h
 *  \brief Header for the pop_tree.cpp file
 *  \author Daniel Garrigan
 *  \version 0.3
 * Much of the code for the NJ algorithm heavily borrows from the PHYLIP
 * package (v 3.6) written by Mary Kuhner, Jon Yamato, Joseph Felsenstein,
 * Akiko Fuseki, Sean Lamont, and Andrew Keefe at the University of Washington
*/

#include "popbam.h"

//
// Define data structures
//

/*!
 * \struct node
 * \brief A structure to represent an individual node in a bifurcating tree
 */
typedef struct _node {
	struct _node *next;               //!< Pointer to the next node in the list
	struct _node *back;               //!< Pointer to the previous node in the list
	bool tip;                         //!< Is the node a tip?
	int index;                        //!< Index number of the node
	double v;                         //!< Branch length above the node
} node;

/*!
 * \var typedef node **ptarray
 * \brief A pointer to a node pointer
 */
typedef node **ptarray;

/*!
 * \struct tree
 * \brief A structure to represent a bifurcating tree
 */
typedef struct _tree {
	node *start;                      //!< Pointer to the start node of the tree
	ptarray nodep;                    //!< Pointer to the linked list of nodes in the tree
} tree;

/*!
 * \class treeData
 * \brief A derived class for passing parameters and data to the tree function
 */
class treeData: public popbamData
{
	public:
		// constructor
		treeData();

		// destructor
		~treeData() {}

		friend void calc_diff_matrix(treeData&);

		// member public variables
		hData_t hap;                            //!< Structure to hold haplotype data (public)
		unsigned int win_size;                  //!< Size of sliding window in kilobases (public)
		unsigned long long *pop_sample_mask;    //!< Bit mask for samples covered from a specific population
		char *refid;                            //!< Pointer to string that holds the reference sequence name
		std::string dist;                       //!< Pointer to the name of the desired distance metric	(-d switch)

		// member public functions
		void make_nj(int);
		void calc_dist_matrix(void);
		std::string parseCommandLine(int, char**);
		void init_tree(void);
		void destroy_tree(void);

	private:
		// member private variables
		unsigned short **diff_matrix;           //!< Array of pairwise sequence differences
		double **dist_matrix;                   //!< Array of divergence calculations
		int min_sites;                          //!< User-specified minimum number of aligned sites to perform analysis
		int ntaxa;                              //!< Total number of tips in the tree
		int *enterorder;                        //!< Array containing the input order of OTUs for the NJ algorithm

		// member private functions
		void join_tree(tree, node**);
		void print_tree(node*, node*);
		void hookup(node*, node*);
		void setup_tree(tree*);
		void tree_init(ptarray *);
		void free_tree(ptarray *);
		static void treeUsage(void);
};

///
/// Function prototypes
///

/*!
 * \fn int make_tree(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Runs the neighbor-joining tree construction procedure
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int make_tree(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);
