/** \file tables.h
 *  \brief Header within definitions of tables
 *  \author Daniel Garrigan
 *  \version 0.3
 */

/*! \def bam_nt16_table[256]
 *  \brief A lookup table
 */
 extern unsigned char bam_nt16_table[256];

/*! \def bam_nt16_nt4_table[16]
 *  \brief A reverse lookup table
 */
extern int bam_nt16_nt4_table[16];

/*! \def iupac[16]
 *  \brief A lookup table for IUPAC codes
 */
extern const char iupac[16];

/*! \def iupac[16]
 *  \brief A reverse lookup table for IUPAC codes
 */
unsigned char iupac_rev[256];
