/** \file tables.h
 *  \brief Header within definitions of tables
 *  \author Daniel Garrigan
 *  \version 0.3
 */
#ifdef __cplusplus 
extern "C" {
#endif
/*! \def bam_nt16_table[256]
 *  \brief A lookup table
 */
 extern const unsigned char bam_nt16_table[256];

/*! \def bam_nt16_nt4_table[16]
 *  \brief A reverse lookup table
 */
extern const int bam_nt16_nt4_table[16];

/*! \def iupac[16]
 *  \brief A lookup table for IUPAC codes
 */
extern const char iupac[16];

/*! \def iupac[16]
 *  \brief A reverse lookup table for IUPAC codes
 */
extern const unsigned char iupac_rev[256];
#ifdef __cplusplus
}
#endif
