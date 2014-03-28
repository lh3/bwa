#ifndef __AC_KSW_H
#define __AC_KSW_H

#include <stdint.h>

#define KSW_XBYTE  0x10000
#define KSW_XSTOP  0x20000
#define KSW_XSUBO  0x40000
#define KSW_XSTART 0x80000

struct _kswq_t;
typedef struct _kswq_t kswq_t;

typedef struct {
	int score; // best score
	int te, qe; // target end and query end
	int score2, te2; // second best score and ending position on the target
	int tb, qb; // target start and query start
} kswr_t;

#ifdef __cplusplus
extern "C" {
#endif

	/**
	 * Aligning two sequences
	 *
	 * @param qlen    length of the query sequence (typically <tlen)
	 * @param query   query sequence with 0 <= query[i] < m
	 * @param tlen    length of the target sequence
	 * @param target  target sequence
	 * @param m       number of residue types
	 * @param mat     m*m scoring matrix in one-dimension array
	 * @param gapo    gap open penalty; a gap of length l cost "-(gapo+l*gape)"
	 * @param gape    gap extension penalty
	 * @param xtra    extra information (see below)
	 * @param qry     query profile (see below)
	 *
	 * @return        alignment information in a struct; unset values to -1
	 *
	 * When xtra==0, ksw_align() uses a signed two-byte integer to store a
	 * score and only finds the best score and the end positions. The 2nd best
	 * score or the start positions are not attempted. The default behavior can
	 * be tuned by setting KSW_X* flags:
	 *
	 *   KSW_XBYTE:  use an unsigned byte to store a score. If overflow occurs,
	 *               kswr_t::score will be set to 255
	 *
	 *   KSW_XSUBO:  track the 2nd best score and the ending position on the
	 *               target if the 2nd best is higher than (xtra&0xffff)
	 *
	 *   KSW_XSTOP:  stop if the maximum score is above (xtra&0xffff)
	 *
	 *   KSW_XSTART: find the start positions
	 *
	 * When *qry==NULL, ksw_align() will compute and allocate the query profile
	 * and when the function returns, *qry will point to the profile, which can
	 * be deallocated simply by free(). If one query is aligned against multiple
	 * target sequences, *qry should be set to NULL during the first call and
	 * freed after the last call. Note that qry can equal 0. In this case, the
	 * query profile will be deallocated in ksw_align().
	 */
	kswr_t ksw_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra, kswq_t **qry);
	kswr_t ksw_align2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int xtra, kswq_t **qry);

	/**
	 * Banded global alignment
	 *
	 * @param qlen    query length
	 * @param query   query sequence with 0 <= query[i] < m
	 * @param tlen    target length
	 * @param target  target sequence with 0 <= target[i] < m
	 * @param m       number of residue types
	 * @param mat     m*m scoring mattrix in one-dimension array
	 * @param gapo    gap open penalty; a gap of length l cost "-(gapo+l*gape)"
	 * @param gape    gap extension penalty
	 * @param w       band width
	 * @param n_cigar (out) number of CIGAR elements
	 * @param cigar   (out) BAM-encoded CIGAR; caller need to deallocate with free()
	 *
	 * @return        score of the alignment
	 */
	int ksw_global(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *n_cigar, uint32_t **cigar);
	int ksw_global2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int *n_cigar, uint32_t **cigar);

	/**
	 * Extend alignment
	 *
	 * The routine aligns $query and $target, assuming their upstream sequences,
	 * which are not provided, have been aligned with score $h0. In return,
	 * region [0,*qle) on the query and [0,*tle) on the target sequences are
	 * aligned together. If *gscore>=0, *gscore keeps the best score such that
	 * the entire query sequence is aligned; *gtle keeps the position on the
	 * target where *gscore is achieved. Returning *gscore and *gtle helps the
	 * caller to decide whether an end-to-end hit or a partial hit is preferred.
	 *
	 * The first 9 parameters are identical to those in ksw_global()
	 *
	 * @param h0      alignment score of upstream sequences
	 * @param _qle    (out) length of the query in the alignment
	 * @param _tle    (out) length of the target in the alignment
	 * @param _gtle   (out) length of the target if query is fully aligned
	 * @param _gscore (out) score of the best end-to-end alignment; negative if not found
	 *
	 * @return        best semi-local alignment score
	 */
	int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);
	int ksw_extend2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);

#ifdef __cplusplus
}
#endif

#endif
