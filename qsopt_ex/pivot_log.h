/* Pivot log: write the exact (rational) LP that QSopt_ex checks optimality
 * against -- constraint matrix A including the logical (slack, surplus and
 * artificial) columns, right-hand side b, objective c and bounds -- followed
 * by the starting basis and every basis change made by the simplex, in every
 * precision (dbl, mpf, mpq).
 *
 * Column indices in the log are canonical and 0-based: structural variable j
 * is column j and the logical of row i is column nstruct + i.  This is the
 * order of QSbasis (cstat then rstat).  QSopt_ex's internal column order can
 * differ between copies of the LP (the scaled copy puts the logicals first),
 * so the template code fills a per-LP map from internal to canonical index
 * before handing basis arrays to this module. */
#ifndef QS_PIVOT_LOG_H
#define QS_PIVOT_LOG_H

#include <gmp.h>

/* nonzero while a log file is open; checked by every hook so that logging
 * costs nothing when disabled */
extern int pivot_log_active;

// open the log and write the header; returns nonzero on failure
int  pivot_log_open (const char *path, const char *probname);
void pivot_log_close (void);

// pause and resume all logging (nests), e.g. around trial factorizations
void pivot_log_suspend (void);
void pivot_log_resume (void);

// write a "% ..." comment line
void pivot_log_comment (const char *fmt, ...)
	__attribute__ ((format (printf, 1, 2)));
// note why the basis is about to change; printed with the next full basis
// written by pivot_log_sync, and discarded if the basis turns out unchanged
void pivot_log_set_reason (const char *fmt, ...)
	__attribute__ ((format (printf, 1, 2)));

// write A, b, c and the bounds of the mpq LP, preceded by comment "what"
struct mpq_ILLlpdata;
void pivot_log_problem (const struct mpq_ILLlpdata *qslp, const char *what);

// write the exact solution in canonical column order: primal x and reduced
// costs d over all ncols columns, duals y over the nrows rows, objective z
void pivot_log_solution (int ncols, int nrows, const mpq_t *x, const mpq_t *y,
		const mpq_t *d, const mpq_t z);

/* map from internal to canonical column index of the LP identified by
 * owner.  On return *fill is 1 if the caller must fill map[0..ncols-1];
 * pass *fill = 1 to force that. */
int *pivot_log_colmap (const void *owner, int ncols, int *fill);

// basis set (internal indices) may have changed without a pivot: write it
// as a full basis if it differs from the last basis written
void pivot_log_sync (int nrows, const int *baz);
// internal column lvar left the basis and evar entered; baz is the basis
// after the change
void pivot_log_pivot (int nrows, const int *baz, int lvar, int evar);

#endif
