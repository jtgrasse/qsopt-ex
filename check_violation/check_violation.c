/* check_violation: given an LP/MPS file and a basis file, use qsopt-ex's
 * exact-arithmetic routines to test whether the *given* basis is an
 * optimality certificate. On failure, report which primal/dual constraint
 * is violated and by how much.
 *
 * usage: ./check_violation <mps_file> <bas_file>
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "QSopt_ex.h"
#include "exact.h"
#include "except.h"
#include "lpdata_mpq.h"
#include "lpdefs_mpq.h"
#include "fct_mpq.h"
#include "basis_mpq.h"

static void usage(const char *prog)
{
	fprintf(stderr,
		"Usage: %s [-L] <lp_or_mps_file> <basis_file>\n"
		"  Reads the LP (MPS by default; -L for LP format) plus a qsopt-ex\n"
		"  basis file.\n",
		prog);
}

int main(int argc, char **argv)
{
	int rval = 0;
	int lpfmt = 0;
	const char *probfile = NULL;
	const char *basisfile = NULL;
	mpq_QSdata *p = NULL;
	QSbasis *basis = NULL;
	mpq_t *x_struct = NULL;
	mpq_t *y_rows = NULL;
	mpq_t *xz_full = NULL;
	int singular = 0;
	int is_optimal = 0;

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-L")) lpfmt = 1;
		else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(argv[0]); return 0; }
		else if (!probfile) probfile = argv[i];
		else if (!basisfile) basisfile = argv[i];
		else { usage(argv[0]); return 1; }
	}
	if (!probfile || !basisfile) { usage(argv[0]); return 1; }

	QSexactStart();
	QSexact_set_precision(128);

	p = mpq_QSread_prob(probfile, lpfmt ? "LP" : "MPS");
	if (p) p->simplex_display = 1;
	if (!p) { fprintf(stderr, "failed to read problem: %s\n", probfile); return 2; }

	rval = mpq_QSread_and_load_basis(p, basisfile);
	if (rval) { fprintf(stderr, "failed to read/load basis: %s\n", basisfile); return 3; }

	basis = mpq_QSget_basis(p);
	if (!basis) { fprintf(stderr, "failed to retrieve basis\n"); return 4; }

	printf("== Checking basis %s against problem %s ==\n", basisfile, probfile);
	fflush(stdout);

	/* prep the exact lpinfo, load and factor the given basis in exact arithmetic.
	 * we test the input basis */
	if (mpq_QSload_basis(p, basis)) {
		fprintf(stderr, "mpq_QSload_basis failed\n"); goto done;
	}
	if (p->cache) {
		mpq_ILLlp_cache_free(p->cache);
		mpq_clear(p->cache->val);
		free(p->cache);
		p->cache = NULL;
	}
	p->qstatus = QS_LP_MODIFIED;
	if (p->qslp->sinfo) { mpq_ILLlp_sinfo_free(p->qslp->sinfo); free(p->qslp->sinfo); p->qslp->sinfo = NULL; }
	if (p->qslp->rA)    { mpq_ILLlp_rows_clear(p->qslp->rA);    free(p->qslp->rA);    p->qslp->rA    = NULL; }

	mpq_free_internal_lpinfo(p->lp);
	mpq_init_internal_lpinfo(p->lp);
	if (mpq_build_internal_lpinfo(p->lp)) {
		fprintf(stderr, "mpq_build_internal_lpinfo failed\n"); goto done;
	}
	mpq_ILLfct_set_variable_type(p->lp);
	if (mpq_ILLbasis_load(p->lp, p->basis)) {
		fprintf(stderr, "mpq_ILLbasis_load failed\n"); goto done;
	}
	if (mpq_ILLbasis_factor(p->lp, &singular) || singular) {
		fprintf(stderr, "basis is %s\n", singular ? "SINGULAR" : "unfactorable");
		is_optimal = 0;
		goto verdict;
	}

	memset(&(p->lp->basisstat), 0, sizeof(mpq_lp_status_info));
	mpq_ILLfct_compute_xbz(p->lp);
	mpq_ILLfct_compute_piz(p->lp);
	mpq_ILLfct_compute_dz(p->lp);

	{
		mpq_ILLlpdata *qslp = p->lp->O;
		const int ncols = qslp->ncols;
		const int nrows = qslp->nrows;
		const int nstruct = qslp->nstruct;

		xz_full  = mpq_EGlpNumAllocArray(ncols);
		/* QSexact_optimal_test writes into p_sol[nstruct + i]
		 * so allocate nstruct + nrows, not just nstruct. */
		x_struct = mpq_EGlpNumAllocArray(nstruct + nrows);
		y_rows   = mpq_EGlpNumAllocArray(nrows);

		for (int i = 0; i < nrows; i++)
			mpq_set(xz_full[p->lp->baz[i]], p->lp->xbz[i]);
		for (int j = 0; j < p->lp->nnbasic; j++) {
			int col = p->lp->nbaz[j];
			if (p->lp->vstat[col] == STAT_UPPER)
				mpq_set(xz_full[col], p->lp->uz[col]);
			else if (p->lp->vstat[col] == STAT_LOWER)
				mpq_set(xz_full[col], p->lp->lz[col]);
			else
				mpq_set_ui(xz_full[col], 0UL, 1UL);
		}
		for (int i = 0; i < nstruct; i++)
			mpq_set(x_struct[i], xz_full[qslp->structmap[i]]);
		// flip sign for maximization
		int flip = (qslp->objsense != QS_MIN);
		for (int i = 0; i < nrows; i++) {
			mpq_set(y_rows[i], p->lp->piz[i]);
			if (flip) mpq_neg(y_rows[i], y_rows[i]);
		}
	}

	/* QSexact_optimal_test returns 1 if (x, y, basis) is a certified optimum,
	 * 0 otherwise. It always prints exact-rational max primal/dual violations. */
	is_optimal = QSexact_optimal_test2(p, x_struct, y_rows, basis);

verdict:
	printf("\n== Verifier verdict ==\n");
	if (is_optimal) {
		mpq_t dobj;
		mpq_init(dobj);
		if (!mpq_QSget_objval(p, &dobj)) {
			char *s = mpq_get_str(NULL, 10, dobj);
			printf("  RESULT: basis is a certified optimum.\n");
			printf("  objective value = %s  (~= %.17g)\n", s, mpq_get_d(dobj));
			void (*freefn)(void *, size_t);
			mp_get_memory_functions(NULL, NULL, &freefn);
			freefn(s, strlen(s) + 1);
		} else {
			printf("  RESULT: basis is a certified optimum.\n");
		}
		mpq_clear(dobj);
	} else {
		printf("  RESULT: basis FAILED verification.\n"
			   "          See 'VERIFICATION SUMMARY' above for exact-rational\n"
			   "          max primal violation and max dual violation.\n");
	}

done:
	if (xz_full)  mpq_EGlpNumFreeArray(xz_full);
	if (x_struct) mpq_EGlpNumFreeArray(x_struct);
	if (y_rows)   mpq_EGlpNumFreeArray(y_rows);
	if (basis) mpq_QSfree_basis(basis);
	if (p) mpq_QSfree_prob(p);
	QSexactClear();
	return is_optimal ? 0 : 1;
}
