/* ========================================================================= */
/* ESolver "Exact Mixed Integer Linear Solver" provides some basic structures
 * and algorithms commons in solving MIP's
 *
 * Copyright (C) 2008 David Applegate, Bill Cook, Sanjeeb Dash, Daniel Espinoza.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 * */
/* ========================================================================= */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <signal.h>
#include <limits.h>
#include <time.h>
#include <sys/resource.h>

#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#include "QSopt_ex.h"

#include "except.h"
#include "logging-private.h"
#include "qs_config.h"
#include "timing_log.h"
#include "basis_snapshot.h"
#include "pivot_log.h"
#include "factor_mpq.h"

/* ========================================================================= */
/** @name static parameters for the main program */
/*@{*/
static char *fname = 0;
static int lpfile = 0;
static int usescaling = 1;
static int showversion = 0;
static int simplexalgo = PRIMAL_SIMPLEX;
static int pstrategy = QS_PRICE_PSTEEP;
static int dstrategy = QS_PRICE_DSTEEP;
static unsigned precision = 128;
static int printsol = 0;
static char *solname = 0;
static char *readbasis = 0;
static char *writebasis = 0;
static int dump_bases = 0;
static char *pivotlog = 0;
/** @brief maximum running time */
static double max_rtime = INT_MAX;
/** @brief maximum memory usage */
static unsigned long memlimit = UINT_MAX;
/*@}*/
/* ========================================================================= */
/** @brief Display options to the screen */
static void usage (char *s)
{
	fprintf (stderr, "Usage: %s [- below -] prob_file\n", s);
	fprintf (stderr, "   -b f  write basis to file f\n");
	fprintf (stderr, "   -B f  read initial basis from file f\n");
	fprintf (stderr, "   -D    dump basis matrix snapshots to <problem>_Bases/\n");
#if 0
	fprintf (stderr, "   -I    solve the MIP using BestBound\n");
	fprintf (stderr, "   -E    edit problem after solving initial version\n");
#endif
	fprintf (stderr, "   -L    input file is in lp format (default: mps)\n");
	fprintf (stderr, "   -N n  use 'n.log' as the log file instead of default\n");
	fprintf (stderr, "   -O    write the final solution to the given file\n");
	fprintf (stderr, "         append .gz/.bz2 to the .sol extension to compress the file\n");
	fprintf (stderr, "   -p #  run primal simplex with pricing rule #\n");
	fprintf (stderr,
					 "         (%d-Dantzig, %d-Devex, %d-Steep (default), %d-Partial\n",
					 QS_PRICE_PDANTZIG, QS_PRICE_PDEVEX, QS_PRICE_PSTEEP,
					 QS_PRICE_PMULTPARTIAL);
	fprintf (stderr,
					 "   -P #  number of bits to use for the float representation (default: 128)\n");
	fprintf (stderr, "   -d #  run dual simplex with pricing rule #\n");
	fprintf (stderr, "         (%d-Dantzig, %d-Steep, %d-Partial, %d-Devex)\n",
					 QS_PRICE_DDANTZIG, QS_PRICE_DSTEEP, QS_PRICE_DMULTPARTIAL,
					 QS_PRICE_DDEVEX);
	fprintf (stderr, "   -S    do NOT scale the initial LP\n");
	fprintf (stderr, "   -T    enable timing output (qsopt_timing.log, time_precision_data, basis_scaling)\n");
	fprintf (stderr, "   -v    print QSopt version number\n");
	fprintf (stderr, "   -W f  write the exact LP (A, b, c, bounds with logicals), the starting\n");
	fprintf (stderr, "         basis and every basis change of the solve to file f\n");
	fprintf (stderr, "   -R n  maximum running time allowed, default %lf\n",
						max_rtime);
	fprintf (stderr, "   -m n  maximum memory usage allowed, default %lu\n",
						memlimit);
}

/* ========================================================================= */
/** @brief decide if a given file is mps or lp (only by extension) */
// AP: ftype = 0 if MPS and 1 if LP
static void get_ftype (char const *const name, int *ftype)
{
	char buff[4096],*argv[128];
	int argc;
	*ftype = 0; /* by default, file is MPS */
	snprintf(buff,4096,"%s",name);
	EGioNParse(buff,128,"."," ",&argc,argv);
	argc-=1;

	if(argc)
	{
		if(strncmp(argv[argc],"gz",3)==0) argc-=1;
		else if(strncmp(argv[argc],"GZ",3)==0) argc-=1;
		else if(strncmp(argv[argc],"bz2",4)==0) argc-=1;
		else if(strncmp(argv[argc],"BZ2",4)==0) argc-=1;
	}
	if(argc)
	{
		if(strncmp(argv[argc],"lp",3)==0) *ftype=1;
		else if(strncmp(argv[argc],"LP",3)==0) *ftype=1;
	}
}
/* ========================================================================= */
/** @brief signal handler for time-limit reached */
static void sighandler(int s)
{
	switch(s)
	{
		case SIGXCPU:
			fprintf(stderr,"TIME_LIMIT_REACHED (ending now)\n");
			exit(EXIT_FAILURE);
		default:
			fprintf(stderr,"Unknown signal %d (ending now)\n",s);
			exit(EXIT_FAILURE);
	}
}
/* ========================================================================= */
/** @brief function to handle resource usage limits */
static int mem_limits(void)
{
	int rval = 0;
	struct rlimit mlim;
	rval = getrlimit(RLIMIT_CPU,&mlim);
	CHECKRVAL(rval);

	// AP: added line below, with correct printing
	fprintf(stderr, "Cur rtime limit %ju, trying to set to %lg\n", (intmax_t) mlim.rlim_cur, max_rtime);
		// fprintf(stderr, "Cur rtime limit %ld, trying to set to %lg\n", mlim.rlim_cur, max_rtime);
	if(max_rtime > mlim.rlim_max) max_rtime = (double)mlim.rlim_max;
	mlim.rlim_cur = (rlim_t)max_rtime;

	// const struct rlimit lim_test = {
    // 	rlim_t rlim_cur = mlim.rlim_cur; /* Soft limit */
    // 	rlim_t rlim_max = mlim.rlim_max; /* Hard limit (ceiling for rlim_cur) */
	// };

	rval = setrlimit(RLIMIT_CPU,&mlim);
	TESTERRNOIF(rval);
	fprintf(stderr, "New rtime limit %ju (%.3lg)\n", (intmax_t) mlim.rlim_cur, max_rtime);
		// fprintf(stderr, "New rtime limit %ld (%.3lg)\n", mlim.rlim_cur, max_rtime);

	rval = getrlimit(RLIMIT_DATA,&mlim);
	TESTERRNOIF(rval);
	fprintf(stderr, "Cur data limit %ld,%ld (soft,hard)\n", mlim.rlim_cur, mlim.rlim_max);

	mlim.rlim_cur = memlimit;
	rval = setrlimit(RLIMIT_DATA,&mlim);
	TESTERRNOIF(rval);
	rval = getrlimit(RLIMIT_DATA,&mlim);
	TESTERRNOIF(rval);
	fprintf(stderr, "New data limit %ld,%ld (soft,hard)\n", mlim.rlim_cur, mlim.rlim_max);

	rval = getrlimit(RLIMIT_AS,&mlim);
	TESTERRNOIF(rval);
	fprintf(stderr, "Cur address space limit %ld,%ld (soft,hard)\n",
					mlim.rlim_cur, mlim.rlim_max);
	mlim.rlim_cur = memlimit;
	rval = setrlimit(RLIMIT_AS,&mlim);
	TESTERRNOIF(rval);
	rval = getrlimit(RLIMIT_AS,&mlim);
	TESTERRNOIF(rval);
	fprintf(stderr, "New address space limit %ld,%ld (soft,hard)\n",
					mlim.rlim_cur, mlim.rlim_max);
	mlim.rlim_cur = 0;
	rval = setrlimit(RLIMIT_CORE,&mlim);
	TESTERRNOIF(rval);
	rval = getrlimit(RLIMIT_CORE,&mlim);
	TESTERRNOIF(rval);
	fprintf(stderr, "New core dump space limit %ld,%ld (soft,hard)\n",
					mlim.rlim_cur, mlim.rlim_max);
	/* set signal handler for SIGXCPU */
	signal(SIGXCPU,sighandler);
	return rval;
}
/* ========================================================================= */
/** @brief parssing options for the program */
static int parseargs (int ac, char **av)
{
	int c;
	int boptind = 1;
	char *boptarg = 0;

	while ((c = ILLutil_bix_getopt (ac, av, "b:B:d:DEILm:O:p:P:R:STvN:W:", &boptind, &boptarg)) != EOF)
		switch (c)
		{
		case 'm':
			memlimit = strtoul(boptarg,0,10);
			break;
		case 'R':
			max_rtime = strtod(boptarg,0);
			break;
		case 'b':
			writebasis = boptarg;
			break;
		case 'B':
			readbasis = boptarg;
			break;
		case 'P':
			precision = atoi (boptarg);
			break;
		case 'd':
			simplexalgo = DUAL_SIMPLEX;
			dstrategy = atoi (boptarg);
			break;
		case 'D':
			dump_bases = 1;
			break;
		case 'L':
			lpfile = 1;
			break;
		case 'O':
			printsol = 1;
			solname = strdup(boptarg);
			break;
		case 'p':
			simplexalgo = PRIMAL_SIMPLEX;
			pstrategy = atoi (boptarg);
			break;
		case 'S':
			usescaling = 0;
			break;
		case 'T':
			set_timing_enabled(1);
			break;
		case 'v':
			showversion = 1;
			break;
		case 'N':
			set_log_file(boptarg);
			break;
		case 'W':
			pivotlog = boptarg;
			break;
		case '?':
		default:
			usage (av[0]);
			return 1;
		}

	if ((boptind == ac) && (showversion))
	{
		char *buf = 0;
		buf = mpq_QSversion ();
		printf ("%s\n", buf);
		mpq_QSfree ((void *) buf);
		exit(0);
	}

	if (boptind != (ac - 1))
	{
		usage (av[0]);
		return 1;
	}

	fname = av[boptind++];
	fprintf (stderr, "Reading problem from %s\n", fname);

	mem_limits();
	return 0;
}

/* get problem name */
static void derive_problem_name (const char *path, char *out, size_t outsz)
{
	const char *base = strrchr (path, '/');
	base = base ? base + 1 : path;
	snprintf (out, outsz, "%s", base);
	size_t n = strlen (out);
	if (n > 3 && strcmp (out + n - 3, ".gz") == 0)        { out[n-3] = 0; n -= 3; }
	else if (n > 4 && strcmp (out + n - 4, ".bz2") == 0)  { out[n-4] = 0; n -= 4; }
	char *dot = strrchr (out, '.');
	if (dot) *dot = 0;
}

/* Write a sparse rational matrix in COO triplet form.
 *
 * File format:
 *   nrows ncols nnz
 *   row col val            (repeated nnz times; val is a GMP rational num/den)
 */
static int write_sparse_mpq_matrix (const char *path, int nrows, int ncols,
		const mpq_t *matval, const int *matbeg, const int *matcnt, const int *matind)
{
	EGioFile_t *f = EGioOpen (path, "w");
	if (!f) {
		fprintf (stderr, "could not open %s for writing\n", path);
		return 1;
	}
	int nnz = 0;
	for (int j = 0; j < ncols; ++j) nnz += matcnt[j];

	EGioPrintf (f, "%d %d %d\n", nrows, ncols, nnz);
	char numbuf[8192];
	for (int j = 0; j < ncols; ++j) {
		int beg = matbeg[j], cnt = matcnt[j];
		for (int k = 0; k < cnt; ++k) {
			int row = matind[beg + k];
			const mpq_t *v = &matval[beg + k];
			gmp_snprintf (numbuf, sizeof (numbuf), "%Qd", *v);
			EGioPrintf (f, "%d %d %s\n", row, j, numbuf);
		}
	}

	EGioClose (f);
	return 0;
}

/* Write a dense rational vector as integer-scaled values, one per entry.
 *
 * File format:
 *   n
 *   L_0 int_val_0
 *   L_1 int_val_1
 *   ...
 * where the original entry is int_val_i / L_i.  L_i is just the denominator
 * of the canonical mpq_t (since each entry is scaled independently).
 */
static int write_dense_mpq_vector (const char *path, int n, const mpq_t *vec)
{
	EGioFile_t *f = EGioOpen (path, "w");
	if (!f) {
		fprintf (stderr, "could not open %s for writing\n", path);
		return 1;
	}
	EGioPrintf (f, "%d\n", n);
	char dbuf[8192], nbuf[8192];
	for (int i = 0; i < n; ++i) {
		gmp_snprintf (dbuf, sizeof (dbuf), "%Zd", mpq_denref (vec[i]));
		gmp_snprintf (nbuf, sizeof (nbuf), "%Zd", mpq_numref (vec[i]));
		EGioPrintf (f, "%s %s\n", dbuf, nbuf);
	}
	EGioClose (f);
	return 0;
}

/* Exact-arithmetic validation of one snapshot: factor B = A(:, baz) with
 * the library's own rational LU (mpq_ILLfactor).  A float simplex run can
 * pivot through EXACTLY singular bases without noticing -- its pivot
 * tolerances hide the dependency (e.g. two basic variables whose A columns
 * are identical), and the run is later abandoned for a higher precision
 * anyway -- so the pivot-trail tail is not automatically usable as test
 * data.  Returns 1 iff the basis is exactly nonsingular. */
static int snapshot_exactly_nonsingular (mpq_ILLlpdata *qslp, const int *baz)
{
	int rval = 0, ok = 0;
	int nsing = 0, *singr = 0, *singc = 0;
	int nrows = qslp->nrows;
	int *bcopy = 0;
	mpq_factor_work f;

	memset (&f, 0, sizeof (f));
	mpq_EGlpNumInitVar (f.fzero_tol);
	mpq_EGlpNumInitVar (f.szero_tol);
	mpq_EGlpNumInitVar (f.partial_tol);
	mpq_EGlpNumInitVar (f.maxelem_orig);
	mpq_EGlpNumInitVar (f.maxelem_factor);
	mpq_EGlpNumInitVar (f.maxelem_cur);
	mpq_EGlpNumInitVar (f.partial_cur);
	mpq_ILLfactor_init_factor_work (&f);
	rval = mpq_ILLfactor_create_factor_work (&f, nrows);
	if (rval) goto CLEANUP;

	// ILLfactor takes a non-const basis; keep the ring slot pristine
	ILL_SAFE_MALLOC (bcopy, nrows, int);
	memcpy (bcopy, baz, sizeof (int) * nrows);

	rval = mpq_ILLfactor (&f, bcopy, qslp->A.matbeg, qslp->A.matcnt,
			qslp->A.matind, qslp->A.matval, &nsing, &singr, &singc);
	ok = (rval == 0 && nsing == 0);

CLEANUP:
	ILL_IFFREE (singr);
	ILL_IFFREE (singc);
	ILL_IFFREE (bcopy);
	mpq_ILLfactor_free_factor_work (&f);
	mpq_EGlpNumClearVar (f.fzero_tol);
	mpq_EGlpNumClearVar (f.szero_tol);
	mpq_EGlpNumClearVar (f.partial_tol);
	mpq_EGlpNumClearVar (f.maxelem_orig);
	mpq_EGlpNumClearVar (f.maxelem_factor);
	mpq_EGlpNumClearVar (f.maxelem_cur);
	mpq_EGlpNumClearVar (f.partial_cur);
	return ok;
}

// dump a single basis, rebuild from A
static int dump_one_basis (const char *dir, int k, int nrows_qs,
		const int *baz, mpq_ILLlpdata *qslp)
{
	char path[2048];
	snprintf (path, sizeof (path), "%s/basis_k%d_B.txt", dir, k);
	EGioFile_t *f = EGioOpen (path, "w");
	if (!f) return 1;

	int nnz = 0;
	for (int j = 0; j < nrows_qs; ++j) nnz += qslp->A.matcnt[baz[j]];

	EGioPrintf (f, "%d %d %d\n", nrows_qs, nrows_qs, nnz);
	char numbuf[8192];
	for (int j = 0; j < nrows_qs; ++j) {
		int col = baz[j];
		int beg = qslp->A.matbeg[col], cnt = qslp->A.matcnt[col];
		for (int k2 = 0; k2 < cnt; ++k2) {
			int row = qslp->A.matind[beg + k2];
			const mpq_t *v = &qslp->A.matval[beg + k2];
			gmp_snprintf (numbuf, sizeof (numbuf), "%Qd", *v);
			EGioPrintf (f, "%d %d %s\n", row, j, numbuf);
		}
	}

	EGioClose (f);
	return 0;
}

/* create the directory <problem>_Bases in the cwd and dump the constraint
 * matrix A and the last k basis matrices in there
 */
static void dump_basis_snapshots (mpq_QSdata *p_mpq, const char *input_fname)
{
	// make sure that we have at least one snapshot
	int n = basis_snapshot_count();
	int total_snaps = get_num_snapshots();
	if (n <= 0) {
		fprintf (stderr, "no basis snapshots captured; nothing to dump\n");
		return;
	}


	if (!p_mpq || !p_mpq->qslp) return;
	mpq_ILLlpdata *qslp = p_mpq->qslp;

	char prob[1024];
	derive_problem_name (input_fname, prob, sizeof (prob));
	char dir[2048];
	snprintf (dir, sizeof (dir), "%s_Bases", prob);

	if (mkdir (dir, 0755) != 0 && errno != EEXIST) {
		fprintf (stderr, "could not create directory %s: %s\n",
				dir, strerror (errno));
		return;
	}

	// meta.txt
	{
		char path[2048];
		snprintf (path, sizeof (path), "%s/meta.txt", dir);
		EGioFile_t *f = EGioOpen (path, "w");
		if (f) {
			EGioPrintf (f, "problem: %s\n", prob);
			EGioPrintf (f, "source_file: %s\n", input_fname);
			EGioPrintf (f, "nrows: %d\n", qslp->nrows);
			EGioPrintf (f, "ncols_total: %d\n",
					qslp->A.matcols);
			EGioPrintf (f, "nstruct: %d\n", qslp->nstruct);
			EGioPrintf (f, "total_snapshots_captured: %d\n", total_snaps);
			EGioClose (f);
		}
	}

	// write the constraint matrix A to <Problem>_Bases
	{
		char path[2048];
		snprintf (path, sizeof (path), "%s/A.txt", dir);
		write_sparse_mpq_matrix (path, qslp->A.matrows, qslp->A.matcols,
				qslp->A.matval, qslp->A.matbeg, qslp->A.matcnt, qslp->A.matind);
	}

	// write the rhs vector b into b.txt
	{
		char path[2048];
		snprintf (path, sizeof (path), "%s/b.txt", dir);
		write_dense_mpq_vector (path, qslp->nrows, qslp->rhs);
	}

	/* dump the most recent NDUMP exactly-nonsingular bases from the ring.
	 * Candidates that fail the exact factorization are skipped with a
	 * note: the tail of a float run's pivot trail can be exactly singular
	 * (that is typically WHY the solver escalated precision), and dumping
	 * those would just make the downstream factorization demos fall over. */
	#define NDUMP 5
	int nvalid = 0, nskipped = 0;
	for (int k = 0; k < n && nvalid < NDUMP; ++k) {
		int nrows_snap = 0;
		const int *baz = NULL;
		if (basis_snapshot_get (k, &nrows_snap, &baz) != 0) continue;
		if (nrows_snap != qslp->nrows) {
			fprintf (stderr,
				"basis snapshot k=%d has nrows=%d but qslp->nrows=%d; skipping\n",
				k, nrows_snap, qslp->nrows);
			continue;
		}
		if (!snapshot_exactly_nonsingular (qslp, baz)) {
			nskipped++;
			continue;
		}
		dump_one_basis (dir, nvalid, nrows_snap, baz, qslp);
		nvalid++;
	}

	// remove stale basis files from earlier runs beyond what we just wrote
	for (int k = nvalid; k < BASIS_SNAPSHOT_CAPACITY; ++k) {
		char path[2048];
		snprintf (path, sizeof (path), "%s/basis_k%d_B.txt", dir, k);
		remove (path);
	}

	// record the validation outcome next to the data
	{
		char path[2048];
		snprintf (path, sizeof (path), "%s/meta.txt", dir);
		EGioFile_t *f = EGioOpen (path, "a");
		if (f) {
			EGioPrintf (f, "ring_candidates: %d\n", n);
			EGioPrintf (f, "exactly_singular_skipped: %d\n", nskipped);
			EGioPrintf (f, "valid_bases_dumped: %d\n", nvalid);
			EGioClose (f);
		}
	}

	fprintf (stderr, "wrote %d exactly-nonsingular basis snapshot(s) + A to "
			"%s/ (%d of %d ring candidates were exactly singular and were "
			"skipped)\n", nvalid, dir, nskipped, n);
	if (nvalid < NDUMP) {
		fprintf (stderr, "note: fewer than %d valid bases existed in the "
				"last %d pivots; the float trajectory was exactly singular "
				"for the rest\n", NDUMP, n);
	}
	#undef NDUMP
}

/* ========================================================================= */
/** @brief the main thing! */
/* ========================================================================= */
int main (int ac, char **av)
{
	// clock start for timing purposes
    clock_t start = clock();
	int rval = 0,
	status = 0;
	mpq_QSdata *p_mpq = 0;
	QSbasis *basis = 0;
	ILLutil_timer timer_solve;
	ILLutil_timer timer_read;
	int ftype = 0;								/* 0 mps, 1 lp */
	mpq_t *y_mpq = 0,
	*x_mpq = 0;
	QSopt_ex_version();
	QSexactStart(); // AP: function in exact.c, calls EGlpNumStart in eg_lpnum.c

	/* parse arguments and initialize EGlpNum related things */
	rval = parseargs (ac, av);
	QSexact_set_precision (precision);
	if (rval)
		goto CLEANUP;
	if (writebasis)
	{
		// basis is pointer to space allocated for final basis, set to all 0s
		basis = EGsMalloc (QSbasis, 1);
		memset (basis, 0, sizeof (QSbasis));
	}

	/* just for the bell's and wistle */
	if (showversion)
	{
		char *buf = 0;
		buf = mpq_QSversion ();
		if (buf == 0)
		{
			ILL_CLEANUP;
		}
		else
		{
			printf ("%s\n", buf);
			mpq_QSfree ((void *) buf);
		}
	}

	/* get the file type */
	if (lpfile)
		ftype = 1;
	else
		get_ftype (fname, &ftype);

	// save file name to data sheet
	if (timing_enabled) {
		EGioFile_t *out = 0;
		out = EGioOpen ("time_precision_data", "a");
		EGioPrintf (out, "%s\n", fname);
		EGioClose (out);
	}

	// adds section header after fname is defined
	log_session_header(fname);

	/* read the mpq problem */
	ILLutil_init_timer (&timer_read, "SOLVER_READ_MPQ");
	ILLutil_start_timer (&timer_read);
	if (ftype == 1)
	{
		p_mpq = mpq_QSread_prob ((const char *) fname, "LP");
		if (p_mpq == 0)
		{
			fprintf (stderr, "Could not read lp file.\n");
			rval = 1;
			ILL_CLEANUP_IF (rval);
		}
	}
	else
	{
		p_mpq = mpq_QSread_prob ((const char *) fname, "MPS");
		if (p_mpq == 0)
		{
			fprintf (stderr, "Could not read mps file.\n");
			rval = 1;
			ILL_CLEANUP_IF (rval);
		}
	}

	/* and get the basis if needed */
	if (readbasis)
	{
		rval = mpq_QSread_and_load_basis (p_mpq, (const char *) readbasis);
		ILL_CLEANUP_IF (rval);
		if (basis)
			mpq_QSfree_basis (basis);
		basis = mpq_QSget_basis (p_mpq);
	}
	ILLutil_stop_timer (&timer_read, 1);
	/* set the readed flags */
	rval = mpq_QSset_param (p_mpq, QS_PARAM_SIMPLEX_DISPLAY, 1)
		|| mpq_QSset_param (p_mpq, QS_PARAM_PRIMAL_PRICING, pstrategy)
		|| mpq_QSset_param (p_mpq, QS_PARAM_DUAL_PRICING, dstrategy)
		|| mpq_QSset_param (p_mpq, QS_PARAM_SIMPLEX_SCALING, usescaling);
	ILL_CLEANUP_IF (rval);
	if (printsol)
	{
		x_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->ncols);
		y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
	}
	ILLutil_init_timer (&timer_solve, "SOLVER");
	ILLutil_start_timer (&timer_solve);
	basis_snapshot_init (p_mpq->qslp->nrows);
	if (pivotlog)
	{
		char prob[1024];
		derive_problem_name (fname, prob, sizeof (prob));
		rval = pivot_log_open (pivotlog, prob);
		ILL_CLEANUP_IF (rval);
	}
	rval = QSexact_solver (p_mpq, x_mpq, y_mpq, basis, simplexalgo, &status);
	pivot_log_close ();
	ILL_CLEANUP_IF (rval);
	ILLutil_stop_timer (&timer_solve, 1);
	if (dump_bases) dump_basis_snapshots (p_mpq, fname);
	if (printsol)
	{
		char out_f_name[1024];
		EGioFile_t *out_f;
		sprintf (out_f_name, "%s", solname);
		out_f = EGioOpen (out_f_name, "w");
		switch (status)
		{
		case QS_LP_OPTIMAL:
			EGioPrintf (out_f, "status = OPTIMAL\n");
			rval = QSexact_print_sol (p_mpq, out_f);
			CHECKRVALG(rval,CLEANUP);
			break;
		case QS_LP_INFEASIBLE:
			EGioPrintf (out_f, "status = INFEASIBLE\n");
			break;
		case QS_LP_UNBOUNDED:
			EGioPrintf (out_f, "status = UNBOUNDED\n");
			break;
		default:
			EGioPrintf (out_f, "status = UNDEFINED\n");
			break;
		}
		EGioClose (out_f);
	}
	/* ending */
CLEANUP:
	if (printsol) EGfree(solname);
	mpq_EGlpNumFreeArray (x_mpq);
	mpq_EGlpNumFreeArray (y_mpq);
	/* free the last allocated basis, and if we wanted to save it, do so */
	if (basis)
	{
		// AP: basis written to writebasis (file provided with -b argument)
		if (writebasis)
			rval = mpq_QSwrite_basis (p_mpq, 0, writebasis);
	}
	mpq_QSfree_basis (basis);
	mpq_QSfree_prob (p_mpq);
	QSexactClear();

	// testing for log file
    	clock_t end = clock();
    	double duration = (double)(end - start) / CLOCKS_PER_SEC;
   	log_timing("Total testing time ", duration);
	log_session_footer(fname);

	return rval;									/* main return */
}
