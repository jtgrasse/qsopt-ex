#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "exact.h"
#include "pivot_log.h"

int pivot_log_active = 0;

static FILE *out = NULL;
// internal -> canonical column map of the LP that last asked for it
static const void *map_owner = NULL;
static int *map = NULL;
static int map_ncols = 0;
// the basis as last written to the log, as a membership array over the
// canonical columns; nbasic < 0 means nothing comparable has been written
static char *in_basis = NULL;
static int in_basis_size = 0;
static int nbasic = -1;
static int nbases_written = 0;
static int nupdates_written = 0;
static char reason[512];
static int suspended = 0;

int pivot_log_open (const char *path, const char *probname)
{
	pivot_log_close ();
	out = fopen (path, "w");
	if (!out)
	{
		fprintf (stderr, "could not open pivot log %s for writing\n", path);
		return 1;
	}
	// rationals and basis lines can be long, buffer generously
	setvbuf (out, NULL, _IOFBF, 1 << 20);
	pivot_log_active = 1;
	nbasic = -1;
	nbases_written = 0;
	nupdates_written = 0;
	reason[0] = 0;

	fprintf (out, "%% problem: %s\n", probname);
	fprintf (out, "%% QSopt_ex pivot log.  Lines starting with %% are comments.\n");
	fprintf (out, "%% All indices are 0-based.  Values are exact rationals written as p/q\n");
	fprintf (out, "%% (or p when integral).  Entries not listed in a sparse section are 0.\n");
	return 0;
}

void pivot_log_close (void)
{
	if (out)
		fclose (out);
	out = NULL;
	pivot_log_active = 0;
	suspended = 0;
	free (map);
	map = NULL;
	map_owner = NULL;
	map_ncols = 0;
	free (in_basis);
	in_basis = NULL;
	in_basis_size = 0;
	nbasic = -1;
}

void pivot_log_suspend (void)
{
	if (out && !suspended++)
		pivot_log_active = 0;
}

void pivot_log_resume (void)
{
	if (out && suspended && !--suspended)
		pivot_log_active = 1;
}

void pivot_log_comment (const char *fmt, ...)
{
	va_list va;
	if (!pivot_log_active) return;
	fputs ("% ", out);
	va_start (va, fmt);
	vfprintf (out, fmt, va);
	va_end (va);
	fputc ('\n', out);
}

void pivot_log_set_reason (const char *fmt, ...)
{
	va_list va;
	if (!pivot_log_active) return;
	va_start (va, fmt);
	vsnprintf (reason, sizeof (reason), fmt, va);
	va_end (va);
}

/* sparse vector section: "name\nlen nnz\n" then "k v" for each nonzero.
 * Values at or beyond +-ILL_MAXDOUBLE are QSopt_ex's infinity and are
 * written as inf / -inf. col maps canonical index to internal index. */
static void write_vector (const char *title, const char *name, int len,
		const mpq_t *v, const int *col)
{
	int k, nnz = 0;
	for (k = 0; k < len; k++)
		if (mpq_sgn (v[col ? col[k] : k]) != 0) nnz++;
	fprintf (out, "\n%% %s\n%s\n%d %d\n", title, name, len, nnz);
	for (k = 0; k < len; k++)
	{
		const mpq_t *x = &v[col ? col[k] : k];
		if (mpq_sgn (*x) == 0) continue;
		if (mpq_cmp (*x, mpq_ILL_MAXDOUBLE) >= 0)
			fprintf (out, "%d inf\n", k);
		else if (mpq_cmp (*x, mpq_ILL_MINDOUBLE) <= 0)
			fprintf (out, "%d -inf\n", k);
		else
			gmp_fprintf (out, "%d %Qd\n", k, *x);
	}
}

void pivot_log_problem (const mpq_ILLlpdata *qslp, const char *what)
{
	int i, j, k, nnz = 0;
	const int nstruct = qslp->nstruct, nrows = qslp->nrows, ncols = qslp->ncols;
	const mpq_ILLmatrix *A = &qslp->A;
	int *col;

	if (!pivot_log_active) return;
	// canonical column -> internal column
	col = (int *) malloc (sizeof (int) * (ncols > 0 ? ncols : 1));
	if (!col) return;
	for (j = 0; j < nstruct; j++) col[j] = qslp->structmap[j];
	for (i = 0; i < nrows; i++) col[nstruct + i] = qslp->rowmap[i];

	for (k = 0; k < ncols; k++)
		for (i = 0; i < A->matcnt[col[k]]; i++)
			if (mpq_sgn (A->matval[A->matbeg[col[k]] + i]) != 0) nnz++;

	fprintf (out, "\n%% %s\n", what);
	fprintf (out, "%% %d rows, %d columns: columns 0..%d are structural variables and\n",
			nrows, ncols, nstruct - 1);
	fprintf (out, "%% column %d+i is the logical variable of row i, so Ax = b holds with\n",
			nstruct);
	fprintf (out, "%% equality in every row.  Logicals: slack (+1, [0,inf)) for <= rows,\n");
	fprintf (out, "%% surplus (-1, [0,inf)) for >= rows, artificial (+1, fixed at 0) for\n");
	fprintf (out, "%% = rows, range (-1, [0,range]) for ranged rows.\n");
	fprintf (out, "%% Objective sense: %s c^T x subject to Ax = b, l <= x <= u.\n",
			qslp->objsense == QS_MAX ? "maximize" : "minimize");

	fprintf (out, "\n%% Constraint matrix\nA\n%d %d %d\n", nrows, ncols, nnz);
	for (k = 0; k < ncols; k++)
	{
		const int beg = A->matbeg[col[k]];
		for (i = 0; i < A->matcnt[col[k]]; i++)
		{
			if (mpq_sgn (A->matval[beg + i]) == 0) continue;
			gmp_fprintf (out, "%d %d %Qd\n", A->matind[beg + i], k,
					A->matval[beg + i]);
		}
	}
	write_vector ("Constraint RHS", "b", nrows, qslp->rhs, NULL);
	write_vector ("Objective coefficients", "c", ncols, qslp->obj, col);
	write_vector ("Variable lower bounds", "l", ncols, qslp->lower, col);
	write_vector ("Variable upper bounds", "u", ncols, qslp->upper, col);
	fputc ('\n', out);
	free (col);
	// a new problem invalidates the logged basis
	nbasic = -1;
}

void pivot_log_solution (int ncols, int nrows, const mpq_t *x, const mpq_t *y,
		const mpq_t *d, const mpq_t z)
{
	if (!pivot_log_active) return;
	fprintf (out, "\n%% Final exact rational solution: Ax = b, l <= x <= u, "
			"A^T y + d = c, z = c^T x = b^T y + d^T x\n");
	write_vector ("Primal solution (all columns, logicals included)", "x",
			ncols, x, NULL);
	write_vector ("Dual solution", "y", nrows, y, NULL);
	write_vector ("Reduced costs", "d", ncols, d, NULL);
	gmp_fprintf (out, "\n%% Objective value\nz\n%Qd\n", z);
}

int *pivot_log_colmap (const void *owner, int ncols, int *fill)
{
	if (owner != map_owner || ncols != map_ncols)
		*fill = 1;
	if (ncols > map_ncols || !map)
	{
		int *nm = (int *) realloc (map, sizeof (int) * (ncols > 0 ? ncols : 1));
		if (!nm)
		{
			pivot_log_comment ("out of memory, pivot log stopped");
			pivot_log_close ();
			return NULL;
		}
		map = nm;
	}
	if (ncols > in_basis_size)
	{
		char *nb = (char *) realloc (in_basis, (size_t) ncols);
		if (!nb)
		{
			pivot_log_comment ("out of memory, pivot log stopped");
			pivot_log_close ();
			return NULL;
		}
		memset (nb + in_basis_size, 0, (size_t) (ncols - in_basis_size));
		in_basis = nb;
		in_basis_size = ncols;
	}
	// indices of an LP of a different size are not comparable
	if (map_ncols && ncols != map_ncols)
		nbasic = -1;
	map_owner = owner;
	map_ncols = ncols;
	return map;
}

static void write_basis (int nrows, const int *baz)
{
	int i;
	if (reason[0])
		fprintf (out, "\n%% %s\n", reason);
	else
		fputc ('\n', out);
	fprintf (out, "%% %s\nB\n", nbases_written ? "New basis" : "Starting basis");
	memset (in_basis, 0, (size_t) in_basis_size);
	for (i = 0; i < nrows; i++)
	{
		const int k = map[baz[i]];
		fprintf (out, i ? " %d" : "%d", k);
		in_basis[k] = 1;
	}
	fputs ("\n\n", out);
	nbasic = nrows;
	nbases_written++;
	reason[0] = 0;
}

void pivot_log_sync (int nrows, const int *baz)
{
	int i;
	if (!pivot_log_active || !map) return;
	if (nbasic == nrows)
	{
		for (i = 0; i < nrows && in_basis[map[baz[i]]]; i++) ;
		if (i == nrows)
		{
			reason[0] = 0;
			return;
		}
	}
	write_basis (nrows, baz);
}

void pivot_log_pivot (int nrows, const int *baz, int lvar, int evar)
{
	int l, e;
	if (!pivot_log_active || !map) return;
	l = map[lvar];
	e = map[evar];
	if (nbasic != nrows || !in_basis[l] || in_basis[e])
	{
		// the basis was changed somewhere we do not hook; resynchronize
		if (!reason[0])
			pivot_log_set_reason ("basis changed outside of a logged pivot");
		write_basis (nrows, baz);
		return;
	}
	// label the columns of the first update in the file
	fprintf (out, nupdates_written++ ? "%% Update\nU\n%d %d\n\n" :
			"%% Update\nU\n%d %d     %%leaving entering\n\n", l, e);
	in_basis[l] = 0;
	in_basis[e] = 1;
}
