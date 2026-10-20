#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([QSopt_ex],[2.5.10.3],[https://github.com/jonls/qsopt-ex/issues])
AC_CONFIG_SRCDIR([qsopt_ex/QSopt_ex.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_REQUIRE_AUX_FILE([tap-driver.sh])
AM_INIT_AUTOMAKE([foreign subdir-objects dist-xz])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_USE_SYSTEM_EXTENSIONS

### Versioning
# This is the internal library version information used by libtool.
# IT DOES NOT FOLLOW THE NORMAL RELEASE VERSIONING!
# Instead follow these rules from the libtool documentation:
# 1. Update the version information only immediately before a public release of your software.
# 2. If the library source code has changed at all since the last update, then increment revision.
#    ("c:r:a" becomes "c:r+1:a")
# 3. If any interfaces have been added, removed, or changed since the last update, increment current,
#    and set revision to 0. ("c+1:0:a")
# 4. If any interfaces have been added since the last public release, then increment age.
# 5. If any interfaces have been removed since the last public release, then set age to 0.

m4_define([lt_current],[3])
m4_define([lt_revision],[0])
m4_define([lt_age],[1])
m4_define([lt_version_info],[lt_current:lt_revision:lt_age])
AC_SUBST([LT_VERSION_INFO], [lt_version_info])

### Checks for programs.
AC_PROG_CC
AC_PROG_SED
LT_INIT


### Check for debug mode
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug], [enable debug mode @<:@default=yes@:>@])],
	[enable_debug=$enable])
AS_IF([test "x$enable_debug" != xno], [
	AC_SUBST([ENABLE_DEBUG],[1])
	AC_DEFINE([DEBUG], [1], [Enable debug mode])
], [
  AC_SUBST([ENABLE_DEBUG],[0])
	AC_DEFINE([DEBUG], [0], [Disable debug mode])
])

### Check for verbosity level
AC_ARG_ENABLE([verbose-level],
	[AS_HELP_STRING([--enable-verbose-level], [fix verbose-level for compile-time routines @<:@default=100@:>@])],
	[],[enable_verbose_level=100])
AC_DEFINE_UNQUOTED([VERBOSE_LEVEL],[$enable_verbose_level],[Compilation-time verbose level])

### Check for libgmp
AC_ARG_VAR([GMP_CFLAGS], [C compiler flags for GMP, overriding automatic detection])
AC_ARG_VAR([GMP_LIBS], [linker flags for GMP, overriding automatic detection])
AS_IF([test -z "$GMP_CFLAGS" -a -z "$GMP_LIBS"], [
	AC_CHECK_LIB([gmp], [__gmpz_init], [
		GMP_CFLAGS=""
		GMP_LIBS="-lgmp"
	], [AC_MSG_ERROR([GNU MP not found, see https://gmplib.org/])])
])

### Check for libz
AC_CHECK_LIB([z], [gzopen], [], AC_MSG_NOTICE([compiling without libz]))

### Check for libbz2
AC_CHECK_LIB([bz2], [BZ2_bzopen], [], AC_MSG_NOTICE([compiling without libbz2]))

### Check for math library
AC_SEARCH_LIBS([pow], [m], [], [
	AC_MSG_ERROR([unable to find library with math functions (pow)])])

### Check for pthread library
AC_CHECK_LIB([pthread], [pthread_mutex_lock])

### Checks for header files.
m4_warn([obsolete],
[The preprocessor macro 'STDC_HEADERS' is obsolete.
  Except in unusual embedded environments, you can safely include all
  C89 headers unconditionally.])dnl
# Autoupdate added the next two lines to ensure that your configure
# script's behavior did not change.  They are probably safe to remove.
AC_CHECK_INCLUDES_DEFAULT
AC_PROG_EGREP

m4_warn([obsolete],
[Update your code to rely only on HAVE_SYS_TIME_H,
then remove this warning and the obsolete code below it.
All current systems provide time.h; it need not be checked for.
Not all systems provide sys/time.h, but those that do, all allow
you to include it and time.h simultaneously.])dnl
AC_CHECK_HEADERS_ONCE([sys/time.h])
# Obsolete code to be removed.
if test $ac_cv_header_sys_time_h = yes; then
  AC_DEFINE([TIME_WITH_SYS_TIME],[1],[Define to 1 if you can safely include both <sys/time.h>
	     and <time.h>.  This macro is obsolete.])
fi
# End of obsolete code.

AC_CHECK_HEADERS([errno.h float.h getopt.h inttypes.h limits.h math.h setjmp.h signal.h stdarg.h stdint.h stdio.h string.h sys/param.h sys/resource.h sys/time.h sys/times.h sys/utsname.h unistd.h])

### Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_TYPEOF
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
m4_warn([obsolete],
[your code may safely assume C89 semantics that RETSIGTYPE is void.
Remove this warning and the 'AC_CACHE_CHECK' when you adjust the code.])dnl
AC_CACHE_CHECK([return type of signal handlers],[ac_cv_type_signal],[AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([#include <sys/types.h>
#include <signal.h>
],
		 [return *(signal (0, 0)) (0) == 1;])],
		   [ac_cv_type_signal=int],
		   [ac_cv_type_signal=void])])
AC_DEFINE_UNQUOTED([RETSIGTYPE],[$ac_cv_type_signal],[Define as the return type of signal handlers
		    ('int' or 'void').])

AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

### Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor memset pow sqrt strdup strerror uname posix_memalign sleep getrusage times clock sigaction signal])

### Print output files
AC_CONFIG_FILES([
	Makefile
])
AC_OUTPUT


echo "
 $PACKAGE_NAME $VERSION

    prefix:		${prefix}
    compiler:		${CC}
    cflags:		${CFLAGS}
    ldflags:		${LDFLAGS}
"
