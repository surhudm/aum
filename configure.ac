AC_INIT([aum], [1.01], [surhud@kicp.uchicago.edu])
AM_INIT_AUTOMAKE([gnu])

AC_CONFIG_SRCDIR([src/cosmology.cpp])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Make sure functions are checked with C++ not C
AC_LANG(C++)

# Checks for libraries.
AC_CHECK_LIB([gslcblas], [cblas_sdot],[],[AC_MSG_ERROR([Error! You need to have gsl cblas library. If you have libgslcblas installed in a non-standard location, please add the path to LDFLAGS])])
AC_CHECK_LIB([gsl], [gsl_linalg_cholesky_invert],[],[AC_MSG_ERROR([Error! You need to have gsl library. If you have libgsl installed in a non-standard location, please add the path to LDFLAGS ])])
AX_PATH_GSL(1.14,[],[AC_MSG_ERROR([Error! You need to have gsl library (version >=1.14). If you have libgsl installed in a non-standard location, please add the path to LDFLAGS ])])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

#Python swig related
AC_DISABLE_STATIC
AC_PROG_LIBTOOL
AM_PATH_PYTHON
AX_SWIG_MULTI_MODULE_SUPPORT
AX_SWIG_PYTHON

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt],[],[AC_MSG_ERROR([Error! Could not find pow and sqrt functions in cmath library])])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 pyext/Makefile])

#AC_OUTPUT(Makefile src/Makefile pyext/Makefile)
AC_OUTPUT
