# Prelude.
AC_PREREQ(2.61)
AC_INIT([POLARES], [1.0], [rabucove@uni-mainz.de],[POLARES])
#Variables:
#	User:
#		GSL_CONFIG_PATH
#       LOOPTOOLS_PATH
#	Intern:
#		LIB_GSL
#		INC_GSL
#		LDD_GSL
#		CUBA_INC_DIR
#		CUBA_LIB_DIR
#		LIB_CUBA
#		INC_CUBA
#		LDD_CUBA
#		LIB_LOOPTOOLS
#		INC_LOOPTOOLS
#		LDD_LOOPTOOLS
AC_ARG_VAR(GSL_CONFIG_PATH,Optional: Path to gsl-config.)
#AC_ARG_VAR(LOOPTOOLS_PATH,Optional: Path to looptools.)
################################################
AC_CONFIG_MACRO_DIR([m4])
################################################
AC_MSG_NOTICE([***** Check for programs  *****])
AM_INIT_AUTOMAKE([-Wall foreign])
AM_PROG_AR([])
#############
AC_MSG_NOTICE([** LIBTOOL])
AC_PROG_LIBTOOL
#############
AC_MSG_NOTICE([** CXX])
AC_LANG(C++)
AC_PROG_CXX
AC_PROG_FC
AC_MSG_NOTICE([CXXFLAGS = $CXXFLAGS])
AC_MSG_NOTICE([LDFLAGS = $LDFLAGS])
################################################
AC_MSG_NOTICE([***** Check for librarys  *****])
#############	gsl
AC_DW_GSL
################################################
AC_MSG_NOTICE([***** Check for headers  *****])
AC_CHECK_HEADERS([sys/time.h])
################################################
AC_MSG_NOTICE([***** Check for  typedefs, structures, and compiler characteristics.  *****])
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME
################################################
AC_MSG_NOTICE([***** Checks for library functions.  *****])
AC_FUNC_ERROR_AT_LINE
AC_FUNC_SELECT_ARGTYPES
AC_CHECK_FUNCS([gettimeofday pow putenv select sqrt])
################################################
AC_MSG_NOTICE([***** Configure cuba *****])
AC_DW_CUBA
################################################
# Choose whether the interface to the LoopTools library should be compiled
AC_ARG_WITH([looptools],
   AS_HELP_STRING([--with-looptools],
                  [include the LoopTools library]),
                  [],
                  [with_looptools=no])
AS_IF([test "x$with_looptools" == xyes],
        AC_DEFINE([POLARES_USE_LOOPTOOLS], [1], [Define if you have LoopTools.])
        [looptools=true],
       [test "x$with_looptools" != xno],
        AC_DEFINE([POLARES_USE_LOOPTOOLS], [1], [Define if you have LoopTools.])
        [looptools=true],
        [looptools=false])
AM_CONDITIONAL(LOOPTOOLS, test "x$looptools" == xtrue)
AM_COND_IF(LOOPTOOLS,AC_DW_LOOPTOOLS)
################################################
#Determine CPU arhitechture (needed for looptools)
#AC_CANONICAL_HOST
################################################
AC_MSG_NOTICE([***** Set Output files *****])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([	Makefile \
			examples/Makefile \
			share/Makefile \
			doc/Makefile \
			src/Makefile \
		])
AC_CONFIG_FILES([POLARES-config])
AC_CONFIG_FILES([POLARES.pc])
AC_OUTPUT
