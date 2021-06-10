#	Function: AC_DW_GSL()
# 	In: GSL_CONFIG_PATH
#	Out: LIB_GSL,INC_GSL,LDD_GSL
AC_DEFUN([AC_DW_GSL],[
	AC_MSG_NOTICE([checking for gsl])
	#test if GSL_CONFIG given by user
	if test -z $GSL_CONFIG_PATH; then
		AC_MSG_NOTICE([GSL_CONFIG_PATH not specified -> try PATH])
		GSL_CONFIG_PATH=$PATH
	fi
	AC_MSG_NOTICE([GSL_CONFIG_PATH = $GSL_CONFIG_PATH])
	#test if GSL_CONFIG/gsl-config exists and set GSLCONFIG
 	AC_PATH_PROG(GSL_CONFIG, gsl-config, [], [$GSL_CONFIG_PATH])
	AC_MSG_NOTICE([gsl-config = $GSL_CONFIG])
	#test if gsl-config is executable
	if test -x "$GSL_CONFIG"; then
		#gsl-config is executable -> set flags
		AC_SUBST([LIB_GSL],`$GSL_CONFIG --libs`)
		AC_SUBST([INC_GSL],`$GSL_CONFIG --cflags`)
		AC_SUBST([LDD_GSL],["-Xlinker -R`$GSL_CONFIG --prefix`/lib"])
		AC_MSG_NOTICE([   INC_GSL = $INC_GSL])
		AC_MSG_NOTICE([   LIB_GSL = $LIB_GSL])
		AC_MSG_NOTICE([   LDD_GSL = $LDD_GSL])
	else
		#gsl-config is NOT executable
		AC_CHECK_LIB([m],[cos],[],[])
		AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],[AC_MSG_ERROR([Cannot find gslcblas library. Please set "GSL_CONFIG_PATH" or install gslcblas in common directorys.])])
		AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR([Cannot find gsl library. Please set "GSL_CONFIG_PATH" or install gsl in common directorys.])])
	fi
	AC_DEFINE_UNQUOTED([HAVE_LIBGSL],[1],[])
	AC_DEFINE_UNQUOTED([HAVE_LIBGSLCBLAS],[1],[])
	AC_DEFINE_UNQUOTED([HAVE_LIBM],[1],[])
#	AC_SUBST([GSL_INCDIR],`$GSL_CONFIG --prefix`/include)
#	AC_CHECK_HEADERS($GSL_INCDIR/gsl/gsl_interp2d.h,[],[AC_MSG_ERROR(Could not find 'gsl_interp2d.h'. Please make sure you have GSL version 2.1 or higher)])
])
