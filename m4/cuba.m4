#	Function: AC_DW_CUBA()
# 	In: CUBA_CONFIG_PATH
#	Out: LIB_CUBA,INC_CUBA,LDD_CUBA
m4_define([ac_cubadir], [src/cuba-4.2])

AC_DEFUN([AC_DW_CUBA],[
	AC_MSG_NOTICE([Configure cuba])
	#AC_CHECK_LIB([vegas],[cubaruninit],[],[AC_MSG_ERROR([Cannot find cuba library. Please set "CUBA_PREFIX" or install cuba in common directorys.])],[])
	AC_CHECK_DECL([LLONG_MAX],[AC_MSG_NOTICE([DECL found llmax])],[AC_MSG_NOTICE([DECL NOT found llmax])],[#include <limits.h>])
	AC_CHECK_DECL([INT_MAX],[AC_MSG_NOTICE([DECL found i_max])],[AC_MSG_NOTICE([DECL NOT found i_max])],[#include <limits.h>])

	AC_SUBST([cubadir],ac_cubadir)
	AC_CONFIG_SUBDIRS(ac_cubadir)
])