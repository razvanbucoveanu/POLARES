#	Function: AC_DW_LOOPTOOLS()
# 	In: LOOPTOOLS_PATH
#	Out: LIB_LOOPTOOLS,INC_LOOPTOOLS,LDD_LOOPTOOLS
m4_define([ac_looptools], [src/LoopTools-2.12])

AC_DEFUN([AC_DW_LOOPTOOLS],[
#code needed if LoopTools is not distributed with RC_PEPS
#	AC_MSG_NOTICE([checking for looptools])
#test if LOOPTOOLS_PATH given by user
#	if test -z $LOOPTOOLS_PATH; then
#		AC_MSG_NOTICE([LOOPTOOLS_PATH not specified -> try PATH])
#		LOOPTOOLS_PATH=$PATH
#	fi
#	AC_MSG_NOTICE([LOOPTOOLS_PATH = $LOOPTOOLS_PATH])
#	AC_CHECK_FILE($LOOPTOOLS_PATH/lib64/libooptools.a,[],[AC_MSG_ERROR(Could not find 'libooptools.a'. Please check looptools path!)])

#	AC_MSG_NOTICE([***** Configure LoopTools *****])

#	AC_SUBST([LIB_LOOPTOOLS],"-L$LOOPTOOLS_PATH/lib64 -looptools -lm -lgfortran")
#	AC_SUBST([INC_LOOPTOOLS],"-I$LOOPTOOLS_PATH/include")
#	AC_SUBST([LDD_LOOPTOOLS],"-static")
#	AC_MSG_NOTICE([   INC_LOOPTOOLS = $INC_LOOPTOOLS])
#	AC_MSG_NOTICE([   LIB_LOOPTOOLS = $LIB_LOOPTOOLS])
#	AC_MSG_NOTICE([   LDD_LOOPTOOLS = $LDD_LOOPTOOLS])

#code needed if LoopTools is distributed with RC_PEPS
	AC_MSG_NOTICE([***** Configure LoopTools *****])
	AC_SUBST([looptoolsdir],ac_looptools)
	AC_CONFIG_SUBDIRS(ac_looptools)

])
