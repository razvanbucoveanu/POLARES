pkginclude_HEADERS = POLARES.h ReadConfigFile.h const.h interpolation.h IO_classes.h form_factors.h \
parameters.h cross_sections.h integrands.h virtual_corrections.h cuba_param.h gsl_rand.h melem.h \
scalar_integrals.h config.h gamma_loop.h melem_pol.h

if LOOPTOOLS
pkginclude_HEADERS +=  looptools_interface.h
endif

#Library

lib_LTLIBRARIES = libPOLARES.la
libPOLARES_la_SOURCES = POLARES.cpp ReadConfigFile.cpp interpolation.cpp form_factors.cpp \
parameters.cpp cross_sections.cpp integrands.cpp virtual_corrections.cpp gsl_rand.cpp \
cuba_param.cpp melem.cpp scalar_integrals.cpp gamma_loop.cpp melem_pol.cpp

if LOOPTOOLS
libPOLARES_la_SOURCES +=  looptools_interface.cpp
endif

libPOLARES_la_CPPFLAGS = $(INC_GSL) -I$(top_srcdir)/$(cubadir)
libPOLARES_la_LIBADD =  $(LIB_GSL) $(top_builddir)/$(cubadir)/libcuba.la

if LOOPTOOLS
#code needed if LoopTools is not distributed with POLARES
#libPOLARES_la_CPPFLAGS += $(INC_LOOPTOOLS)
#libPOLARES_la_LIBADD += $(LIB_LOOPTOOLS)
#libPOLARES_la_LDFLAGS = $(LDD_LOOPTOOLS)
#code needed if LoopTools is distributed with POLARES
libPOLARES_la_CPPFLAGS += -I$(top_srcdir)/$(looptoolsdir)/build
libPOLARES_la_LIBADD += -L$(top_builddir)/$(looptoolsdir)/build -looptools -lm -lgfortran
libPOLARES_la_LDFLAGS = -static -fPIC
endif

#bin_PROGRAMS = POLARES
#noinst_PROGRAMS = POLARES

#POLARES_SOURCES = POLARES_main.cpp

#POLARES_CPPFLAGS = -I$(top_srcdir)/src
#POLARES_LDADD = libPOLARES.la
