# needed for boost.m4
#ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = POLARES.in POLARES_basic.in
SUBDIRS = $(looptoolsdir) $(cubadir) src examples share doc

bin_SCRIPTS = POLARES-config
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA= POLARES.pc

.PHONY: examples
examples:
	make examples -C examples
