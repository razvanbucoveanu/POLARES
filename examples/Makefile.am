EXAMPLES = main_example multiple_random_E_test sigma_diff_Omega_l_test shiftQ2

EXTRA_PROGRAMS = $(EXAMPLES)

multiple_random_E_test_SOURCES = multiple_random_E_test.cpp
sigma_diff_Omega_l_test_SOURCES = sigma_diff_Omega_l_test.cpp
main_example_SOURCES = main_example.cpp
shiftQ2_SOURCES = shiftQ2.cpp

multiple_random_E_test_CPPFLAGS = -I$(top_srcdir)/src
multiple_random_E_test_LDADD = -L$(top_srcdir)/src -lPOLARES
sigma_diff_Omega_l_test_CPPFLAGS = -I$(top_srcdir)/src
sigma_diff_Omega_l_test_LDADD = -L$(top_srcdir)/src -lPOLARES
main_example_CPPFLAGS = -I$(top_srcdir)/src
main_example_LDADD = -L$(top_srcdir)/src -lPOLARES
#main_example_LDFLAGS = -static -fPIC
shiftQ2_CPPFLAGS = -I$(top_srcdir)/src
shiftQ2_LDADD = -L$(top_srcdir)/src -lPOLARES


EXTRA_DIST = POLARES.out POLARES_basic.out Makefile_example

.PHONY: examples
examples: $(EXAMPLES)
