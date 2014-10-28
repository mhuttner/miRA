ACLOCAL_AMFLAGS = -I m4 --install
bin_PROGRAMS = miRA
miRA_SOURCES = src/main.c src/help.c src/cluster.c src/parse_sam.c src/errors.c
miRA_CFLAGS = -std=c11
miRA_CPPFLAGS = -DDEBUG

check_PROGRAMS = miRAtest
TESTS = $(check_PROGRAMS)

miRAtest_SOURCES = test/main.c
miRAtest_CFLAGS = -std=c11



EXTRA_DIST = m4/NOTES
