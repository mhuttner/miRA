ACLOCAL_AMFLAGS = -I m4 --install
bin_PROGRAMS = miRA
miRA_SOURCES = src/main.c src/help.c src/cluster.c src/parse_sam.c src/errors.c
miRA_CFLAGS = -std=c11
miRA_CPPFLAGS = -DDEBUG

noinst_PROGRAMS = miRAtest


miRAtest_SOURCES = test/main.c test/testerino.c test/test_cluster.c test/test_parse_sam.c src/errors.c src/parse_sam.c src/cluster.c
miRAtest_CFLAGS = -std=c11

.PHONY: clean test

test:	miRAtest
	./miRAtest




EXTRA_DIST = m4/NOTES
