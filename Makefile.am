ACLOCAL_AMFLAGS = -I m4 --install
bin_PROGRAMS = miRA
miRA_SOURCES = src/main.c src/help.c src/cluster.c src/parse_sam.c src/errors.c src/vfold.c src/bed.c src/fasta.c
miRA_CFLAGS = -std=c11 -I src/ViennaRNA-2.1.8/H/
miRA_CPPFLAGS = -DDEBUG
miRA_LDFLAGS = -L src/ViennaRNA-2.1.8/lib/




EXTRA_PROGRAMS = miRAtest


miRAtest_SOURCES = test/main.c test/testerino.c test/test_cluster.c test/test_parse_sam.c test/test_bed_file_io.c src/errors.c src/parse_sam.c src/cluster.c src/vfold.c src/bed.c src/fasta.c test/test_fasta.c
miRAtest_CFLAGS = -std=c11

.PHONY: clean test

test:	miRAtest
	./miRAtest




EXTRA_DIST = m4/NOTES
