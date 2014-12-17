ACLOCAL_AMFLAGS = -I m4 --install
bin_PROGRAMS = miRA
miRA_SOURCES = src/main.c src/help.c src/cluster.c src/parse_sam.c src/errors.c src/vfold.c src/bed.c src/fasta.c
miRA_CFLAGS = -std=c11 -I src/Lfold/
miRA_CPPFLAGS = -DDEBUG
miRA_LDADD = libLfold.a




lib_LIBRARIES = libLfold.a


libLfold_a_SOURCES  = src/Lfold/fold_vars.c \
                    src/Lfold/energy_par.c \
                    src/Lfold/utils.c \
                    src/Lfold/params.c \
                    src/Lfold/Lfold.c \
                    src/Lfold/gquad.c

libRNA_a_CFLAGS = $(AM_CFLAGS)


noinst_HEADERS  = src/Lfold/intl11.h src/Lfold/intl11dH.h\
                  src/Lfold/intl21.h src/Lfold/intl21dH.h\
                  src/Lfold/intl22.h src/Lfold/intl22dH.h
                  






EXTRA_PROGRAMS = miRAtest


miRAtest_SOURCES = test/main.c test/testerino.c test/test_cluster.c test/test_parse_sam.c test/test_bed_file_io.c src/errors.c src/parse_sam.c src/cluster.c src/vfold.c src/bed.c src/fasta.c test/test_fasta.c test/test_vfold.c
miRAtest_CFLAGS = -std=c11 -I src/ViennaRNA-2.1.8/H/

.PHONY: clean test

test:	miRAtest
	./miRAtest




EXTRA_DIST = m4/NOTES src/Lfold/model_avg.inc src/Lfold/model_sd.inc
