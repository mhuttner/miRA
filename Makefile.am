ACLOCAL_AMFLAGS = -I m4 --install
bin_PROGRAMS = miRA
miRA_SOURCES = src/main.c src/help.c src/cluster.c src/parse_sam.c src/errors.c src/vfold.c src/bed.c src/fasta.c src/util.c src/structure_evaluation.c src/candidates.c src/coverage.c src/reporting.c
miRA_CFLAGS = -std=c11
miRA_CPPFLAGS = -DDEBUG
miRA_LDADD = libLfold.a




lib_LIBRARIES = libLfold.a


libLfold_a_SOURCES  = src/Lfold/fold_vars.c \
                    src/Lfold/energy_par.c \
                    src/Lfold/utils.c \
                    src/Lfold/params.c \
                    src/Lfold/Lfold.c \
                    src/Lfold/gquad.c \
                    src/Lfold/fold.c

libRNA_a_CFLAGS = -std=c11 


noinst_HEADERS  = src/Lfold/intl11.h src/Lfold/intl11dH.h\
                  src/Lfold/intl21.h src/Lfold/intl21dH.h\
                  src/Lfold/intl22.h src/Lfold/intl22dH.h
                  






EXTRA_PROGRAMS = miRAtest


miRAtest_SOURCES = test/main.c test/testerino.c test/test_cluster.c test/test_parse_sam.c test/test_bed_file_io.c src/errors.c src/parse_sam.c src/cluster.c src/vfold.c src/bed.c src/fasta.c test/test_fasta.c test/test_vfold.c src/util.c test/test_util.c src/structure_evaluation.c src/candidates.c src/coverage.c
miRAtest_CFLAGS = -std=c11 
miRAtest_LDADD = libLfold.a

.PHONY: clean test

test:	miRAtest
	./miRAtest




EXTRA_DIST = m4/NOTES src/Lfold/model_avg.inc src/Lfold/model_sd.inc src/Lfold/circfold.inc
