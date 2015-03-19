ACLOCAL_AMFLAGS = -I m4 --install
bin_PROGRAMS = miRA
miRAdir = src
miRA_SOURCES = src/main.c src/help.c src/cluster.c src/parse_sam.c src/errors.c src/vfold.c src/bed.c src/fasta.c src/util.c src/structure_evaluation.c src/candidates.c src/coverage.c src/reporting.c src/full.c
miRA_HEADERS = src/help.h src/cluster.h src/parse_sam.h src/errors.h src/vfold.h src/bed.h src/fasta.h src/util.h src/structure_evaluation.h src/candidates.h src/coverage.h src/reporting.h src/full.h src/defs.h src/uthash.h
miRA_CFLAGS = -std=c99 $(OPENMP_CFLAGS)
miRA_CPPFLAGS = -DDEBUG
miRA_LDADD = libLfold.a




lib_LIBRARIES = libLfold.a

libLfold_adir = src/Lfold
libLfold_a_SOURCES  = src/Lfold/fold_vars.c \
                    src/Lfold/energy_par.c \
                    src/Lfold/utils.c \
                    src/Lfold/params.c \
                    src/Lfold/Lfold.c \
                    src/Lfold/gquad.c \
                    src/Lfold/fold.c

libLfold_a_HEADERS = src/Lfold/Lfold.h \
src/Lfold/energy_par.h \
src/Lfold/loop_energies.h \
src/Lfold/aln_util.h \
src/Lfold/fold.h \
src/Lfold/pair_mat.h \
src/Lfold/config_old.h \
src/Lfold/fold_vars.h \
src/Lfold/params.h \
src/Lfold/data_structures.h \
src/Lfold/gquad.h \
src/Lfold/utils.h \
src/Lfold/energy_const.h


libLfold_a_CFLAGS = -std=c99 $(OPENMP_CFLAGS)


noinst_HEADERS  = src/Lfold/intl11.h src/Lfold/intl11dH.h\
                  src/Lfold/intl21.h src/Lfold/intl21dH.h\
                  src/Lfold/intl22.h src/Lfold/intl22dH.h
                  






EXTRA_PROGRAMS = miRAtest

miRAtestdir = test
miRAtest_SOURCES = test/main.c test/testerino.c test/test_cluster.c test/test_parse_sam.c test/test_bed_file_io.c src/errors.c src/parse_sam.c src/cluster.c src/vfold.c src/bed.c src/fasta.c test/test_fasta.c test/test_vfold.c src/util.c test/test_util.c src/structure_evaluation.c src/candidates.c src/coverage.c src/reporting.c src/full.c
miRAtest_HEADERS = test/testerino.h test/test_cluster.h test/test_parse_sam.h test/test_bed_file_io.h src/errors.h src/parse_sam.h src/cluster.h src/vfold.h src/bed.h src/fasta.h test/test_fasta.h test/test_vfold.h src/util.h test/test_util.h src/structure_evaluation.h src/candidates.h src/coverage.h src/reporting.h src/full.h src/defs.h src/uthash.h
miRAtest_CFLAGS = -std=c99
miRAtest_LDADD = libLfold.a

.PHONY: clean test

test:	miRAtest
	./miRAtest




EXTRA_DIST = m4/NOTES src/Lfold/model_avg.inc src/Lfold/model_sd.inc src/Lfold/circfold.inc VARNAv3-91.jar README.md LICENSE src/Lfold/COPYING test/data/test.config test/data/test.fasta example/sample_reads.sam example/sample_sequence.fasta example/sample_configuration.config example/sample_output
