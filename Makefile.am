ACLOCAL_AMFLAGS = -I m4 --install
bin_PROGRAMS = miRA
miRA_SOURCES = src/main.c src/help.c src/help.h src/cluster.c src/cluster.h src/parse_sam.c src/parse_sam.h src/errors.c src/errors.h
miRA_CPPFLAGS = -DDEBUG

TESTS = check_sum check_parse_sam
check_PROGRAMS = check_sum check_parse_sam
check_sum_SOURCES = test/check_sum.c
check_sum_CFLAGS = @CHECK_CFLAGS@
check_sum_LDADD = @CHECK_LIBS@

check_parse_sam_SOURCES = test/check_parse_sam.c src/parse_sam.c src/parse_sum.h src/errors.h
check_parse_sam_CFLAGS = @CHECK_CFLAGS@
check_parse_sam_LDADD = @CHECK_LIBS@

EXTRA_DIST = m4/NOTES
