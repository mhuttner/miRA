ACLOCAL_AMFLAGS = -I m4 --install
bin_PROGRAMS = hello
hello_SOURCES = src/hello.c src/sum.c src/sum.h
hello_CPPFLAGS = -DDEBUG

TESTS = check_sum
check_PROGRAMS = check_sum
check_sum_SOURCES = test/check_sum.c src/sum.c src/sum.h
check_sum_CFLAGS = @CHECK_CFLAGS@
check_sum_LDADD = @CHECK_LIBS@

EXTRA_DIST = m4/NOTES
