AC_INIT([miRA], [0.1], [michael@mhuttner.com],
        [miRA], [http://mhuttner.com/])
AC_PREREQ([2.68])

AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([1.11 subdir-objects -Wall -Werror])
AC_CONFIG_FILES([Makefile])

AM_WITH_DMALLOC

AC_CHECK_LIB([m], [exp])


AC_HEADER_STDC
AC_HEADER_STDBOOL
AC_CHECK_HEADERS([malloc.h float.h limits.h stdlib.h string.h strings.h unistd.h unistd.h math.h])


AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor strdup strstr strchr strrchr strstr strtol strtoul pow rint sqrt erand48 memset memmove erand48])

AC_C_CONST
AC_TYPE_SIZE_T
AC_C_INLINE([])

AC_PROG_RANLIB
AM_PROG_AR
AC_PROG_CC
AM_PROG_CC_C_O
AX_PROG_JAVA
if test x$ac_cv_prog_java_works = xyes; then
	AC_DEFINE([HAVE_JAVA],[1],[Does Java work])
fi
AC_DEFINE_UNQUOTED([JAVA_VM_COMMAND],"$ac_cv_prog_JAVA",[Command to call java])



AC_OUTPUT
