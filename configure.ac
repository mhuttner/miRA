AC_INIT([hello], [0.01], [ss@sobering.ca],
        [hello], [http://tankapp.com/])
AC_PREREQ([2.68])

AC_CONFIG_SRCDIR([src/hello.c])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([1.11 subdir-objects -Wall -Werror])
AC_CONFIG_FILES([Makefile])

AC_PROG_CC
AM_PROG_CC_C_O

PKG_CHECK_MODULES([CHECK], [check >= 0.9.4])

AC_OUTPUT
