#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([bwa], [0.7.5a], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([QSufSort.c])
AM_INIT_AUTOMAKE([foreign])
#AC_CONFIG_HEADERS([config.h])
# Checks for programs.
AC_PROG_CC
AC_PROG_CC_STDC

AC_PROG_RANLIB
# Checks for libraries.
# FIXME: Replace `main' with a function in `-lbwa':
AC_CHECK_LIB([bwa], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h fcntl.h limits.h netdb.h stddef.h stdint.h stdlib.h string.h sys/socket.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_CHECK_HEADER_STDBOOL
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_FORK
AC_CHECK_FUNCS([dup2 gettimeofday memmove memset select socket sqrt strdup strerror strstr strtol])

AC_DEFINE([USE_MALLOC_WRAPPERS], [], [])
AC_DEFINE([HAVE_PTHREAD], [], [])


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
