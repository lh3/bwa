#               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([bwa], [0.7.15], [kmzhou4@yahoo.com])
AM_INIT_AUTOMAKE([-Wall -Werror])
AM_PROG_AR
AC_CONFIG_SRCDIR([src/is.c])
AC_CONFIG_HEADERS([config.h])
LT_PREREQ([2.4])
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_LN_S
AC_PROG_MAKE_SET

# harded coded in the original makefile
AC_DEFINE([USE_MALLOC_WRAPPERS], [1], [Hard coded in bwa makefile])

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lbwa':
#AC_CHECK_LIB([bwa], [main])
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([rt], [shm_open])
AC_CHECK_LIB([z], [deflate])

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
AC_FUNC_ALLOCA
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_FORK
AC_CHECK_FUNCS([dup2 ftruncate gettimeofday memmove memset munmap select socket sqrt strdup strerror strstr strtol])

AC_CONFIG_FILES([Makefile
					  src/Makefile])
AC_OUTPUT
