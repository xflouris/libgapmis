#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([libgapmis], [0.0.9], [Tomas.Flouri@h-its.org])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/gapmis_cpu.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
# AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h float.h stdint.h stdlib.h string.h sys/time.h unistd.h])

AC_ARG_ENABLE([opt], [ --enable-opt=sse|gpu    disable sse or gpu optimization], [opt=${enableval}], [opt=cpu])
have_pthreads=no
AC_SEARCH_LIBS([pthread_create], [pthread], [have_pthreads=yes])

if test "x${have_pthreads}" = xyes; then
  AC_CHECK_HEADERS([pthread.h], [], [have_pthreads=no])
fi

have_opencl=yes
AC_CHECK_HEADERS([CL/opencl.h], [], [have_opencl=no])

if test "x${have_opencl}" = xno; then
  if test "x${opt}" = xgpu; then
    AC_MSG_WARN([
    ---------------------------------------------------------
    Unable to find the OpenCL header files which are required
    for building the GPU version of libgapmis.
    Building only the non-optimized version.
    ---------------------------------------------------------])
    opt=cpu
  fi
fi

if test "x${have_pthreads}" = xno; then
  if test "x${opt}" = xsse; then
    AC_MSG_WARN([
    --------------------------------------------------------
    Unable to find pthreads on this system which is required
    for building the SSE version of libgapmis.
    Building only the non-optimized version.
    --------------------------------------------------------])
    opt=cpu
  fi
fi
AM_CONDITIONAL(OPT_SSE, test "$opt" = "sse")
AM_CONDITIONAL(OPT_GPU, test "$opt" = "gpu")
AM_CONDITIONAL(OPT_CPU, test "$opt" = "cpu")


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_CHECK_FUNCS([clock_gettime getpagesize memmove strchr strstr])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 examples/Makefile
                 man/Makefile])
AC_OUTPUT
