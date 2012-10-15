ACLOCAL_AMFLAGS = -I m4 --install

lib_LTLIBRARIES = libgapmis.la
libgapmis_la_SOURCES = gapmis_cpu.c EDNAFULL.h EBLOSUM62.h
libgapmis_la_CFLAGS = -g -Wall -msse3 -O3 -fomit-frame-pointer -funroll-loops

if OPT_SSE
lib_LTLIBRARIES += libgapmis_sse.la
libgapmis_sse_la_SOURCES = gapmis_sse.cc sse_utils/thread.h sse_utils/cycle.h sse_utils/sse_unit.h sse_utils/aligned_buffer.h EDNAFULL.h EBLOSUM62.h
libgapmis_sse_la_CPPFLAGS = -O3 -Wall -march=native -D_USE_SSE -lstdc++ -lpthread
endif
if OPT_GPU
lib_LTLIBRARIES += libgapmis_gpu.la
libgapmis_gpu_la_SOURCES = gapmis_gpu.c EDNAFULL.h EBLOSUM62.h
libgapmis_gpu_la_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -Wall -msse3 -D_USE_GPU -lOpenCL
endif

pkginclude_HEADERS = gapmis.h errors.h

