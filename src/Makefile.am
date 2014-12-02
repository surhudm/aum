bin_PROGRAMS=aumtest
lib_LTLIBRARIES = libcosmology.la libhod.la
libcosmology_la_SOURCES = cosmology.cpp powerspectrum.cpp haloes.cpp gauleg.cpp 
libhod_la_SOURCES = hod.cpp cosmology.cpp powerspectrum.cpp haloes.cpp gauleg.cpp 
pkginclude_HEADERS = cosmology.h gauleg.h hod.h
aumtest_SOURCES=test.cpp
aumtest_LDADD=libhod.la
