AUTOMAKE_OPTIONS = subdir-objects
BUILT_SOURCES = $(srcdir)/cosmology_wrap.cpp  $(srcdir)/hod_wrap.cpp

COS_SOURCES = cosmology.i 
HOD_SOURCES = hod.i 

pkgpython_PYTHON = cosmology.py hod.py
pkgpyexec_LTLIBRARIES = _cosmology.la _hod.la

_cosmology_la_SOURCES = $(srcdir)/cosmology_wrap.cpp $(COS_SOURCES)
_cosmology_la_CPPFLAGS = $(AX_SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_cosmology_la_LDFLAGS = -module
_cosmology_la_LIBADD = ../src/libcosmology.la

_hod_la_SOURCES = $(srcdir)/hod_wrap.cpp $(HOD_SOURCES)
_hod_la_CPPFLAGS = $(AX_SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_hod_la_LDFLAGS = -module
_hod_la_LIBADD = ../src/libhod.la

$(srcdir)/cosmology_wrap.cpp : $(COS_SOURCES)
	$(SWIG) -c++ $(AX_SWIG_PYTHON_OPT) -I$(top_srcdir)/src -o $@ $<

$(srcdir)/hod_wrap.cpp : $(HOD_SOURCES)
	$(SWIG) -c++ $(AX_SWIG_PYTHON_OPT) -I$(top_srcdir)/src -o $@ $<

clean:
	rm -f $(BUILT_SOURCES) *.py *.lo *.la
