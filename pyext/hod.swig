%module hod
%{
    #define SWIG_FILE_WITH_INIT
    #include "hod.h"
%}

%include "carrays.i"
%array_functions(double, dArray);
%array_class(double, doubleArray);
%include "cosmology.h"
%include "hod.h"
