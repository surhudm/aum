%module cosmology
%include cpointer.i
%pointer_class(double,dp)
%{
    #define SWIG_FILE_WITH_INIT
    #include "cosmology.h"
%}
%include "cosmology.h"
%apply double& INOUT { double& a };
