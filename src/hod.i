%module hod
%include cpointer.i
%pointer_class(double,dp)
%apply double& INOUT { double& a };
%feature("autodoc", 1);
%include "carrays.i"
%include "cosmology.i"
%array_functions(double, dArray);
%array_class(double, doubleArray);
%{
    #define SWIG_FILE_WITH_INIT
    #include "hod.h"
%}

%feature("docstring") hod::hod
"Initializes hod and cosmology object

:Parameters:

-  cosmo structure:
    Om0 : Matter density parameter
    Omk : Curvature parameter
    w0 : Dark energy equation of state parameter
    wa : Dark energy equation of state parameter
    Omb : Baryon density parameter
    h : Hubble parameter
    th : CMB temperature
    s8 : sigma8
    nspec : power spectrum index
    ximax : Parameter psi defined in van den Bosch 2013, only relevant for halo model calculation
    cfac : Constant multiplicative factor for the c-M relation
-  hodp structure:
    Mmin : Minimum halo mass in the central HOD
    siglogM : Scatter in halo masses in the central HOD
    Msat : Satellite halo occupation mass scale
    alpsat : Slope of the satellite halo occupation
    Mcut : Cut off mass scale of the satellite halo occupation
    fac : Unused parameter
    csbycdm : Multiplicative factor for satellite concentrations


:Returns:

-   HOD object

Without any inputs, initializes to flat WMAP3 LCDM cosmology, cfac=1.0, ximax=log10(8.0).

:Examples:

    >>> import hod as h
    >>> a = h.cosmo()
    >>> h = h.hodpars()
    >>> a.Om0=0.3
    >>> help(a)

"

%include "cosmology.h"
%include "hod.h"
