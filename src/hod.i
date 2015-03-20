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
>>> p = h.cosmo()
>>> q = h.hodpars()
>>> p.Om0 = 0.307115
>>> p.w0 = -1
>>> p.wa = 0
>>> p.Omk = 0.0
>>> p.hval = 0.6777
>>> p.Omb = 0.048206
>>> p.th = 2.726
>>> p.s8 = 0.8228
>>> p.nspec = 0.96
>>> p.ximax = log10(8.0)
>>> p.cfac = 1.0
>>> q.Mmin = 13.0
>>> q.siglogM = 0.5
>>> q.Msat = 14.0
>>> q.alpsat = 1.0
>>> q.Mcut = 13.5
>>> q.csbycdm = 1.0
>>> q.fac = 1.0
>>> a = h.hod(p, q)
>>> help(a)

"

%feature("docstring") hod::ncen
"Average number of central galaxies in halo of mass M

:Parameters:

-   log M200: logarithm of the mass of the halo

:Returns:

-   Average number of central galaxies in the halo


:Examples:

>>> import hod as h
>>> a = h.hod()
>>> a.ncen(12.0)

"

%feature("docstring") hod::nsat
"Average number of satellite galaxies in halo of mass M

:Parameters:

-   log M200: logarithm of the mass of the halo

:Returns:

-   Average number of satellite galaxies in the halo


:Examples:

>>> import hod as h
>>> a = h.hod()
>>> a.nsat(12.0)

"

%feature("docstring") hod::nsatz
"Average number density of satellite galaxies at redshift z

:Parameters:

-   z: Redshift

:Returns:

-   Average number density of satellite galaxies at redshift z


:Examples:

>>> import hod as h
>>> a = h.hod()
>>> a.nsatz(12.0)

"

%feature("docstring") hod::ncenz
"Average number density of central galaxies at redshift z

:Parameters:

-   z: Redshift

:Returns:

-   Average number density of central galaxies at redshift z


:Examples:

>>> import hod as h
>>> a = h.hod()
>>> a.ncenz(12.0)

"

%feature("docstring") hod::avmass_tot
"Average halo mass of all galaxies

:Parameters:

-   z: Redshift

:Returns:

-   Average halo mass of all galaxies at redshift z


:Examples:

>>> import hod as h
>>> a = h.hod()
>>> a.avmass_tot(12.0)

"

%feature("docstring") hod::avmass_cen
"Average halo mass of central galaxies

:Parameters:

-   z: Redshift

:Returns:

-   Average halo mass of central galaxies at redshift z


:Examples:

>>> import hod as h
>>> a = h.hod()
>>> a.avmass_cen(12.0)

"

%feature("docstring") hod::galaxy_bias
"Average halo bias of all galaxies at redshift z

:Parameters:

-   z: Redshift

:Returns:

-   Average halo bias of all galaxies at redshift z


:Examples:

>>> import hod as h
>>> a = h.hod()
>>> a.galaxy_bias(12.0)

"

%feature("docstring") hod::gets8
"Output value of sigma8

:Parameters:

-   None : No inputs

:Returns:

-   sigma8 : sigma8

:Examples:

>>> a.gets8()
"

%feature("docstring") hod::getOmb
"Output value of Omegab

:Parameters:

-   None : No input parameters

:Returns:

-   Omegab : Baryon density parameter

:Examples:

>>> a.getOmb()
"

%feature("docstring") hod::geth
"Output value of h value

:Parameters:

-   None : No input parameters

:Returns:

-   h : Hubble parameter

:Examples:

>>> a.geth()
"


%feature("docstring") hod::getOmk
"Output value of curvature parameter

:Parameters:

-   None : No input parameters

:Returns:

-   Omk : Curvature parameter

:Examples:

>>> a.getOmk()
"

%feature("docstring") hod::set_cfactor
"Set multiplicative constant to multiply all dark matter concentrations

:Parameters:

-   cfac : multiplicative constant

:Returns:

-   None : No return value

:Examples:

>>> a.set_cfactor(1.0)
"

%feature("docstring") cosmology::D2gg_num
"Power per logarithmic k interval in the galaxy galaxy power spectrum Delta^2(k,z)

:Parameters:

-   k: Wavenumber (in h Mpc^{-1})
-   z: Redshift

:Returns:

-   Delta_gg^2(k,z)

:Examples:

>>> a.D2gg_num(0.1,0.0)

"

%feature("docstring") cosmology::D2gd_num
"Power per logarithmic k interval in the galaxy matter power spectrum Delta^2(k,z)

:Parameters:

-   k: Wavenumber (in h Mpc^{-1})
-   z: Redshift

:Returns:

-   Delta_gd^2(k,z)

:Examples:

>>> a.D2gd_num(0.1,0.0)

"

%feature("docstring") cosmology::xigg_num
"Galaxy-galaxy correlation function at distance radius and redshift z

:Parameters:

-   r: Wavenumber (in h Mpc^{-1})
-   z: Redshift

:Returns:

-   xi_gg(r,z)

:Examples:

>>> a.xigg_num(0.1,0.0)

"

%feature("docstring") cosmology::xigd_num
"Galaxy-matter correlation function at distance radius and redshift z

:Parameters:

-   r: Wavenumber (in h Mpc^{-1})
-   z: Redshift

:Returns:

-   xi_gd(r,z)

:Examples:

>>> a.xigd_num(0.1,0.0)

"


%include "cosmology.h"
%include "hod.h"
