%module cosmology
%include cpointer.i
%pointer_class(double,dp)
%apply double& INOUT { double& a };
%feature("autodoc", 1);
%{
    #define SWIG_FILE_WITH_INIT
    #include "cosmology.h"
%}

%feature("docstring") cosmology::cosmology
"Initializes cosmology object

:Parameters:

-   Omega0 : Matter density parameter
-   OmegaK : Curvature parameter
-   w0 : Dark energy equation of state parameter
-   wa : Dark energy equation of state parameter
-   Omegab : Baryon density parameter
-   h : Hubble parameter
-   ThetaCMB : CMB temperature
-   sigma8 : sigma8
-   ns : power spectrum index
-   psi : Parameter psi defined in van den Bosch 2013, only relevant for halo model calculation
-   cfac : Constant multiplicative factor for the c-M relation

:Returns:

-   Cosmology object

Without any inputs, initializes to flat WMAP3 LCDM cosmology, cfac=1.0, ximax=log10(8.0).

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> help(a)

"

%feature("docstring") cosmology::cosmo_free
"Frees all memory associated with cosmology object

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.cosmo_free()

"


%feature("docstring") cosmology::Dcofz
" Comoving distance as a function of redshift

:Parameters:

-   z : Redshift

:Returns:

-   Comoving distance

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.Dcofz(0.5)

"

%feature("docstring") cosmology::get_sinsqang
" Calculate the square of the sin of the angle between two galaxies

:Parameters:

-   x1 : Cartesian x for the unit vector pointing at Galaxy 1
-   y1 : Cartesian y for the unit vector pointing at Galaxy 1
-   z1 : Cartesian z for the unit vector pointing at Galaxy 1
-   x2 : Cartesian x for the unit vector pointing at Galaxy 2
-   y2 : Cartesian y for the unit vector pointing at Galaxy 2
-   z2 : Cartesian z for the unit vector pointing at Galaxy 2

:Returns:

-   Square of the sin of the angle between two galaxies

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.get_sinsqang(0.0, 1.0, 0.0, 1.0, 0.0, 0.0)

"
%feature("docstring") cosmology::get_logrp
" Calculate the projected separation between two galaxies

:Parameters:

-   x1 : Cartesian x for the unit vector pointing at Galaxy 1
-   y1 : Cartesian y for the unit vector pointing at Galaxy 1
-   z1 : Cartesian z for the unit vector pointing at Galaxy 1
-   x2 : Cartesian x for the unit vector pointing at Galaxy 2
-   y2 : Cartesian y for the unit vector pointing at Galaxy 2
-   z2 : Cartesian z for the unit vector pointing at Galaxy 2

:Returns:

-   Log of projected separation between two galaxies

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.get_logrp(0.0, 1.0, 0.0, 1.0, 0.0, 0.0)

"


%feature("docstring") cosmology::get_deltapi
" Calculate the line of sight separation between two galaxies

:Parameters:

-   z1 : Redshift
-   z2 : Redshift

:Returns:

-   Line of sight separation between two galaxies

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.get_deltapi(0.0, 0.2)

"

%feature("docstring") cosmology::Dlofz
" Luminosity distance as a function of redshift

:Parameters:

-   z : Redshift

:Returns:

-   Luminosity distance

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.Dlofz(0.5)

"


%feature("docstring") cosmology::Daofz
" Angular diameter distance as a function of redshift

:Parameters:

-   z : Redshift

:Returns:

-   Angular diameter distance

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.Daofz(0.5)

"

%feature("docstring") cosmology::Daofzlh
" Angular diameter distance as a function of redshift of lens and source

:Parameters:

-   zl : Redshift of lens
-   zh : Redshift of source

:Returns:

-   Angular diameter distance between two redshifts

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.Daofzlh(0.5,1.0)

"

%feature("docstring") cosmology::growthfactor_num
" Growth factor as a function of redshift (normalized to unity at redshift zero)

:Parameters:

-   z : Redshift

:Returns:

-   Growth factor at redshift z

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.growthfactor_num(0.5)

"

%feature("docstring") cosmology::dlnDdln1pz
" Negative of the logarithmic derivative of growth factor with scale factor

:Parameters:

-   z : Redshift

:Returns:

-   Negative of the logarithmic derivative of growth factor with scale factor

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.dlnDdln1pz(0.5)

"


%feature("docstring") cosmology::Omega
" Matter density parameter at redshift z

:Parameters:

-   z : Redshift

:Returns:

-   Matter density parameter at redshift z

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.Omega(0.5)

"


%feature("docstring") cosmology::Omegaw
" Dark energy density parameter at redshift z

:Parameters:

-   z : Redshift

:Returns:

-   Dark energy density parameter at redshift z

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.Omegaw(0.5)

"


%feature("docstring") cosmology::Delta_crit
" Virial density contrast at redshift z a'la Bryan and Norman '98

:Parameters:

-   z : Redshift

:Returns:

-   Virial density contrast (with respect to critical density at redshift z)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.Delta_crit(0.5)

"

%feature("docstring") cosmology::set_optmf
" Set mass function option

:Parameters:

-   option = 1: Tinker et al. 2010 mass function (well tested and consistent
    with the bias prescription
-   option = 2: Sheth Tormen mass function
-   option = 3: Bhattacharya et al. 2010 mass function

:Returns:

-   Set mass function choice (default is equal to 1 if this function is not called)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.set_optmf(1)

"

%feature("docstring") cosmology::Delta2_L_num
" Power per logarithmic k interval in the linear matter power spectrum Delta^2(k,z)

:Parameters:

-   k: Wavenumber (in h Mpc^{-1})
-   z: Redshift

:Returns:

-   Delta^2(k,z)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.Delta2_L_num(0.1,0.0)

"

%feature("docstring") cosmology::Delta2_NL_num
" Power per logarithmic k interval in the nonlinear matter power spectrum Delta^2_NL(k,z)

:Parameters:

-   k: Wavenumber (in h Mpc^{-1})
-   z: Redshift

:Returns:

-   Nonlinear Delta^2(k,z)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.Delta2_NL_num(0.1,0.0)

"

%feature("docstring") cosmology::xi_L_num
" Linear matter correlation function

:Parameters:

-   r: Scale (in hinv Mpc)
-   z: Redshift

:Returns:

-   Linear matter correlation function (r,z)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.xi_L_num(0.1,0.0)

"

%feature("docstring") cosmology::xi_NL_num
" Non-Linear matter correlation function

:Parameters:

-   r: Scale (in hinv Mpc)
-   z: Redshift

:Returns:

-   Non-linear matter correlation function (r,z)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.xi_NL_num(0.1,0.0)

"

%feature("docstring") cosmology::nofm
" Mass function as a function of mass and redshift

:Parameters:

-   M: Mass (in hinv Msun)
-   z: Redshift

:Returns:

-   dN(>M)/dM of halos, where N(>M) is the cumulative number density of halos with mass larger than M, commonly referred to as mass function

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.nofm(1e12,0.0)

"

%feature("docstring") cosmology::bias
" Halo bias function as a function of mass and redshift

:Parameters:

-   M: Mass (in hinv Msun)
-   z: Redshift

:Returns:

-   Large scale halo bias

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.bias(1e12,0.0)

"

%feature("docstring") cosmology::varM_TH_num
" Variance of fluctuations on a given mass scale [sigma^2(M,z)]

:Parameters:

-   M: Mass (in hinv Msun)
-   z: Redshift

:Returns:

-   Variance of fluctuations when density field is smoothed on the lagrangian radius corresponding to a given mass scale 

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.varM_TH_num(1e12,0.0)

"

%feature("docstring") cosmology::varM_TH_num_deriv
" dln sigma^2/dln M

:Parameters:

-   M: Mass (in hinv Msun)
-   z: Redshift

:Returns:

- dln sigma^2/dln M  (M,z)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.varM_TH_num_deriv(1e12,0.0)

"


%feature("docstring") cosmology::Nplus
" Number density of halos with mass above a given mass at a given redshift

:Parameters:

-   M: Mass (in hinv Msun)
-   z: Redshift

:Returns:

- N(>M,z)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.Nplus(1e12,0.0)

"

%feature("docstring") cosmology::getM
" Find mass such that halos with mass larger than it have a given number density at a given redshift

:Parameters:

-   Nplus: Target number density (in h^3 Mpc^{-3})
-   z: Redshift

:Returns:

-   M: Mass (hinv Msun)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.getM(1e-6,0.0)

"

%feature("docstring") cosmology::modelNFWhalo
" Output the virial mass, physical virial radius, virial
concentration of a halo, its physical radius with density
contrast 200m and the corresponding concentration c200m
given a mass M200m at redshift z

:Parameters:

-   M200m: Mass (hinv Msun) defined 200 times overdense with respect to the background
-   z: Redshift

:Returns:

-   Mvir : The virial mass (hinv Msun)
-   Rvir : The physical virial radius (hinv Mpc)
-   cvir : The virial concentration
-   R200m : The physical boundary of halo 200 times overdense with respect to background density (hinv Mpc)
-   c200m : The concentration of halo 200 times overdense with respect to background density

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.modelNFWhalo(1e12,0.0)

"

%feature("docstring") cosmology::modelNFWhalo_com
" Output the virial mass, comoving virial radius, virial
concentration of a halo, its comoving radius with density
contrast 200m and the corresponding concentration c200m
given a mass M200m at redshift z

:Parameters:

-   M200m: Mass (hinv Msun) defined 200 times overdense with respect to the background
-   z: Redshift

:Returns:

-   Mvir : The virial mass (hinv Msun)
-   Rvir : The comoving virial radius (hinv Mpc)
-   cvir : The virial concentration
-   R200m : The comoving boundary of halo 200 times overdense with respect to background density (hinv Mpc)
-   c200m : The concentration of halo 200 times overdense with respect to background density

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.modelNFWhalo_com(1e12,0.0)

"

%feature("docstring") cosmology::conc
" 
Concentration of halos

:Parameters:

-   Mvir: Virial mass (hinv Msun) 
-   z: Redshift

:Returns:

-   cvir : The virial concentration

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
    >>> a.conc(1e12,0.0)

"

%feature("docstring") cosmology::Eofz
"Returns the cosmological expansion function E(z)

:Parameters:

-   z : Redshift

:Returns:

-   Eofz: Expansion function

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.Eofz(0.5)
    1.28111279753
"

%feature("docstring") cosmology::setnew_z
" Reset the global redshift at which many of the splines in the cosmology code are initialized. This is rarely used function.

:Parameters:

-   z : Redshift

:Returns:

-   None: None

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.setnew_z(0.5)
"


%feature("docstring") cosmology::gets8
" Output value of sigma8

:Parameters:

-   None : No inputs

:Returns:

-   sigma8 : sigma8

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.gets8()
"

%feature("docstring") cosmology::getOmb
" Output value of Omegab

:Parameters:

-   None : No input parameters

:Returns:

-   Omegab : Baryon density parameter

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.getOmb()
"

%feature("docstring") cosmology::geth
" Output value of h value

:Parameters:

-   None : No input parameters

:Returns:

-   h : Hubble parameter

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.geth()
"


%feature("docstring") cosmology::getns
" Output value of spectral index

:Parameters:

-   None : No input parameters

:Returns:

-   ns : Spectral index of initial density fluctuations

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.getns()
"

%feature("docstring") cosmology::getxinlzetamax
" Returns the value of psi from van den Bosch 2013

:Parameters:

-   None : No input parameters

:Returns:

-   psi : Psi defined in van den Bosch 2013

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.getxinlzetamax()
"

%feature("docstring") cosmology::get_cfac
" Returns factor multiplying all concentrations calculated by the code

:Parameters:

-   None : No input parameters

:Returns:

-   cfac : Factor multiplying all concentrations output by the code

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.get_cfac()
"

%feature("docstring") cosmology::set_cfac
" Sets factor multiplying all concentrations calculated by the code

:Parameters:

-   cfac : Factor multiplying all concentrations output by the code

:Returns:

-   None : None

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.set_cfac(1.0)
"

%feature("docstring") cosmology::getzmax
" Get the maximum redshift to which a galaxy can be observed by SDSS spectroscopic survey

:Parameters:

-   xL: Luminosity in h^{-2} Lsun

:Returns:

-   z : Redshift

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.getzmax(1e12)
"


%feature("docstring") cosmology::getLmin
" Get the minimum luminosity of galaxies that can be observed by SDSS spectroscopic survey at a given redshift

:Parameters:

-   z : Redshift

:Returns:

-   xL: Luminosity in h^{-2} Lsun

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.getLmin(0.1)
"


%feature("docstring") cosmology::Time
" Get the time in units of 1/(H_0 km/s/Mpc/yr)

:Parameters:

-   z : Redshift

:Returns:

-   Time : in units of 1/(H_0 km/s/Mpc/yr)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.Time(0.1)
"


%feature("docstring") cosmology::Lookback
" Get the lookback time in units of 1/(H_0 km/s/Mpc/yr)

:Parameters:

-   z : Redshift

:Returns:

-   Lookback Time : in units of 1/(H_0 km/s/Mpc/yr)

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.Lookback(0.1)
"


%feature("docstring") cosmology::wpnl
" Get the projected non-linear matter correlation function

:Parameters:

-   r : Separation of galaxies
-   z : Redshift
-   pimax : Line-of-sight integration limit

:Returns:

-   wpnl : Projected non-linear matter correlation function

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.wpnl(0.1,0.1,100.0)
"

%feature("docstring") cosmology::wpl
" Get the projected linear matter correlation function

:Parameters:

-   r : Separation of galaxies
-   z : Redshift
-   pimax : Line-of-sight integration limit

:Returns:

-   wpl : Projected linear matter correlation function

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.wpl(0.1,0.1,100.0)
"


%feature("docstring") cosmology::wpnl_kaiser
" Get the projected non-linear matter correlation function accounting for Kaiser effects

:Parameters:

-   r : Separation of galaxies
-   z : Redshift
-   pimax : Line-of-sight integration limit
-   fkai : Kaiser factor f/b

:Returns:

-   wpnl_kaiser : Projected non-linear matter correlation function

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.wpnl_kaiser(0.1,0.1,100.0,1.0)
"

%feature("docstring") cosmology::wpl_kaiser
" Get the projected linear matter correlation function accounting for Kaiser effect

:Parameters:

-   r : Separation of galaxies
-   z : Redshift
-   pimax : Line-of-sight integration limit
-   fkai : Kaiser factor f/b

:Returns:

-   wpl : Projected linear matter correlation function

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.wpl_kaiser(0.1,0.1,100.0,1.0)
"


%feature("docstring") cosmology::xi_NL_kaiser
" Get the non-linear matter correlation function accounting for Kaiser effects

:Parameters:

-   r : Separation of galaxies
-   z : Redshift
-   mu : Cosine of angle between separation vector and line of sight
-   fkai : Kaiser factor f/b

:Returns:

-   xi_nl_kaiser : Non-linear matter correlation function

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.xi_nl_kaiser(0.1,0.1,0.5,1.0)
"

%feature("docstring") cosmology::xi_L_kaiser
" Get the linear matter correlation function accounting for Kaiser effect

:Parameters:

-   r : Separation of galaxies
-   z : Redshift
-   mu : Cosine of angle between separation vector and line of sight
-   fkai : Kaiser factor f/b

:Returns:

-   xi_l : Linear matter correlation function

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.xi_L_kaiser(0.1,0.1,0.5,1.0)
"

%feature("docstring") cosmology::rsound
" Get the comoving sound horizon at the drag epoch

:Parameters:

-   None : None

:Returns:

-   rsound : Comoving sound horizon at drag epoch a'la Eisenstein and Hu 98

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.rsound()
"

%include "cosmology.h"
