///Cosmology routines
#include "cosmology.h"

///Default constructor WMAP3 flat cosmology
cosmology::cosmology()
{
    Omega0 = 0.238;
    Omegab = 0.02233/pow(0.734,2.);
    Omegak=0.0;
    Omegal = 1.-Omega0-Omegak;
    h      = 0.734 ;
    theta  = 1.0093;///2.7;
    sigma8 = 0.744;
    ns     = 0.951;
    w0     =-1.0;
    wa     = 0.0;
    xiNLzetamax=log10(8.0);
    cfactor=1.0;
    initialize();
    if(verbose) std::cout<<"# Cosmo constructor 1 called\n";
}

/** Initialize the cosmology dependent routines
 *  and output to the log
 * */
void cosmology::initialize()
{
    /// Initialize options
    opt_mf=1;
    opt_b=1;
    opt_c=1;
    opt_ps_L=1;
    opt_ps_NL=1;

    /// Initialize the time today
    t0=Time(0.0);

    /// Init z_glob
    z_glob=0.0;
    gf_glob=1.0;

    /// Initialise the boolean variables
    bool_init_Chi=false;
    bool_initPS=false;
    bool_init_PSL0=false;
    bool_init_PSNL=false;
    bool_init_xiL0=false;
    bool_init_xiNL=false;
    bool_initSmith=false;
    bool_init_varM_TH_spline=false;
    bool_init_GF=false;
    bool_inituk=false;
    init_Tink=false;
    bool_initQk=false;
    bool_init_xiNL_bar=false;
    bool_init_xiNL_barbar=false;
    bool_init_xiL_bar=false;
    bool_init_xiL_barbar=false;

    /// Fix some constants
    d2norm=1.0;
    rho_crit_0=3.0e4/(8.*M_PI*gee);
    facmtor=(1./(4./3.*M_PI*Omega0*rho_crit_0));

    /// Output cosmology
    if(verbose){
        std::cout<<"# "<<"============================================================="<<std::endl;
        std::cout<<"# "<<"Cosmology initialised to: "<<std::endl;
        std::cout<<"# Om0="<<Omega0<<"\n" \
            <<"# Oml="<<Omegal<<"\n" \
            <<"# Omk="<<Omegak<<"\n" \
            <<"#  w0="<<w0<<"\n" \
            <<"#  wa="<<wa<<"\n" \
            <<"# Omb="<<Omegab<<"\n" \
            <<"# hub="<<h     <<"\n" \
            <<"# tht="<<theta <<"\n" \
            <<"# s8 ="<<sigma8<<"\n" \
            <<"# ns ="<<ns<<"\n" \
            <<"# xinlzetamax ="<<xiNLzetamax<<"\n" \
            <<"# cfactor ="<<cfactor            \
            <<std::endl;
        std::cout<<"# t0 ="<<t0<<", ,rho_crit_0="<<rho_crit_0<<std::endl;
        std::cout<<"# "<<"============================================================="<<std::endl;
    }

    /// Fix the normalisation of the power spectrum
    fixd2norm();

    /// Initialise the growth factor on a spline. I think this will be done by default by d2norm
    init_growthfactor();

    /// Get gauleg masses
    if(!mock){
        gauleg(9.0,16.0,x9_16,w9_16,N9_16);
    }else{
        /// Get gauleg masses for mock L250
        //    gauleg(12.0,14.5,x9_16,w9_16,N9_16);
        /// get gauleg masses for mock lgpc
        //    gauleg(13.0,14.5,x9_16,w9_16,n9_16);
        /// get gauleg masses for mock Takahiro
        gauleg(log10(1.7724E12),16.0,x9_16,w9_16,N9_16);
        //gauleg(log10(1.7724E9),16.0,x9_16,w9_16,N9_16);

    }

    /// Get gauleg pi integrals
    gauleg(0.0,2.0*M_PI,x0_2p,w0_2p,N0_2p);

    init_powerspectra_L();

    // f_g(M) swindle
    fgm_m0=-99.0;
    fgm_slp=-99.0;

}

///Constructor for cosmology by passing a cosmo struct
cosmology::cosmology(cosmo p)
{
    Omega0 = p.Om0;
    Omegab = p.Omb;
    w0     = p.w0 ;
    wa     = p.wa ;
    Omegak = p.Omk;
    Omegal = 1.-p.Om0-p.Omk;
    h      = p.hval;
    theta  = p.th/2.7;
    sigma8 = p.s8;
    ns     = p.nspec;
    xiNLzetamax=p.ximax;
    cfactor=p.cfac;
    initialize();
    if(verbose) std::cout<<"# Cosmo constructor 2 called\n";
}

///Constructor for cosmology by passing all cosmological parameters explicitly
cosmology::cosmology(double om0, double omk, double xw0, double xwa,
        double omb, double xh ,double th, double s8, double nspec,double ximax,
        double cfac)
{
    Omega0 = om0;
    Omegab = omb;
    Omegak = omk;
    Omegal = 1.0-om0-omk;
    w0     = xw0;
    wa     = xwa;
    h      = xh;
    theta  = th/2.7;
    sigma8 = s8;
    ns     = nspec;
    xiNLzetamax=ximax;
    cfactor=cfac;
    initialize();
    if(verbose) std::cout<<"# Cosmo constructor 3 called\n";
}

///Destructor for cosmology object
cosmology::~cosmology()
{
    cosmo_free();
}


/// Integral for time from Big Bang: Valid only late in the matter dominated era.
double dTime(double x, void * params)
{
    c_params c1 = *(c_params *) params;
    //    std::cout<<1<<std::endl;
    cosmology *c2;
    c2 = c1.cptr;
    double Om0=(*c2).Omega0;
    double w0 =(*c2).w0;
    double wa =(*c2).wa;
    double scale=1./(1.+x);
    double Oml=(*c2).Omegal*exp(3.*( -(1+w0+wa)*log(scale) - wa*(1-scale)  ));
    double temp=Om0*pow(x,-3.0)+Oml+(1.-Om0-(*c2).Omegal)*pow(x,-2.0);
    if(temp<=0.0){
        std::cout<<"ERROR: Calculating time for this cosmology not valid.\n";
        std::cout<<"ERROR: wa is "<<wa<<"\n";
        std::cout<<"ERROR: factors are "<<x<<" "<<Om0*pow(x,-3.0)<<" "<<Oml<<" "<<(1.-Om0-Oml)*pow(x,-2.0)<<"\n";
        std::cout<<"ERROR: temp is "<<temp<<" for scale factor "<<x<<"\n";
        exit(0);
    }
    return 1./(x*sqrt(temp));
}

/// Lookback time in units of 1/H0
double cosmology::Lookback(double z)
{
        /*
    /// Int_{(1+z)^{-1}}^{1} 1/E(a) da/a
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    double result, error;

    gsl_function F;
    F.function = &(dTime);

    //Initialize parameter block for dlookback
    c_params p;
    p.cptr = this;

    F.params = &p;

    gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
    gsl_integration_qags (&F, 1./(1.+z), 1., 0, 1e-6, 1000, w, &result, &error); 
    gsl_set_error_handler(oldhand);

    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);

    gsl_integration_workspace_free (w);

    return result;
    */
    return t0-Time(z);
}

/// Time in units of 1/H0
double cosmology::Time(double z)
{
    /// Int_0^{(1+z)^{-1}} 1/E(a) da/a
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    double result, error;

    gsl_function F;
    F.function = &(dTime);

    //Initialize parameter block for dlookback
    c_params p;
    p.cptr = this;

    F.params = &p;

    gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
    gsl_integration_qags (&F, 0., 1./(1.+z), 0, 1e-6, 1000, w, &result, &error); 
    gsl_set_error_handler(oldhand);

    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);

    gsl_integration_workspace_free (w);

    return result;
}

/// Expansion factor (z)
double cosmology::Eofz(double z)
{
    return sqrt(Omega0*pow(1.+z,3.0)+Omegal*exp( 3.*( (1+w0+wa)*log(1.+z) - wa*z/(1.+z)  ) )+(1-Omega0-Omegal)*pow(1.+z,2.0));
}

/// Omega(z)
double cosmology::Omega(double z)
{
    return Omega0*pow(1+z,3.)/pow(Eofz(z),2.);
}

/// Delta_crit(z) a'la Bryan and Norman '98
/// This function is only defined for flat cosmologies and for cosmologies with Omegal=0.0
/// If this function is required then please find a substitute
double cosmology::Delta_crit(double z)
{
    double result;

    double x=Omega(z)-1.;

    if (Omega0<1.&&Omegal==0.0)
    {
        result=18.0 * M_PI*M_PI + 60.0*x - 32.0 * x*x;
    } else if (fabs(Omega0+Omegal-1.0)<1.e-4)
    {
        result=18.0 * M_PI*M_PI + 82.0*x - 39.0 * x*x;
    } else
    {
        std::cout<<"# Aborting as Delta_Crit not defined for this cosmology."<<std::endl;
        exit(0);
    }

    return result;
}

/// Integrand for comoving distance
double dChi(double x, void * params)
{
    c_params c1 = *(c_params *) params;
    //    std::cout<<1<<std::endl;
    cosmology *c2;
    c2=c1.cptr;
    double Om0=(*c2).Omega0;
    double w0 =(*c2).w0;
    double wa =(*c2).wa;
    double scale=1./(1.+x);
    double Oml=(*c2).Omegal*exp(3.*( -(1+w0+wa)*log(scale) - wa*(1-scale)  ));
    return 1./sqrt(Om0*pow(1.+x,3.0)+Oml+(1.-Om0-(*c2).Omegal)*pow(1.+x,2.0));
}

/// Chi(z) : Comoving distance in units of h^{-1} Mpc
double cosmology::Chiofz(double z)
{
    /// Int_0^{z} 1/E(z) dz
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    double result, error;

    gsl_function F;
    F.function = &(dChi);

    // Initialize parameter block for dChi
    c_params p;
    p.cptr = this;

    F.params = &p;

    gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
    gsl_integration_qags (&F, 0., z, 0, 1e-6, 1000, w, &result, &error); 
    gsl_set_error_handler(oldhand);

    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);

    gsl_integration_workspace_free (w);

    return result*c/100.0;
}

/// Chiofz initialization
void cosmology::init_Chi(){
    if(!bool_init_Chi){
        if(verbose){
            std::cout<<"# "<<"============================================================="<<std::endl;
            std::cout<<"# "<<"The spline for numerical calculation of comoving distance was not initialized. Initializing..."<<std::endl;
        }
        //Define the limits for which to calculate the variance
        double zmin=0.0;
        double zmax=6.0;
        double zdiff=(zmax-zmin)/(10*Nsigma-1.);
        double xx[10*Nsigma],yy[10*Nsigma];
        for(int i=0;i<10*Nsigma;i++){
            xx[i]=zmin+zdiff*i;
            yy[i]=Chiofz(xx[i]);
        }
        Chi_acc = gsl_interp_accel_alloc ();
        Chi_spline = gsl_spline_alloc (gsl_interp_cspline, 10*Nsigma);
        gsl_spline_init (Chi_spline,xx,yy,10*Nsigma);

        bool_init_Chi=true;
        if(verbose){
            std::cout<<"# "<<"The spline for numerical calculation of comoving distance is now initialized"<<std::endl;
            std::cout<<"# "<<"============================================================="<<std::endl;
        }
    }
}

double cosmology::Chiofz_num(double z){
    if(z>6.0 || z<0.0){
        return Chiofz(z);
    }
    if(!bool_init_Chi){
        init_Chi();
    }
    double result=gsl_spline_eval (Chi_spline,z, Chi_acc);
    return result;
}

// Comoving distance in units of h^{-1} Mpc
double cosmology::Dcofz(double z){
    return (1+z)*Daofz(z);
}

/// Dl(z) : Luminosity distance in units of h^{-1} Mpc
double cosmology::Dlofz(double z)
{
    return (1+z)*(1+z)*Daofz(z);
}

/// Da(z) : Angular diameter distance in units of h^{-1} Mpc
double cosmology::Daofz(double z)
{
    double scale=fabs(pow((100.0/c),2.0)*Omegak);
    double Chi=Chiofz_num(z);
    double da=0.0;
    if(Omegak==0.0)
    {
        da=Chi/(1.+z);
    } else if ( Omegak < 0.0 )
    {
        da=1./sqrt(scale)*sin(sqrt(scale)*Chi)/(1.+z);
    } else if ( Omegak > 0.0 )
    {
        da=1./sqrt(scale)*sinh(sqrt(scale)*Chi)/(1.+z);
    }
    return da;
}

// Da(z1,z2) : Angular diameter distance in units of h^{-1} Mpc
double cosmology::Daofzlh(double zl,double zh)
{
    double Omegak=1-Omega0-Omegal;
    double scale=fabs(pow((100.0/c),2.0)*Omegak);
    double Chi1=Chiofz(zl);
    double Chi2=Chiofz(zh);
    double da;
    if(Omegak==0.0)
    {
        da=(Chi2-Chi1)/(1.+zh);
    }  else if ( Omegak < 0.0 )
    {
        da=1./sqrt(scale)*sin(sqrt(scale)*(Chi2-Chi1))/(1.+zh);
    } else if ( Omegak > 0.0 )
    {
        da=1./sqrt(scale)*sinh(sqrt(scale)*(Chi2-Chi1))/(1.+zh);
    }
    return da;
}


/// New friend functions for growth factor
double E_sq(gf_par &cp, double &z){
    return (cp.Omega0*pow(1.+z,3.) + cp.OmegaL*exp( 3*(1.+cp.w0+cp.wa)*log(1.+z) - 3*cp.wa + 3.*cp.wa/(1.+z) ) + (1-cp.Omega0-cp.OmegaL)*pow(1.+z,2.) );
}

double dE_sqdz(gf_par &cp, double &z){
    double temp=exp( 3*(1.+cp.w0+cp.wa)*log(1.+z) - 3*cp.wa + 3.*cp.wa/(1.+z) );
    return (cp.Omega0*3.*pow(1.+z,2.) +  cp.OmegaL*temp*( 3.*(1.+cp.w0+cp.wa)*1./(1.+z) -  3.*cp.wa/pow(1.+z,2.)  )  + (1-cp.Omega0-cp.OmegaL)*2.*(1.+z)  );
}

void getall(gf_par &cp, double &z, double &Esqt, double &dEsqdzt,double &d2lnEsqdz2t){
    Esqt=E_sq(cp,z);
    double temp=exp( 3*(1.+cp.w0+cp.wa)*log(1.+z) - 3*cp.wa + 3.*cp.wa/(1.+z) );
    dEsqdzt=(cp.Omega0*3.*pow(1.+z,2.) +  cp.OmegaL*temp*( 3.*(1.+cp.w0+cp.wa)*1./(1.+z) -  3.*cp.wa/pow(1.+z,2.)  )  + (1-cp.Omega0-cp.OmegaL)*2.*(1.+z)  );

    double temp2=( 3.*(1.+cp.w0+cp.wa)*1./(1.+z) -  3.*cp.wa/pow(1.+z,2.)  );
    double d2Esqdz2t= cp.Omega0*6.*(1.+z) +  cp.OmegaL*temp*( temp2*temp2 - 3.*(1.+cp.w0+cp.wa)*1./pow(1.+z,2.) + 6*cp.wa/pow(1.+z,3.) )  + (1-cp.Omega0-cp.OmegaL)*2. ;

    d2lnEsqdz2t=1./Esqt*( d2Esqdz2t - pow(dEsqdzt,2.)/Esqt) ;
}

double d2lnE_sqdz2(gf_par &cp, double &z){
    double Esqt=E_sq(cp,z);
    double dEsqdzt=dE_sqdz(cp,z);
    double temp=exp( 3*(1.+cp.w0+cp.wa)*log(1.+z) - 3*cp.wa + 3.*cp.wa/(1.+z) );
    double temp2=( 3.*(1.+cp.w0+cp.wa)*1./(1.+z) -  3.*cp.wa/pow(1.+z,2.)  );
    double d2Esqdz2t= cp.Omega0*6.*(1.+z) +  cp.OmegaL*temp*( temp2*temp2 - 3.*(1.+cp.w0+cp.wa)*1./pow(1.+z,2.) + 6*cp.wa/pow(1.+z,3.) )  + (1-cp.Omega0-cp.OmegaL)*2. ;

    return 1./Esqt*( d2Esqdz2t - pow(dEsqdzt,2.)/Esqt) ;
}

int gf_func (double z, const double y[], double f[], void *params)
{
    gf_par cp=*(gf_par *)params;
    double Esqt=E_sq(cp,z);
    f[0] = y[1];
    f[1] = 3.*cp.Omega0*(1.+z)*y[0]/(2.*Esqt)  - y[1]*dE_sqdz(cp,z)/(2.*Esqt) + y[1]/(1.+z);
    return GSL_SUCCESS;
}

int gf_jac (double z, const double y[], double *dfdy, double dfdz[], void *params)
{
    gf_par cp=*(gf_par *)params;
    double Esqt,dEsqdzt,d2lnEsqdz2t;

    getall(cp,z,Esqt,dEsqdzt,d2lnEsqdz2t);

    gsl_matrix_view dfdy_mat
        = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0,  3.*cp.Omega0*(1.+z)/(2.*Esqt) );
    gsl_matrix_set (m, 1, 1, 1./(1.+z) - dEsqdzt/(2.*Esqt) );
    dfdz[0] = 0.0;
    dfdz[1] = -3./2.*cp.Omega0*y[0]*1./Esqt*( (1.+z)*dEsqdzt/Esqt - 1.)
        - y[1]/pow(1.+z,2.)
        - y[1]/2.*d2lnEsqdz2t;
    return GSL_SUCCESS;
}
/// End of friend functions for growth factor

/// Growth factor (z) numerical implementation. Initializing spline. 
void cosmology::init_growthfactor()
{
    if(!bool_init_GF){
        if(verbose){
            std::cout<<"# "<<"============================================================="<<std::endl;
            std::cout<<"# "<<"The spline for numerical calculation of growth factor was not initialized. Initializing..."<<std::endl;
        }
        //Define the limits for which to calculate the variance
        double zmin=0.0;
        double zmax=10.0;

        /// Integrator
        const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;

        gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
        gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-6, 0.0);
        gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (2);

        gf_par cp;
        cp.Omega0=Omega0;
        cp.OmegaL=Omegal;
        cp.w0=w0;
        cp.wa=wa;
        gsl_odeiv_system sys = {gf_func, gf_jac, 2, &cp};

        double z=1088.0;
        double h=-1.e-6;

        // Normalising the growth factor to redshift 1088 for the time being
        // g(z)=1089./(1.+z) is a good approx. at high z such as z>1088
        // y[1]=g'(z) is then defined as -1089./pow((1.+z),2.)
        // Initialise!
        double y[2] = { 1.0, -1./1089. };

        double xx[Ngf],yy[Ngf];
        for (int i=Ngf-1;i>=0;i--)
        {
            double zz=zmin+i/(Ngf-1.)*(zmax-zmin);
            xx[i]=zz;

            while (z > xx[i])
            {
                int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &z, xx[i], &h, y);

                if (status != GSL_SUCCESS)
                {
                    std::cout<<"Error encountered while calculating growthfactor"<<std::endl;
                    exit(0);
                    break;
                }

                //printf ("%.5e %.5e %.5e\n", z, y[0], y[1]);
            }
            yy[i]=y[0];

        }
        //exit(0);

        /// Free up the ode solver
        gsl_odeiv_evolve_free (e);
        gsl_odeiv_control_free (c);
        gsl_odeiv_step_free (s);

        /// Normalise back to growth factor at redshift 0
        for (int i=Ngf-1;i>=0;i--)
        {
            yy[i]=yy[i]/yy[0];
        }
        /// All the growth factors are properly normalized with respect to growth factor at redshift zero now! 
        /// So continue old style from here on!

        GF_acc = gsl_interp_accel_alloc ();
        GF_spline = gsl_spline_alloc (gsl_interp_cspline, Ngf);
        gsl_spline_init (GF_spline,xx,yy,Ngf);

        bool_init_GF=true;

        /// Output for checking
        /*
           for (int i=0;i<Ngf;i++)
           {
           std::cout<<xx[i]<<" "<<yy[i]<<std::endl;
           }
           exit(0);
         */

        if(verbose){
            std::cout<<"# "<<"The spline for numerical calculation of growth factor is now initialized"<<std::endl;
            std::cout<<"# "<<"============================================================="<<std::endl;
        }
    }

}

/// Growthfactor
double cosmology::growthfactor(double z)
{
    {
        std::cout<<"Error consider increasing zmax from 10 to"<<z<<std::endl;
        exit(0);
    }
}

/// Calculate the growth factor at redshift z numerically by interpolating
double cosmology::growthfactor_num( double z)
{
    if(!bool_init_GF)
    {
        init_growthfactor();
    }
    double result;
    if(z<=10.0&&z>=0.0)
    {
        result=gsl_spline_eval (GF_spline,z, GF_acc);
    } else
    {
        //std::cout<<"# "<<"GF_spline: Interpolation not possible z not within range "<<z<<std::endl;
        result=growthfactor(z);
        //exit(0);
    }
    return result;

}

double cosmology::dlnDdln1pz(double z){
    double result=gsl_spline_eval_deriv(GF_spline,z, GF_acc);
    result=result*(1.+z)/gsl_spline_eval (GF_spline,z, GF_acc);
    return result;
}

/// Access to private variables
double cosmology::gets8()
{
    return sigma8;
}
double cosmology::getOmb()
{
    return Omegab;
}
double cosmology::geth()
{
    return h;
}
double cosmology::getns()
{
    return ns;
}
double cosmology::getzetamax()
{
    return zetamax;
}
double cosmology::getzeta_rmax()
{
    return zeta_rmax;
}
double cosmology::getxinlzetamax()
{
    return xiNLzetamax;
}
double cosmology::get_cfac()
{
    return cfactor;
}
double cosmology::set_cfac(double cfac)
{
    cfactor=cfac;
    return cfactor;
}

void cosmology::cosmo_free(){
    if(verbose){
        std::cout<<"# Cosmo free destructor\n";
    }
    if(bool_init_Chi){
        gsl_interp_accel_free(Chi_acc);
        gsl_spline_free(Chi_spline);
        bool_init_Chi=false;
    }
    if(bool_init_PSL0){
        gsl_interp_accel_free(PSL0_acc);
        gsl_spline_free(PSL0_spline);
        bool_init_PSL0=false;
    }
    if(bool_init_PSNL){
        gsl_interp_accel_free(PSNL_acc);
        gsl_spline_free(PSNL_spline);
        bool_init_PSNL=false;
    }
    if(bool_init_xiNL){
        gsl_interp_accel_free(xiNL_acc);
        gsl_spline_free(xiNL_spline);
        bool_init_xiNL=false;
    }
    if(bool_init_xiL0){
        gsl_interp_accel_free(xiL0_acc);
        gsl_spline_free(xiL0_spline);
        bool_init_xiL0=false;
    }
    if(bool_init_varM_TH_spline){
        gsl_interp_accel_free(varM_TH_num_acc);
        gsl_spline_free(varM_TH_num_spline);
        bool_init_varM_TH_spline=false;
    }
    if(bool_init_GF){
        gsl_interp_accel_free(GF_acc);
        gsl_spline_free(GF_spline);
        bool_init_GF=false;
    }
    if(bool_inituk){
        gsl_interp_accel_free(uk_c_acc);
        gsl_interp_accel_free(uk_krs_acc);
        bool_inituk=false;
    }
    if(bool_initQk)
    {
        delete []Qk1;
        delete []Qk2;
        bool_initQk=false;
    }
}

double cosmology::getLmin(double z, double L1){
    // m - M = 5.0 log10 (dLum(z)/10pc)
    // M - 5 log h = 17.77 - 5.0 log10 (dLum(z) h/Mpc) - 25.0
    double Mmax=17.77-5.0*log10(Dlofz(z))-25.0;
    double Lmin=(Mmax-4.76)/(-2.5);
    if(Lmin>L1){
        return Lmin;
    }else{
        return L1;
    }
}

double findzmax(double z, void * params){
    z_params c1 = *(z_params *) params;
    cosmology *c2;
    double *mag;
    c2=c1.cptr;
    mag=c1.mag;
    return 17.77-( *mag + 5.0*log10((*c2).Dlofz(z)) + 25.0  );
}

/// Get maximum redshift to which a galaxy with luminosity L can be seen
double cosmology::getzmax(double xL){
    double mag=-2.5*xL+4.76;
    /// Find the redshift where the apparent magnitude equals 17.77
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    // m is in log_{10}(m)
    double z_lo = 1E-5, z_hi = 10.0;

    gsl_function F;
    z_params zpar;
    zpar.cptr=this;
    zpar.mag=&mag;

    F.function = &findzmax;
    F.params = &zpar;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, z_lo, z_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        z_lo = gsl_root_fsolver_x_lower (s);
        z_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (z_lo, z_hi,0, 1e-6);

        if (status == GSL_SUCCESS)
        {
            //std::cout<<"# "<<"getM:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return res;

}

double cosmology::Omegaw(double z){
    return Omegal*exp( 3.*( (1+w0+wa)*log(1.+z) - wa*z/(1.+z)  ) )/pow(Eofz(z),2.);
}
void cosmology::set_optmf(int opt){
    opt_mf=opt;
}

void cosmology::renew(cosmo p){
    // First free up the original constructors
    cosmo_free();

    // Now change the cosmology and reinitialize
    Omega0 = p.Om0;
    Omegab = p.Omb;
    w0     = p.w0 ;
    wa     = p.wa ;
    Omegak = p.Omk;
    Omegal = 1.-p.Om0-p.Omk;
    h      = p.hval;
    theta  = p.th/2.7;
    sigma8 = p.s8;
    ns     = p.nspec;
    xiNLzetamax=p.ximax;
    cfactor=p.cfac;
    initialize();
    if(verbose) std::cout<<"# Cosmo parameter renewal\n";
}

double cosmology::rsound(){
    return sd;
}

double cosmology::get_deltapi(double z1, double z2){
    return c*fabs(z1-z2)/100./Eofz((z1+z2)/2.);
}

double cosmology::get_sinsqang(double x1, double x2, double x3, double y1, double y2, double y3){
    double cosang=x1*y1+x2*y2+x3*y3;
    return 1.-cosang*cosang;
}

double cosmology::get_logrp(double x1, double x2, double x3, double y1, double y2, double y3,double Chisq){
    double cosang=x1*y1+x2*y2+x3*y3;
    if(cosang>1.0) return -100.0;
    return 0.5*log10(Chisq*(1.-cosang*cosang));
}
