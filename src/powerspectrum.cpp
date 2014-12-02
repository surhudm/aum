/// Power spectrum routines part of cosmology class
#include "cosmology.h"

/// Linear Powerspectrum wrappers
double cosmology::Delta2_L(double k,double z)
{
    double result;

    if(opt_ps_L==1)
    {
        result=Delta2_EH(k, z);
    } else
    {
        std::cout<<"Linear power spectrum option not supported yet."<<std::endl;
        exit(0);
    }

    return result;

}

double cosmology::Pk_L(double k,double z)
{
    double result;

    if(opt_ps_L==1)
    {
        result=Pk_EH(k, z);
    } else
    {
        std::cout<<"Linear power spectrum option not supported yet."<<std::endl;
        exit(0);
    }

    return result;

}

/// Nonlinear Powerspectrum wrappers
double cosmology::Pk_NL(double k,double z)
{
    double result;

    if(opt_ps_NL==1)
    {
        result=PkNL_S(k, z);
    } else
    {
        std::cout<<"Non-linear power spectrum option not supported yet."<<std::endl;
        exit(0);
    }

    return result;

}

double cosmology::Delta2_NL(double k,double z)
{
    double result;

    if(opt_ps_NL==1)
    {
        result=Delta2NL_S(k, z);
    } else
    {
        std::cout<<"Non-linear power spectrum option not supported yet."<<std::endl;
        exit(0);
    }

    return result;

}

/// Initialize power spectrum variables Eisenstein and Hu 1998
double cosmology::initPS_EH()
{

    double O0h2=Omega0*h*h;
    double Obh2=Omegab*h*h;

    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<" O0h2="<<O0h2<<" "<<" Obh2="<<Obh2<<std::endl;
    }

    /// Equations refer to Eisenstein and Hu '98.
    // Eq. 2
    zeq = 2.50*10000*O0h2*pow(theta,-4.);
    //std::cout<<"# "<<" zeq="<<zeq<<" "<<O0h2<<" "<<theta<<std::endl;
    //std::cout<<"# "<<" zeq="<<zeq<<std::endl;

    // Eq. 3
    //keq = sqrt(2*O0h2*pow(100./c,2.)*zeq);
    keq = 0.0746*Omega0*h*h*pow(theta,-2.);
    //std::cout<<"# "<<" keq="<<keq<<std::endl;//" "<<7.46e-2*Omega0*h*h*pow(theta,-2.)<<std::endl;

    //Eq. 4, define this to be 1+zdrag, needs a fix in all releases.
    double b1=0.313*pow(O0h2,-0.419)*(1.+0.607*pow(O0h2,0.674));
    double b2=0.238*pow(O0h2,0.223);
    zd = 1+1291*pow(O0h2,0.251)*(1.+b1*pow(Obh2,b2))/(1.+0.659*pow(O0h2,0.828));
    //std::cout<<"# "<<" zd="<<zd<<std::endl;

    //Eq. 5
    Rd = 31.5*Obh2*pow(theta,-4.)*(1000./zd);
    Req = 31.5*Obh2*pow(theta,-4.)*(1000./zeq);
    //std::cout<<"# "<<" Rd="<<Rd<<std::endl;
    //std::cout<<"# "<<" Req="<<Req<<std::endl;

    //Eq.6
    sd = 2./(3.*keq)*sqrt(6./Req)*log( ( sqrt(1.+Rd) + sqrt(Rd+Req) )/(1.+sqrt(Req))  );
    if(verbose){
    std::cout<<"# "<<" Sound horizon at drag epoch in comoving Mpc="<<sd<<std::endl;
    std::cout<<"# "<<" Angular sound horizon at drag epoch="<<sd/(1+zd)*h/Daofz(zd)<<std::endl;
    }
    //std::cout<<"# "<<" sd="<<sd<<std::endl;

    //Eq. 7
    //double factor=pow(10.4*O0h2,-0.95);
    ksilk = 1.6*pow(Obh2,0.52)*pow(O0h2,0.73)*(1.+pow(10.4*O0h2,-0.95));
    //std::cout<<"# "<<" ksilk="<<ksilk<<std::endl;

    //Eq. 11
    double a1 = pow(46.9*O0h2,0.670)*(1.+pow(32.1*O0h2,-0.532));
    double a2 = pow(12.0*O0h2,0.424)*(1.+pow(45.0*O0h2,-0.582));
    alphac = pow(a1,-Omegab/Omega0)*pow(a2,-pow(Omegab/Omega0,3.));
    //std::cout<<"# "<<" alphac="<<alphac<<std::endl;

    //Eq. 12
    double Omegac=Omega0-Omegab;
    double bb1 = 0.944/(1.+pow(458*O0h2,-0.708));
    double bb2 = pow(0.395*O0h2,-0.0266);
    betacinv=1.+bb1*(pow(Omegac/Omega0,bb2)-1.);
    //std::cout<<"# "<<" betac="<<1./betacinv<<std::endl;

    //Eq. 14
    double y = zeq/zd;
    double Gy= y*( -6.*sqrt(1.+y) + (2.+3.*y)*log( (sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)  )  );
    alphab=2.07*keq*sd*pow(1.+Rd,-0.75)*Gy;
    //std::cout<<"# "<<" alphab="<<alphab<<" "<<Gy<<" "<<log( (sqrt(1.+y)+1.)/(sqrt(1.+y)-1.))<<std::endl;
    //std::cout<<"# "<<-6.*sqrt(1.+y)<<" "<<(2.+3.*y)*log( (sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)  )<<std::endl;

    //Eq. 23
    betanode=8.41*pow(O0h2,0.435);
    //std::cout<<"# "<<" betanode="<<betanode<<std::endl;

    //Eq. 24
    betab = 0.5 + Omegab/Omega0 + (3. - 2.*Omegab/Omega0)*sqrt(1.+pow(17.2*O0h2,2.));
    //std::cout<<"# "<<" betab="<<betab<<std::endl;

    bool_initPS=true;

    if(verbose){
    std::cout<<"# "<<"Eisenstein and Hu Power spectrum variables initialised."<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }

    return 1;

}

/// T0master for Eisenstein and Hu power spectrum
double cosmology::T0master_EH(double q, double alpc, double betc)
{
    // Eq. 20
    // The value q = k/(keq*13.41) is passed instead of k to avoid repititive
    // calculation;
    double CT0 = 14.2/alpc+386./(1+69.9*pow(q,1.08));

    //std::cout<<" Check master "<<CT0<<" "<<alpc<<" "<<betc<<" "<<q<<std::endl;
    
    // Eq. 21 This is also called TF_pressureless
    return log( e + 1.8*betc*q )/( log(e+1.8*betc*q)+ CT0*q*q);

}

/*
// k should be supplied in Mpc^{-1}
double cosmology::TCold_EH(double k)
{
    //Eq. 18
    double f= 1./(1.+pow(k*sd/5.4,4.));
    double q=k/(keq*13.41);
    
    //Eq. 17
    return f*T0master_EH(q,1.,1./betacinv)+(1-f)*T0master_EH(q,alphac,1./betacinv);

}

// k should be supplied in Mpc^{-1}
double cosmology::Tb_EH(double k)
{
    double q=k/(keq*13.41);
    // Eq. 21 first part
    double befj0= T0master_EH(q,1.,1.)/(1.+pow(k*sd/5.2,2.)) + alphab/(1.+pow(betab/(k*sd),3.))*exp(-pow(k/ksilk,1.4) );

    // Eq. 22
    double stilde=sd/pow((1.+pow(betanode/(k*sd),3.)),1./3.);

    //std::cout<<" Check baryonic part "<<befj0<<" "<<stilde<<" "<<q<<std::endl;
    //Eq. 21 remaining part
    double out;
    //if(fabs(k*stilde)>1.0e-4)
    //{ 
        out=befj0*gsl_sf_bessel_j0(k*stilde);
        //out=befj0*sin(k*stilde)/(k*stilde);
    //}else
    //{
    //    out=befj0;
    //}
    return out;
}*/

/// Transfer function a'la Eisenstein and Hu, k should be supplied in Mpc^{-1}
double cosmology::TF_EH(double kpassed)
{

    //double k = h*kpassed;
    double k = kpassed;
    if (!bool_initPS) 
    {
	if(verbose){
        std::cout<<"# "<<" Power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = initPS_EH(); 

	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }

    double q=k/(keq*13.41);
    // Eq. 21 first part
    double befj0= T0master_EH(q,1.,1.)/(1.+pow(k*sd/5.2,2.)) + alphab/(1.+pow(betab/(k*sd),3.))*exp(-pow(k/ksilk,1.4) );
    // Eq. 22
    double stilde=sd/pow((1.+pow(betanode/(k*sd),3.)),1./3.);
    // Eq. 21 remaining part
    double Tb=befj0*gsl_sf_bessel_j0(k*stilde);

    // Eq. 18
    double f= 1./(1.+pow(k*sd/5.4,4.));
    // Eq. 17
    double TC=f*T0master_EH(q,1.,1./betacinv)+(1-f)*T0master_EH(q,alphac,1./betacinv);
    double TF=Omegab/Omega0*Tb + (Omega0-Omegab)/Omega0*TC;
    return TF;

    //return Omegab/Omega0*Tb_EH(k) + (Omega0-Omegab)/Omega0*TCold_EH(k);
}

/// Delta^2(k) a'la Eisenstein and Hu, k should be in Mpc^{-1}
/// kpassed should be supplied as h Mpc^{-1}.
/// Result is dimensionless as it should be
double cosmology::Delta2_EH(double kpassed, double z)
{
    double k = h*kpassed;

    if (!bool_initPS) 
    {
	if(verbose){
        std::cout<<"# "<<" Power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = initPS_EH(); 
	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }
    //Normalize the scale to 100./c; If d2norm is 1, then the absolute value as
    //such has no meaning. d2norm should be fixed using sigma8. Check routine,
    //fix_d2norm()
    //if(d2norm==1.){std::cout<<"Did you normalise the power spectrum?";}
    double res=d2norm*pow(c*kpassed/(100.*h),3.+ns)*pow(TF_EH(k),2.);
    if(z!=z_glob){
        res=res*pow(growthfactor_num(z),2.);
    }else{
        res=res*pow(gf_glob,2.);
    }

    return res;
}

/// P(k) a'la Eisenstein and Hu, k should be supplied in Mpc^{-1}
/// kpassed should be supplied as h Mpc^{-1}.
/// Result is in the units of (h^{-1} Mpc)^3
double cosmology::Pk_EH(double kpassed, double z)
{
    double k = kpassed;

    double res=Delta2_EH(k,z);

    return res*(2*M_PI*M_PI)/(pow(k,3.));
}

/// Integrand for var_TH: 1/k \Delta^2(k) (3.*j1(kR)/(kR))^2  k in hMpc^{-1}
double dvar_TH(double x, void * params)
{
    cvar_params c1 = *(cvar_params *) params;
//    std::cout<<1<<std::endl;
    cosmology *c2;
    double *R;
    double *z;
    bool *psi;
    c2=c1.cptr;
    R=c1.R;
    z=c1.z;
    psi=c1.psinit;
    double arg=(*R)*x;
    double d2;
    if(*psi){
	d2=(*c2).Delta2_L_num(x,*z);
    }else{
	d2=(*c2).Delta2_L(x,*z);
    }
    double bes= pow(3.*gsl_sf_bessel_j1(arg)/(arg) ,2.);
    //std::cout<<1./x*d2*bes<<"Is this quick"<<std::endl;
    return 1./x*d2*bes;
}

/*
// Integrand for var_TH: 1/k \Delta^2(k) (3.*j1(kR)/(kR))^2  k in hMpc^{-1}
double dvar_TH2(double x, void * params)
{
    cvar_params c1 = *(cvar_params *) params;
//    std::cout<<1<<std::endl;
    cosmology *c2;
    double *R;
    double *z;
    c2=c1.cptr;
    R=c1.R;
    z=c1.z;
    
    double arg1=(*R)*x;
    double arg2=1./(*R)*x;
    
    double d2_1=(*c2).Delta2_EH(x,*z);
    double d2_2=(*c2).Delta2_EH(1./x,*z);

    double bes_1= pow(3.*gsl_sf_bessel_j1(arg1)/(arg1) ,2.);
    double bes_2= pow(3.*gsl_sf_bessel_j1(arg2)/(arg2) ,2.);

    //std::cout<<1./x*(d2_1*bes_1+d2_2*bes_2)<<std::endl;
    return 1./x*(d2_1*bes_1+d2_2*bes_2);
}*/

/// Calculate the variance of the power spectrum on a scale R in h^{-1} Mpc
double cosmology::var_TH(double R, double z)
{
    ///Int_{0}^{\infty} dk/k Delta^2(k) (3j1(kR)/kR)^2
    double result1, result2, error;

    gsl_function F;
    F.function = &(dvar_TH);

    //Initialize parameter block
    cvar_params p;
    p.cptr = this;
    p.R=&R;
    p.z=&z;
    p.psinit=&bool_init_PSL0;

    F.params = &p;

    // Do it in two steps as suggested by Frank. First upto 1/Rf and then from
    // 1/Rf to infty.
    int status=0;
    gsl_error_handler_t *oldhand=gsl_set_error_handler_off();

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    status=gsl_integration_qags (&F, 0., 1./R, 0, 1e-6, 1000, w, &result1, &error); 
    if(status!=0 and verbose){
        std::cout<<"Proceeding despite convergence problems"<<std::endl;
    }
    gsl_integration_workspace_free (w);
    
    gsl_integration_workspace * v = gsl_integration_workspace_alloc (1000);
    status=gsl_integration_qagiu (&F, 1./R, 0, 1e-6, 1000, v, &result2, &error); 
    if(status!=0 and verbose){
        std::cout<<"Proceeding despite convergence problems"<<std::endl;
    }
    gsl_integration_workspace_free (v);
    gsl_set_error_handler(oldhand);
    
    //gsl_integration_qagi (&F, 0, 1e-4, 1000, w, &result, &error); 

    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);


    return result1+result2;
    
}

/// Calculate the variance of the power spectrum on a mass scale M in h^{-1} Msun
double cosmology::varM_TH(double M, double z)
{
    ///Calculate R of top hat for M
    double r=pow(M*facmtor,1./3.);
    //std::cout<<M<<" "<<r<<std::endl;
    return var_TH(r,z);
    
}

/// Initialize the spline for numerically calculating the variance
/// This is valid for redshift z=0.0. Should be appropriately multiplied by the
/// growth factor to extend to redshift z.
void cosmology::init_varM_TH_spline()
{
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for numerical calculation of variance was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance
    double mmin=5.0;
    double mmax=18.0;

    double xx[Nsigma],yy[Nsigma];
    for (int i=0;i<Nsigma;i++)
    {
        double m=mmin+i/(Nsigma-1.)*(mmax-mmin);
        xx[i]=m;
        m=pow(10.,m);
        yy[i]=log(varM_TH(m,0.0))/log(10.);
        if(verbose){
        std::cout<<"Variance: mass, sigma "<<xx[i]<<" "<<yy[i]<<std::endl;
        }
    }
    varM_TH_num_acc = gsl_interp_accel_alloc ();
    varM_TH_num_spline = gsl_spline_alloc (gsl_interp_cspline, Nsigma);
    gsl_spline_init (varM_TH_num_spline,xx,yy,Nsigma);

    bool_init_varM_TH_spline=true;
    if(verbose){
    std::cout<<"# "<<"The spline for numerical calculation of variance is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
}

/// Calculate the variance of the power spectrum on a mass scale M in h^{-1} Msun
/// numerically by interpolating
double cosmology::varM_TH_num(double M, double z)
{
    //
    if(!bool_init_varM_TH_spline)
    {
        //std::cout<<"# I called it: varM_TH_spline"<<std::endl;
        init_varM_TH_spline();
    }
    double logm=log(M)/log(10.);
    double result;
    if(logm>5.0&&logm<18.0)
    {
        result=gsl_spline_eval (varM_TH_num_spline,logm, varM_TH_num_acc);
    } else
    {
        std::cout<<"# "<<"varM_TH_num_spline: Interpolation not possible M not within range "<<M<<std::endl;
        //exit(0);
        result=log10(varM_TH(M,0.0));
    }
    result = pow(10.,result);
    if(z!=z_glob){
        result=result*pow(growthfactor_num(z),2.);
    }else{
        result=result*pow(gf_glob,2.);
    }
    return result;

}

/// Calculate dlog\sigma/dlogM for a mass M in h^{-1} Msun
/// numerically using the spline
double cosmology::varM_TH_num_deriv(double M, double z)
{
    //
    if(!bool_init_varM_TH_spline)
    {
        //std::cout<<"# I called it: varM_TH_num_deriv"<<std::endl;
        init_varM_TH_spline();
    }
    double logm=log(M)/log(10.);
    double result;
    if(logm>5&&logm<18.0)
    {
        result=gsl_spline_eval_deriv (varM_TH_num_spline,logm, varM_TH_num_acc);
    } else
    {
        std::cout<<"# "<<"varM_TH_num_spline: Interpolation for derivative not possible, M not within range "<<M<<std::endl;
        exit(0);
        result=varM_TH(M,z);
    }
    result = pow(10.,result);
    if(z!=z_glob){
        result=result*pow(growthfactor_num(z),2.);
    }else{
        result=result*pow(gf_glob,2.);
    }
    return result;

}

/// Integrand for var_G: 1/k \Delta^2(k) exp(-(kR)^2)  k in hMpc^{-1}
double dvar_G(double x, void * params)
{
    cvar_params c1 = *(cvar_params *) params;
//    std::cout<<1<<std::endl;
    cosmology *c2;
    double *R;
    double *z;
    c2=c1.cptr;
    R=c1.R;
    z=c1.z;
    double d2=(*c2).Delta2_L(x,*z);
    double res=1./x*d2*exp(-pow(x*(*R),2));
    return res;
}

/// Calculate the variance of the power spectrum on a scale R in h^{-1} Mpc
/// using a Gaussian filter. 
double cosmology::var_G(double R, double z)
{
    //Smith et al. 2003 Eq. 54
    //Int_0^{\infty} dk/k Delta^2(k) exp(-(kR)^2)
    double result1, result2, error;

    gsl_function F;
    F.function = &(dvar_G);

    //Initialize parameter block
    cvar_params p;
    p.cptr = this;
    p.R=&R;
    p.z=&z;

    F.params = &p;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F, 0., 1./R, 0, 1e-6, 1000, w, &result1, &error); 
    gsl_integration_workspace_free (w);

    gsl_integration_workspace * v = gsl_integration_workspace_alloc (1000);
    gsl_integration_qagiu (&F, 1./R, 0, 1e-6, 1000, v, &result2, &error); 
    gsl_integration_workspace_free (v);

    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);


    return result1+result2;
    
}

/// Fix the normalization of the power spectrum using the value of \f$\sigma_8\$f supplied by the user.
void cosmology::fixd2norm()
{
    double sig82=var_TH(8.,0.0);
    d2norm=pow(sigma8,2)/sig82;
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The power spectrum has been normalized. Checking:";
    }
    if(fabs(var_TH(8.,0.0)-pow(sigma8,2))>1.0e-4)
    {
        std::cout<<"# "<<"Problem in fixing the power spectrum normalization. Quitting"<<sigma8<<" "<<var_TH(8.,0.0)<<" "<<bool_initPS<<std::endl;
        exit(0);
    }
    if(verbose){
    std::cout<<"# "<<"Check successful."<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
}

/// Root of this function yields k_{\sigma} for Smith et al. nonlinear power spectrum
double findksig(double xk, void *params)
{
    ksig_params c1 = *(ksig_params *) params;
    cosmology *c2;
    double *z;
    c2=c1.cptr;
    z=c1.z;
    double r=1./pow(10.0,xk);
    double var=(*c2).var_G(r,*z);
    // Smith et al. 2003 Eq. 53
    //std::cout<<xk<<" "<<r<<" "<<var-1.;
    return var-1.;

}

/// Find k_{\sigma} for Smith et al. nonlinear power spectrum
double cosmology::getksigma(double z)
{
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    // k is in log_{10}(k)
    double xk_lo = -1.7, xk_hi = 1.0;

    gsl_function F;
    ksig_params p;
    p.cptr = this;
    p.z = &z;

    F.function = &findksig;
    F.params = &p;
   
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, xk_lo, xk_hi);

    do
    {
        //std::cout<<iter<<std::endl;
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        xk_lo = gsl_root_fsolver_x_lower (s);
        xk_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (xk_lo, xk_hi,0, 1e-6);

        if (status == GSL_SUCCESS)
        {
	if(verbose){
            std::cout<<"# "<<"init_Smith:ksigma:Brent converged after "<< iter<<" iterations"<<std::endl;
	}
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return pow(10.,res);

}

/// Integrand for neff for Smith et al. nonlinear power spectrum
double dneffint(double k, void* params)
{
    /// Int_0^{\infty} dk 1./k Delta^2(k) y^2 exp(-y^2)
    ksig_params c1 = *(ksig_params *) params;
//    std::cout<<1<<std::endl;
    cosmology *c2;
    double *z;
    c2=c1.cptr;
    z=c1.z;
    
    //double k=pow(10.0,logk);
    double y=k/(*c2).ksigma;
    
    return 1./k*(*c2).Delta2_L(k,*z)*pow(y,2.)*exp(-pow(y,2.));

}


/// neff for Smith et al. nonlinear power spectrum
double cosmology::getneff(double z)
{
    /// Smith et al. 2003 Eq. 59
    /// neff= 2*(Int_0^{\infty} dk/k Delta^2(k) y^2 exp(-y^2))- 3. (with y=k/ksig)
    // First get the integral
    double result1, result2, error;

    gsl_function F;
    F.function = &(dneffint);

    //Initialize parameter block
    ksig_params p;
    p.cptr = this;
    p.z=&z;

    F.params = &p;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F, 0., exp(1.), 0, 1e-6, 1000, w, &result1, &error); 
    gsl_integration_workspace_free (w);
    gsl_integration_workspace * u = gsl_integration_workspace_alloc (1000);
    gsl_integration_qagiu (&F, exp(1.), 0, 1e-6, 1000, u, &result2, &error); 
    gsl_integration_workspace_free (u);
    /*
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F, -1., 1., 0, 1e-6, 1000, w, &result2, &error); 
    gsl_integration_workspace_free (w);
    gsl_integration_workspace * u = gsl_integration_workspace_alloc (1000);
    gsl_integration_qagil (&F, -1., 0, 1e-6, 1000, u, &result1, &error); 
    gsl_integration_workspace_free (u);
    gsl_integration_workspace * v = gsl_integration_workspace_alloc (1000);
    gsl_integration_qagiu (&F, 1., 0, 1e-6, 1000, v, &result3, &error); 
    gsl_integration_workspace_free (v);*/

    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);


    //return 2.*(result1+result2+result3)-3.;
    return 2.*(result1+result2)-3.;
    

}

/// Integrand for Curvature C for Smith et al. nonlinear power spectrum
//double dCint(double logk, void* params)
double dCint(double k, void* params)
{
    // Int_0^{\infty} dk 1./k Delta^2(k) (y^2-y^4) exp(-y^2)
    ksig_params c1 = *(ksig_params *) params;
    cosmology *c2;
    double *z;
    c2=c1.cptr;
    z=c1.z;
    //double k=pow(10.0,logk);
    
    double y=k/(*c2).ksigma;
    //std::cout<<" ksigma"<<(*c2).ksigma<<" "<<*z<<" "<<(*c2).growthfactor_num(*z)<<std::endl;
    
    return (1./k)*(*c2).Delta2_L(k,*z)*(pow(y,2.)-pow(y,4.))*exp(-pow(y,2.));

}


/// Curvature for Smith et al. nonlinear power spectrum
double cosmology::getC(double z)
{
    //Smith et al. 2003 Eq. 60
    //C= (3.+neff)^2+4*(Int_0^{\infty} dk/k Delta^2(k) (y^2-y^4) exp(-y^2) 
    //(with y=k/ksig)
    // First get the integral
    double result1, result2, error;

    gsl_function F;
    F.function = &(dCint);

    //Initialize parameter block
    ksig_params p;
    p.cptr = this;
    p.z=&z;

    F.params = &p;

    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F, 0., exp(1.), 0, 1e-6, 1000, w, &result1, &error); 
    gsl_integration_workspace_free (w);
    gsl_integration_workspace * u = gsl_integration_workspace_alloc (1000);
    gsl_integration_qagiu (&F, exp(1.), 0, 1e-6, 1000, u, &result2, &error); 
    gsl_integration_workspace_free (u);
    /*
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F, -1., 1., 0, 1e-6, 1000, w, &result2, &error); 
    gsl_integration_workspace_free (w);
    std::cout<<"I did the first integral"<<std::endl;
    gsl_integration_workspace * u = gsl_integration_workspace_alloc (1000);
    gsl_integration_qagil (&F, -1., 0, 1e-6, 1000, u, &result1, &error); 
    gsl_integration_workspace_free (u);
    gsl_integration_workspace * v = gsl_integration_workspace_alloc (1000);
    gsl_integration_qagiu (&F, 1., 0, 1e-6, 1000, v, &result3, &error); 
    gsl_integration_workspace_free (v);*/

    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);
    //std::cout<<"# Smith Curvature integral"<<(result1+result2+result3)<<std::endl;
    //return pow(3.+smithneff,2.)+4.*(result1+result2+result3);

    return pow(3.+smithneff,2.)+4.*(result1+result2);
    

}

/// Initialize the Smith et al. power spectrum variables.
double cosmology::init_Smith(double z)
{
   //Find ksigma first
    ksigma=getksigma(z);
    //std::cout<<"# "<<"ksigma ="<<ksigma<<std::endl;

    //Now initialize other variables, C, neff
    smithneff=getneff(z);
    smithC=getC(z);
    //std::cout<<"# "<<"neff ="<<smithneff<<std::endl;
    //std::cout<<"# "<<"C ="<<smithC<<std::endl;

    if(!takahashicorr){
        //Eq. C9 to C16
        an=1.4861+1.8369*smithneff+1.6762*pow(smithneff,2.)+0.7940*pow(smithneff,3.)+0.1670*pow(smithneff,4.)-0.6206*smithC;
        an=pow(10.,an);
        //std::cout<<"# "<<"an="<<an<<std::endl;

        bn=0.9463+0.9466*smithneff+0.3084*pow(smithneff,2.)-0.9400*smithC;
        bn=pow(10.,bn);
        //std::cout<<"# "<<"bn="<<bn<<std::endl;

        cn=-0.2807+0.6669*smithneff+0.3214*pow(smithneff,2.)-0.0793*smithC;
        cn=pow(10.,cn);
        //std::cout<<"# "<<"cn="<<cn<<std::endl;

        gamman=0.8649+0.2989*smithneff+0.1631*smithC;
        alphan=1.3884+0.3700*smithneff-0.1452*pow(smithneff,2.);
        betan =0.8291+0.9854*smithneff+0.3401*pow(smithneff,2.);
        //std::cout<<"# "<<"alphan="<<alphan<<std::endl;
        //std::cout<<"# "<<"betan="<<betan<<std::endl;
        //std::cout<<"# "<<"gamman="<<gamman<<std::endl;

        mun=-3.5442+0.1908*smithneff;
        mun=pow(10.,mun);
        //std::cout<<"# "<<"mun="<<mun<<std::endl;

        nun= 0.9589+1.2857*smithneff;
        nun=pow(10.,nun);
        //std::cout<<"# "<<"nun="<<nun<<std::endl;
    }else{
        double omegawz=Omegaw(z)*(1+w0+wa*z/(1.+z));
        //Eq. A6 to A13
        an=1.5222+2.8553*smithneff+2.3706*pow(smithneff,2.)+0.9903*pow(smithneff,3.)+0.2250*pow(smithneff,4.)-0.6038*smithC+0.1749*omegawz;
        an=pow(10.,an);
        //std::cout<<"# "<<"an="<<an<<std::endl;

        bn=-0.5642+0.5864*smithneff+0.5716*pow(smithneff,2.)-1.5474*smithC+0.2279*omegawz;
        bn=pow(10.,bn);
        //std::cout<<"# "<<"bn="<<bn<<std::endl;

        cn=0.3698+2.0404*smithneff+0.8161*pow(smithneff,2.)+0.5869*smithC;
        cn=pow(10.,cn);
        //std::cout<<"# "<<"cn="<<cn<<std::endl;

        gamman=0.1971-0.0843*smithneff+0.8460*smithC;
        alphan=fabs(6.0835+1.3373*smithneff-0.1959*pow(smithneff,2.)-5.5274*smithC);
        betan =2.0379-0.7354*smithneff+0.3157*pow(smithneff,2.)+1.2490*pow(smithneff,3.)+0.3980*pow(smithneff,4.)-0.1682*smithC;
        //std::cout<<"# "<<"alphan="<<alphan<<std::endl;
        //std::cout<<"# "<<"betan="<<betan<<std::endl;
        //std::cout<<"# "<<"gamman="<<gamman<<std::endl;

        mun=0.0;
        //std::cout<<"# "<<"mun="<<mun<<std::endl;

        nun= 5.2105+3.6902*smithneff;
        nun=pow(10.,nun);
        //std::cout<<"# "<<"nun="<<nun<<std::endl;
    }

    // f1,f2,f3 from Eq. C17 and C18
    double Om=pow(1.+z,3.)/pow(Eofz(z),2.)*Omega0;
    //double Om=Omega0;
    //std::cout<<"# "<<Om<<" Om "<<Omega0<<std::endl;

    if(fabs(Omega0+Omegal-1.)<1.e-4 && (w0==-1.) && (wa==0.0))
    {
        f1=pow(Om,-0.0307);
        f2=pow(Om,-0.0585);
        f3=pow(Om,0.0743);
    } else if(Omega0<1.&&Omegal==0.0 && (w0==-1.) && (wa==0.0))
    {
        f1=pow(Om,-0.0732);
        f2=pow(Om,-0.1423);
        f3=pow(Om,0.0725);
    }else
    {
	std::cout<<"ERROR: Smith et al power spectrum: This cosmology is not supported."<<std::endl;
	std::cout<<"ERROR: Smith et al power spectrum: Cosmology should be flat or without"<<std::endl; 
	std::cout<<"ERROR: Smith et al power spectrum: cosmological constant and with equation of state w=-1."<<std::endl;
	std::cout<<"ERROR: Try using the cosmic emulator from Heitman et al ."<<std::endl;
	exit(0);
    }
    //std::cout<<"# "<<f1<<"=f1 "<<f2<<"=ff2 "<<f3<<"=f3 "<<std::endl;

    bool_initSmith=true;

    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    }

    return 1;

}

/// Delta^2(k,z) a'la Smith et al. 2003
/// kpassed should be supplied as h Mpc^{-1}.
double cosmology::Delta2NL_S(double k, double z)
{
    if (!bool_initPS) 
    {
	if(verbose){
        std::cout<<"# "<<" Power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = initPS_EH(); 
	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }
    if (!bool_initSmith) 
    {
	if(verbose){
        std::cout<<"# "<<" Nonlinear power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = init_Smith(z); 
	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }
    double d2lk=Delta2_L(k,z);

    double y=k/ksigma;
    double fy=y/4.0+pow(y,2.)/8.0;
    //Eq. C2
    double Delta2_Q=d2lk*pow(1.+d2lk,betan)/(1.+alphan*d2lk)*exp(-fy);

    // Eq. C4
    double Delta2_Hp=an*pow(y,3.*f1)/(1.+bn*pow(y,f2)+pow(cn*f3*y,3.-gamman));
    // Eq. C3
    double Delta2_H=Delta2_Hp/(1.+mun/y+nun/pow(y,2.));

    return Delta2_H+Delta2_Q;
}

/// P(k,z) a'la Smith et al. 2003
/// kpassed should be supplied as h Mpc^{-1}.
/// Result is in the units of (h^{-1} Mpc)^3
double cosmology::PkNL_S(double kpassed, double z)
{
    if (!bool_initPS) 
    {
	if(verbose){
        std::cout<<"# "<<" Power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = initPS_EH(); 
	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }
    if (!bool_initSmith) 
    {
	if(verbose){
        std::cout<<"# "<<" Nonlinear power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = init_Smith(z); 
	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }
    double k = kpassed;

    double res=Delta2NL_S(k,z);

    return res*(2*M_PI*M_PI)/(pow(k,3.));
}

/// Delta^2(k,z) a'la Eisenstein and Hu (linear) & Smith et al. 2003 (nonlinear)
/// kpassed should be supplied as h Mpc^{-1}.
void cosmology::Delta2_EH_S(double k, double z, double &lin, double &nonlin)
{
    if (!bool_initPS) 
    {
	if(verbose){
        std::cout<<"# "<<" Power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = initPS_EH(); 
	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }
    if (!bool_initSmith) 
    {
	if(verbose){
        std::cout<<"# "<<" Nonlinear power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = init_Smith(z); 
	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }
    double d2lk=Delta2_EH(k,z);
    //std::cout<<"#d2lk "<<d2lk<<std::endl;

    double y=k/ksigma;
    double fy=y/4.0+pow(y,2.)/8.0;
    //Eq. C2
    double Delta2_Q=d2lk*pow(1.+d2lk,betan)/(1.+alphan*d2lk)*exp(-fy);

    // Eq. C4
    double Delta2_Hp=an*pow(y,3.*f1)/(1.+bn*pow(y,f2)+pow(cn*f3*y,3.-gamman));
    // Eq. C3
    double Delta2_H=Delta2_Hp/(1.+mun/y+nun/pow(y,2.));

    // Eq. C1
    nonlin=Delta2_H+Delta2_Q;
    lin=d2lk;
}

/// P(k,z) a'la Eisenstein and Hu (linear) & Smith et al. 2003 (nonlinear)
/// kpassed should be supplied as h Mpc^{-1}.
void cosmology::Pk_EH_S(double k,double z,double &lin ,double &nonlin)
{
    if (!bool_initPS) 
    {
	if(verbose){
        std::cout<<"# "<<" Power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = initPS_EH(); 
	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }
    if (!bool_initSmith) 
    {
	if(verbose){
        std::cout<<"# "<<" Nonlinear power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = init_Smith(z); 
	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }
    double d2lin,d2nl;
    Delta2_EH_S(k,z,d2lin,d2nl);
    lin=d2lin*(2*M_PI*M_PI)/(pow(k,3.));
    nonlin=d2nl*(2*M_PI*M_PI)/(pow(k,3.));

}

/// Fourier transforms

double dxi_NL(double x, void* params)
{
    /// Int_{0}^{infty} Delta(x/r)/(x*x)*sin(x) dx
    xi_params c1 = *(xi_params *)params;

    cosmology *c2;
    double *r;
    double *z;

    c2=c1.cptr;
    r=c1.r;
    z=c1.z;

    double arg=x/(*r);
    double res=(*c2).Delta2_NL_num(arg,*z)/(x*x);

    //std::cout<<"Check "<<x<<" "<<res<<" "<<*r<<std::endl;

    return res;
}

/// Non linear correlation function. Fourier transform the power spectrum
double cosmology::xi_NL(double r,double z)
{
    /// Setup the integral and execute
    double result,error;

    gsl_function F;
    F.function = &(dxi_NL);

    //Initialize parameter block
    xi_params p;
    p.cptr = this;
    p.r=&r;
    p.z=&z;

    F.params = &p;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (100);
    gsl_integration_workspace * cw = gsl_integration_workspace_alloc (100);
    gsl_integration_qawo_table * qt = gsl_integration_qawo_table_alloc (1.0,1.0,GSL_INTEG_SINE,100);
    gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
    int status=gsl_integration_qawf (&F, 1.e-30, 1.e-4, 100, w, cw, qt, &result, &error); 
    if(status!=0.0) {
	std::cout<<"Error occurred in xi NL qawf: "<<r<<std::endl;
	result=0.0;
    }

    gsl_integration_workspace_free (w);
    gsl_integration_workspace_free (cw);
    gsl_integration_qawo_table_free (qt);
    gsl_set_error_handler(oldhand);

    //std::cout<<r<<" "<<result<<" "<<error<<std::endl;

    return result;
}

double dxi_L(double x, void* params)
{
    /// Int_{0}^{infty} Delta(x/r)/(x*x)*sin(x) dx
    xi_params c1 = *(xi_params *)params;

    cosmology *c2;
    double *r;
    double *z;

    c2=c1.cptr;
    r=c1.r;
    z=c1.z;

    double arg=x/(*r);
    double res=(*c2).Delta2_L_num(arg,*z)/(x*x);

    //std::cout<<"Check "<<x<<" "<<res<<" "<<*r<<std::endl;

    return res;
}

/// Linear correlation function. Fourier transform the power spectrum
double cosmology::xi_L(double r,double z)
{
    //std::cout<<"Check: "<<r<<std::endl;
    /// Setup the integral and execute
    double result,error;
    result=0.0;
    //if(r<1.0e-10) return result;

    gsl_function F;
    F.function = &(dxi_L);

    //Initialize parameter block
    xi_params p;
    p.cptr = this;
    p.r=&r;
    p.z=&z;

    F.params = &p;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (100);
    gsl_integration_workspace * cw = gsl_integration_workspace_alloc (100);
    gsl_integration_qawo_table * qt = gsl_integration_qawo_table_alloc (1.0,1.0,GSL_INTEG_SINE,100);
    //gsl_set_error_handler_off();
    gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
    int status=gsl_integration_qawf (&F, 1.0e-30, 1.e-4, 100, w, cw, qt, &result, &error); 
    gsl_set_error_handler(oldhand);
    if(status!=0.0) {
	std::cout<<"Error occurred in xi_L qawf: "<<r<<std::endl;
	result=0.0;
    }

    gsl_integration_workspace_free (w);
    gsl_integration_workspace_free (cw);
    gsl_integration_qawo_table_free (qt);

    ///std::cout<<"Check: "<<r<<" "<<result<<std::endl;

    return result;
}

/// Linear correlation function. Fourier transform the power spectrum
double cosmology::xi_NL_old(double r,double z)
{
    /// Setup the integral and execute
    double result;

    result=0.0;
    int Nloop=1000;

    for (int i=0;i<Nloop;i++){
        double xc=i*2.0*M_PI;
        for (int j=0;j<N0_2p;j++){
            double x=xc+x0_2p[j];
            result=result+w0_2p[j]*Delta2_NL(x/r,z)/(x*x)*sin(x);
        }
    }

    return result;
}

/// Linear correlation function. Fourier transform the power spectrum
double cosmology::xi_L_old(double r,double z)
{
    /// Setup the integral and execute
    double result;

    result=0.0;
    int Nloop=240;

    for (int i=0;i<Nloop;i++){
	double xc=i*2*M_PI;
	for (int j=0;j<N0_2p;j++){
	    double x=xc+x0_2p[j];
	    result=result+w0_2p[j]*Delta2_L(x/r,z)/(x*x)*sin(x);
	}
    }

    return result;
}

/// Power spectra numerical implementation. Initializing spline. 
void cosmology::init_powerspectra_L()
{
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for linear power spectra was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance

    double xx[Npower],yy[Npower];
    for (int i=0;i<Npower;i++)
    {
        double k=kmin+i/(Npower-1.)*(kmax-kmin);
        xx[i]=k;
	k=pow(10.0,k);
	/// Improvement possible here
        yy[i]=log10(Delta2_L(k,0.0));
    }
    PSL0_acc = gsl_interp_accel_alloc ();
    PSL0_spline = gsl_spline_alloc (gsl_interp_cspline, Npower);
    gsl_spline_init (PSL0_spline,xx,yy,Npower);
    PSL0_dlow=(yy[2]-yy[0])/(xx[2]-xx[0]);
    PSL0_xlow=xx[1]; PSL0_ylow=yy[1];
    PSL0_dhigh=(yy[Npower-3]-yy[Npower-1])/(xx[Npower-3]-xx[Npower-1]);
    PSL0_xhigh=xx[Npower-2]; PSL0_yhigh=yy[Npower-2];

    bool_init_PSL0=true;
    if(verbose){
    std::cout<<"# "<<"The spline for numerical calculation of linear power spectra is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
       
}


/// Power spectra numerical implementation. Initializing spline. 
void cosmology::init_powerspectra_NL(double z)
{
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for nonlinear power spectra was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance

    double xx[Npower],zz[Npower];
    for (int i=0;i<Npower;i++)
    {
        double k=kmin+i/(Npower-1.)*(kmax-kmin);
        xx[i]=k;
        k=pow(10.0,k);
        double ycorr=1.0;
        if(!takahashicorr && peacockcorr){
            ycorr=k/10;
            ycorr=(1+2*ycorr*ycorr)/(1+ycorr*ycorr);
        }
	/// Improvement possible here
        zz[i]=log10(Delta2_NL(k,z)*ycorr);
    }

    PSNL_acc = gsl_interp_accel_alloc ();
    PSNL_spline = gsl_spline_alloc (gsl_interp_cspline, Npower);
    gsl_spline_init (PSNL_spline,xx,zz,Npower);
    PSNL_dlow=(zz[2]-zz[0])/(xx[2]-xx[0]);
    PSNL_xlow=xx[1]; PSNL_ylow=zz[1];
    PSNL_dhigh=(zz[Npower-3]-zz[Npower-1])/(xx[Npower-3]-xx[Npower-1]);
    PSNL_xhigh=xx[Npower-2]; PSNL_yhigh=zz[Npower-2];

    bool_init_PSNL=true;
    if(verbose){
    std::cout<<"# "<<"The spline for numerical calculation of power spectra is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
       
}

/// Calculate Delta^2L numerically with k in h Mpc^{-1}
double cosmology::Delta2_L_num(double k, double z)
{
    //
    if(!bool_init_PSL0)
    {
        init_powerspectra_L();
    }
    double logk=log10(k);
    double result;
    if(logk>kmin&&logk<kmax)
    {
        result=gsl_spline_eval (PSL0_spline,logk, PSL0_acc);
    } else if (logk<=kmin)
    {
        //std::cout<<"# "<<"PSL_spline: Interpolation not possible k not within range "<<k<<std::endl;
        //exit(0);
        //result=log10(Delta2_L(k,z));
        //result=1.0e-30;
	//Extrapolate now
	result=PSL0_ylow+PSL0_dlow*(logk-PSL0_xlow);
    } else if (logk>=kmax)
    {
	result=PSL0_yhigh+PSL0_dhigh*(logk-PSL0_xhigh);
    }

    result = pow(10.,result);
    if(z!=z_glob){
        result=result*pow(growthfactor_num(z),2.);
    }else{
        result=result*pow(gf_glob,2.);
    }
    return result;

}

/// Calculate Delta^2L numerically with k in h Mpc^{-1}
double cosmology::Delta2_NL_num(double k, double z)
{
    if(z!=z_glob){
	setnew_z(z);
    }
    //
    if(!bool_init_PSNL)
    {
        init_powerspectra_NL(z);
    }
    double logk=log10(k);
    double result;
    if(logk>kmin&&logk<kmax)
    {
        result=gsl_spline_eval (PSNL_spline,logk, PSNL_acc);
    } else if (logk <= kmin)
    {
        //std::cout<<"# "<<"PSL_spline: Interpolation not possible k not within range "<<k<<std::endl;
        //exit(0);
        //result=log10(Delta2_NL(k,z));
        //result=1.0e-30;
	result=PSNL_ylow+PSNL_dlow*(logk-PSNL_xlow);
    } else if (logk>=kmax)
    {
	result=PSNL_yhigh+PSNL_dhigh*(logk-PSNL_xhigh);
    }
    result = pow(10.,result);
    return result;

}

/// linear correlation numerical implementation. Initializing spline. 
void cosmology::init_xi_L()
{
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for linear correlations was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance

    double xx[Nxi],yy[Nxi];
    for (int i=0;i<Nxi;i++)
    {
        double r=rmin+i/(Nxi-1.)*(rmax-rmin);
        xx[i]=r;
	r=pow(10.0,r);
	/// Store 1+xi as that is a positive number
        yy[i]=log10(1.+xi_L(r,0.0));
    }
    xiL0_acc = gsl_interp_accel_alloc ();
    xiL0_spline = gsl_spline_alloc (gsl_interp_cspline, Nxi);
    gsl_spline_init (xiL0_spline,xx,yy,Nxi);
    xiL0_dlow=(yy[2]-yy[0])/(xx[2]-xx[0]);
    xiL0_xlow=xx[1]; xiL0_ylow=yy[1];
    xiL0_dhigh=(yy[Nxi-3]-yy[Nxi-1])/(xx[Nxi-3]-xx[Nxi-1]);
    xiL0_xhigh=xx[Nxi-2]; xiL0_yhigh=yy[Nxi-2];

    bool_init_xiL0=true;
    if(verbose){
    std::cout<<"# "<<"The spline for linear correlation is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
       
}

/// Calculate xi_L numerically with r in h^{-1} Mpc
double cosmology::xi_L_num(double r, double z)
{
    //
    if(!bool_init_xiL0)
    {
        init_xi_L();
    }
    double logr=log10(r);
    double result;
    if(logr>=rmin&&logr<=rmax)
    {
        result=gsl_spline_eval (xiL0_spline,logr, xiL0_acc);
    } else if (logr<rmin)
    {
        //std::cout<<"# "<<"xiL_spline: r not within range "<<r<<std::endl;
        //result=xi_L(r,z);
	result=xiL0_ylow+xiL0_dlow*(logr-xiL0_xlow);
    } else if (logr>rmax)
    {
	result=xiL0_yhigh+xiL0_dhigh*(logr-xiL0_xhigh);
    }
    result=pow(10.0,result)-1.;
    if(z!=z_glob){
        result=result*pow(growthfactor_num(z),2.);
    }else{
        result=result*pow(gf_glob,2.);
    }
    return result;

}

/// Nonlinear correlation numerical implementation. Initializing spline. 
void cosmology::init_xi_NL(double z)
{
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for non-linear correlations was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance

    double xx[Nxi],yy[Nxi];
    for (int i=0;i<Nxi;i++)
    {
        double r=rmin+i/(Nxi-1.)*(rmax-rmin);
        xx[i]=r;
	r=pow(10.0,r);
	/// Store 1+xi as that is a positive number
        double xinl=xi_NL(r,z);
        yy[i]=log10(1.+xinl);
        //std::cout<<r<<" "<<xinl<<std::endl;
    }
    xiNL_acc = gsl_interp_accel_alloc ();
    xiNL_spline = gsl_spline_alloc (gsl_interp_cspline, Nxi);
    gsl_spline_init (xiNL_spline,xx,yy,Nxi);
    xiNL_dlow=(yy[2]-yy[0])/(xx[2]-xx[0]);
    xiNL_xlow=xx[1]; xiNL_ylow=yy[1];
    xiNL_dhigh=(yy[Nxi-3]-yy[Nxi-1])/(xx[Nxi-3]-xx[Nxi-1]);
    xiNL_xhigh=xx[Nxi-2]; xiNL_yhigh=yy[Nxi-2];
    bool_init_xiNL=true;

    /// Also initialize rzetamax here!
    setrzeta(z);
    //exit(100);

    if(verbose){
    std::cout<<"# "<<"The spline for nonlinear correlation is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
       
}

/// Calculate xi_NL numerically with r in h^{-1} Mpc
double cosmology::xi_NL_num(double r, double z)
{
    if(z!=z_glob){
	setnew_z(z);
    }
    //
    if(!bool_init_xiNL)
    {
        init_xi_NL(z);
    }
    double logr=log10(r);
    double result;
    if(logr>=rmin&&logr<=rmax)
    {
        result=gsl_spline_eval (xiNL_spline,logr, xiNL_acc);
    } else if (logr<rmin)
    {
        //std::cout<<"# "<<"xiNL_spline: r not within range "<<r<<std::endl;
        //result=xi_NL(r,z);
	result=xiNL_ylow+xiNL_dlow*(logr-xiNL_xlow);
    } else if (logr>rmax)
    {
	result=xiNL_yhigh+xiNL_dhigh*(logr-xiNL_xhigh);
	//result=0.0;
    }
    result=pow(10.0,result)-1.;
    return result;

}

double dPktest_L(double x, void* params)
{
    /// 4 pi / (k**3.0)*Int_{0}^{infty} x*sin(x)*xi(x/k) dx
    pk_params c1 = *(pk_params *)params;

    cosmology *c2;
    double *k;
    double *z;

    c2=c1.cptr;
    k=c1.k;
    z=c1.z;

    double arg=x/(*k);
    double res=(*c2).xi_L_num(arg,*z)*x;

    //std::cout<<"Check "<<x<<" "<<res<<" "<<*k<<std::endl;

    return res;
}

/// Linear correlation function. Fourier transform the power spectrum
double cosmology::Pktest_L(double k,double z)
{
    /// Setup the integral and execute
    double result,error;

    gsl_function F;
    F.function = &(dPktest_L);

    //Initialize parameter block
    pk_params p;
    p.cptr = this;
    p.k=&k;
    p.z=&z;

    F.params = &p;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_workspace * cw = gsl_integration_workspace_alloc (1000);
    gsl_integration_qawo_table * qt = gsl_integration_qawo_table_alloc (1.0,1.0,GSL_INTEG_SINE,1000);
    int status=gsl_integration_qawf (&F, 1.e-30, 1.e-2, 1000, w, cw, qt, &result, &error); 
    if(status!=0.0) {
	std::cout<<"Error occurred in qawf: "<<k<<std::endl;
	result=0.0;
    }

    gsl_integration_workspace_free (w);
    gsl_integration_workspace_free (cw);
    gsl_integration_qawo_table_free (qt);

    return result*4*M_PI/pow(k,3.0);
}

double dPktest_NL(double x, void* params)
{
    /// 4 pi / (k**3.0)*Int_{0}^{infty} x*sin(x)*xi(x/k) dx
    pk_params c1 = *(pk_params *)params;

    cosmology *c2;
    double *k;
    double *z;

    c2=c1.cptr;
    k=c1.k;
    z=c1.z;

    //std::cout<<"Check "<<x<<" "<<" "<<*k<<std::endl;
    double arg=x/(*k);
    double res=(*c2).xi_NL_num(arg,*z)*x; //*4.*M_PI;

    //std::cout<<"Check out "<<x<<" "<<res<<" "<<*k<<std::endl;

    return res;
}

/// Non Linear correlation function. Fourier transform the power spectrum
double cosmology::Pktest_NL(double k,double z)
{
    /// Setup the integral and execute
    double result,error;

    gsl_function F;
    F.function = &(dPktest_NL);

    //Initialize parameter block
    pk_params p;
    p.cptr = this;
    p.k=&k;
    p.z=&z;

    F.params = &p;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_workspace * cw = gsl_integration_workspace_alloc (1000);
    gsl_integration_qawo_table * qt = gsl_integration_qawo_table_alloc (1.0,1.0,GSL_INTEG_SINE,1000);
    int status=gsl_integration_qawf (&F, 1.e-30, 1.e-2, 1000, w, cw, qt, &result, &error); 
    if(status!=0.0) {
	std::cout<<"Error occurred in qawf: "<<k<<std::endl;
	result=0.0;
    }

    gsl_integration_workspace_free (w);
    gsl_integration_workspace_free (cw);
    gsl_integration_qawo_table_free (qt);

    return result*4.0*M_PI/pow(k,3.0);
}

double dPktest_zetaNL(double x, void* params)
{
    /// 4 pi / (k**3.0)*Int_{0}^{infty} x*sin(x)*xi(x/k) dx
    pk_params c1 = *(pk_params *)params;

    cosmology *c2;
    double *k;
    double *z;

    c2=c1.cptr;
    k=c1.k;
    z=c1.z;

    //std::cout<<"Check "<<x<<" "<<" "<<*k<<std::endl;
    double arg=x/(*k);
    double xiNL=(*c2).xi_NL_num(arg,*z);
    double res=xiNL*pow(1.+1.17*xiNL,1.49)/pow(1.+0.69*xiNL,2.09)*x;

    //std::cout<<"Check out "<<x<<" "<<res<<" "<<*k<<std::endl;

    return res;
}

/// Non Linear correlation function. Fourier transform the power spectrum
double cosmology::Pktest_zetaNL(double k,double z)
{
    /// Setup the integral and execute
    double result,error;

    gsl_function F;
    F.function = &(dPktest_zetaNL);

    //Initialize parameter block
    pk_params p;
    p.cptr = this;
    p.k=&k;
    p.z=&z;

    F.params = &p;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_workspace * cw = gsl_integration_workspace_alloc (1000);
    gsl_integration_qawo_table * qt = gsl_integration_qawo_table_alloc (1.0,1.0,GSL_INTEG_SINE,1000);
    int status=gsl_integration_qawf (&F, 1.e-30, 1.e-2, 1000, w, cw, qt, &result, &error); 
    if(status!=0.0) {
	std::cout<<"Error occurred in qawf: "<<k<<std::endl;
	result=0.0;
    }

    gsl_integration_workspace_free (w);
    gsl_integration_workspace_free (cw);
    gsl_integration_qawo_table_free (qt);

    return result*4.0*M_PI/pow(k,3.0);
}

/// Set new redshift
/// Linear power spectrum is not required to be reinited
void cosmology::setnew_z(double z)
{
    bool_initSmith=false;

    if(bool_init_PSNL)
    {
	gsl_interp_accel_free(PSNL_acc);
	gsl_spline_free(PSNL_spline);
	bool_init_PSNL=false;
    }

    if(bool_init_xiNL)
    {
	gsl_interp_accel_free(xiNL_acc);
	gsl_spline_free(xiNL_spline);
	bool_init_xiNL=false;
    }

    if(bool_init_xiNL_bar)
    {
	gsl_interp_accel_free(xiNL_bar_acc);
	gsl_spline_free(xiNL_bar_spline);
	bool_init_xiNL_bar=false;
    }

    if(bool_init_xiNL_barbar)
    {
	gsl_interp_accel_free(xiNL_barbar_acc);
	gsl_spline_free(xiNL_barbar_spline);
	bool_init_xiNL_barbar=false;
    }

    if(bool_init_xiL_bar)
    {
	gsl_interp_accel_free(xiL_bar_acc);
	gsl_spline_free(xiL_bar_spline);
	bool_init_xiL_bar=false;
    }

    if(bool_init_xiL_barbar)
    {
	gsl_interp_accel_free(xiL_barbar_acc);
	gsl_spline_free(xiL_barbar_spline);
	bool_init_xiL_barbar=false;
    }

    init_Tink=false;
    z_glob=z;
    gf_glob=growthfactor_num(z);
    if(verbose){
        std::cout<<"# NEW REDSHIFT SET"<<std::endl;
    }

    if(bool_initQk)
    {
        if(verbose){
                std::cout<<"# Deleting Qk's"<<std::endl;
        }
        delete []Qk1;
        delete []Qk2;
    }
    bool_initQk=false;
}

double findrz(double x, void *params){
    rz_params c1 = *(rz_params *) params;
    cosmology *c2;
    double *z;
    double *tgt;

    c2=c1.cptr;
    z=c1.z;
    tgt=c1.tgt;

    double xinl=(*c2).xi_NL_num(x,*z);
    double rzeta=pow(1.+1.17*xinl,1.49)/pow(1.+0.69*xinl,2.09);

    return xinl*rzeta-(*tgt);

}

double cosmology::findrzfn(double x, double tgt, double z){

    double xinl=xi_NL_num(x,z);
    double rzeta=pow(1.+1.17*xinl,1.49)/pow(1.+0.69*xinl,2.09);

    return xinl*rzeta-tgt;

}


/// Set the value of rzeta
void cosmology::setrzeta(double z){
    /// Gets called by the init_xi_NL routine
    
    // Finds the radius where rzeta*xiNL=10.0**(xiNLzetamax)
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;


    gsl_function F;
    rz_params p;
    p.cptr = this;
    p.z=&z;
    double tgt=pow(10.,xiNLzetamax);
    p.tgt=&tgt;

    double r_lo =0.1, r_hi =1.5;
    bool bracket=false;
    while(findrzfn(r_hi,tgt,z)>0.0){
        r_lo=r_hi;
        r_hi=r_hi*1.2;
        bracket=true;
        if(verbose) std::cout<<"# In first bracket loop,"<<r_lo<<" "<<r_hi<<std::endl;
    }

    if(!bracket){
        r_lo=r_hi;
        while(findrzfn(r_lo,tgt,z)<0.0){
            r_hi=r_lo;
            r_lo=r_lo/1.2;
            bracket=true;
            if(verbose) std::cout<<"# In second bracket loop,"<<r_lo<<" "<<r_hi<<std::endl;
        }
    }

    if(!bracket){
        /// Problem in bracketing! Try to exit gracefully!
        std::cout<<"# Problem in bracketing in setrzeta"<<std::endl;
        exit(1001);
    }

    F.function = &findrz;
    F.params = &p;
   
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, r_lo, r_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        r_lo = gsl_root_fsolver_x_lower (s);
        r_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (r_lo, r_hi,0, 1e-6);

        if (status == GSL_SUCCESS)
        {
            if(verbose) std::cout<<"# "<<"rzeta:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    /// Set value of r and zetamax
    zeta_rmax=res;
    double xinl=xi_NL_num(zeta_rmax,z);
    zetamax=pow(1.+1.17*xinl,1.49)/pow(1.+0.69*xinl,2.09);
    if(verbose) std::cout<<"# xinl "<<xinl<<" "<<zeta_rmax<<" "<<zetamax<<std::endl;
}

double dwpnl(double x, void* params)
{
    /// Wp(R) = 2.* R Int_{0}^{xmax} xi( \sqrt(x^2+1) R ) dx
    /// R is multiplied outside the integral
    projwp_params c1 = *(projwp_params *)params;

    cosmology *c2;
    double *R;
    double *z;

    c2=c1.cptr;
    R=c1.R;
    z=c1.z;

    double arg=sqrt(x*x+1.)*(*R);
    double res=(*c2).xi_NL_num(arg,*z);

    return res;
}

double cosmology::wpnl(double z, double rad,double projmax){

    /// QAGW integral
    ///Int_{0}^{ymax} dy 2  R \xi_gg( \sqrt(y^2+1) R)
    /// ymax=pimax/R
    double result, error;
    double ymax=(projmax/rad);

    gsl_function F;
    F.function = &(dwpnl);

    projwp_params p;
    p.cptr = this;
    p.R = &(rad);
    p.z = &z;

    F.params = &p;

    int stat;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    stat=gsl_integration_qag (&F, 0.0, ymax, 0, 1.e-3, 1000, 5, w, &result, &error);
    gsl_integration_workspace_free (w);

        result=2.0*result*rad;
        return result;
}

double dwpl(double x, void* params)
{
    /// Wp(R) = 2.* R Int_{0}^{xmax} xi( \sqrt(x^2+1) R ) dx
    /// R is multiplied outside the integral
    projwp_params c1 = *(projwp_params *)params;

    cosmology *c2;
    double *R;
    double *z;

    c2=c1.cptr;
    R=c1.R;
    z=c1.z;

    double arg=sqrt(x*x+1.)*(*R);
    double res=(*c2).xi_L_num(arg,*z);

    return res;
}

double cosmology::wpl(double z, double rad, double projmax){

	/// QAGW integral
	///Int_{0}^{ymax} dy 2  R \xi_gg( \sqrt(y^2+1) R)
	/// ymax=pimax/R
	double result, error;
	double ymax=(projmax/rad);

	gsl_function F;
	F.function = &(dwpl);

	projwp_params p;
	p.cptr = this;
	p.R = &(rad);
	p.z = &z;

	F.params = &p;

	int stat;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	stat=gsl_integration_qag (&F, 0.0, ymax, 0, 1.e-3, 1000, 5, w, &result, &error);
	gsl_integration_workspace_free (w);

        result=2.0*result*rad;
        return result;
}
/// Friend function for calculating Qk2
double dQk(double x, void * params){
    qk_params c1 = *(qk_params *)params;

    cosmology *c2;
    double *k;
    double *z;
    double *xmax;
    double zetamax;
    double zetarmax;
    double opt;

    c2=c1.cptr;
    k=c1.k;
    z=c1.z;
    xmax=c1.xmax;
    opt=c1.opt;
    double r=x/(*k);
    double xinl=(*c2).xi_NL_num(r,*z);
    zetamax=(*c2).getzetamax();
    zetarmax=(*c2).getzeta_rmax();
    double kfac=(4.0*M_PI)/(*k);

    double res=0.0;
    if(r<(zetarmax)){
        res=kfac*x*(zetamax)*xinl;
    }else{
        double zeta=pow(1.+1.17*xinl,1.49)/pow(1.+0.69*xinl,2.09);
        res=kfac*x*zeta*xinl;
    }

    if(opt==1){
        if(res<0.0 && x< *xmax){
            *xmax=x;
        }
        if(res<0.0) res=0;
    }else if(opt==2){
        res=-res;
    }else if(opt==3){
        printf("%e %e %e\n",x,res,r);
    }
    return res;

}

/**/
void cosmology::initQk(double z, double r200[]){
    if(z!=z_glob){
	setnew_z(z);
    }
    // Allocate the arrays for Qk
    Qk1=new double[kbins*N9_16*N9_16];
    Qk2=new double[kbins*N9_16*N9_16];

    double karr[kbins];
    double hod_kdiff=(hod_kmax-hod_kmin)/(kbins-1.);
    for(int jk=0;jk<kbins;jk++){
        karr[jk]=hod_kmin+jk*hod_kdiff;
        karr[jk]=pow(10.,karr[jk]);
        double k=karr[jk];
        if(k>pow(10.0,-2.0))
        {
            const static double qk_rmin=-2.0;
            const static double qk_rmax=1.0;
            const static int Nqk=100;
            double qk_xx[Nqk];
            double qk_yy[Nqk];
            gsl_interp_accel *qk_acc;
            gsl_spline *qk_spline;

            for(int i=0;i<Nqk;i++)
            {
                qk_xx[i]=pow(10.,qk_rmin+(qk_rmax-qk_rmin)*i/(Nqk-1.));
                qk_yy[i]=0.0;
            }

            /// Initialize xi_NL in case it is not inited yet. This sets
            /// up zetarmax too!
            if(!bool_init_xiNL)
            {
                init_xi_NL(z);
            }

            /// Set up the oscillation table
            gsl_integration_qawo_table * qt0 = gsl_integration_qawo_table_alloc (1.0,k*(qk_xx[1]-qk_xx[0]),GSL_INTEG_SINE,1000);

            /// Define integrand
                gsl_function F;
                F.function=&(dQk);

                double xmax;
                qk_params p;
                p.cptr =this;
                p.k=&k;
                p.z=&z;
                p.xmax=&xmax;

                F.params=&p;

            /// Attempt to parallelize this loop
            for(int i=0;i<Nqk;i++)
            {
                //xmax=1.0E30;
                p.opt=0;
                double result,error;
                result=error=0.0;

                gsl_integration_workspace * w =
                gsl_integration_workspace_alloc (1000);
                gsl_integration_workspace * cw =
                gsl_integration_workspace_alloc (1000);
                gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
                int status=gsl_integration_qawf(&F,k*qk_xx[i],1.e-3, 1000, w, cw, qt0, &result, &error); 
                gsl_set_error_handler(oldhand);
                if(status!=0 && verbose){
                    //printf("# Problems in convergence: initQk, k=%e z=%e rmin=%e, result=%e error=%e \n", k, z, qk_xx[i],result,error);
                    //printf("# Problems in convergence: initQk, xmax=%e \n", xmax);
                    printf("# Problems in convergence: initQk \n");
                    //exit(101);
                }
                gsl_integration_workspace_free (w);
                gsl_integration_workspace_free (cw);

                qk_yy[i]=result;

            }

            gsl_integration_qawo_table_free(qt0);
             
            double kfac=1.0/k/k;
            double kfac2=(4.0*M_PI)/k/k/k;
            for(int i=0;i<Nqk;i++){
                qk_xx[i]=log10(qk_xx[i]);
                qk_yy[i]=qk_yy[i];
            }

            qk_acc = gsl_interp_accel_alloc ();
            qk_spline = gsl_spline_alloc (gsl_interp_cspline, Nqk);
            gsl_spline_init (qk_spline,qk_xx,qk_yy,Nqk);

            /// Fill Qk1 and Qk2 using the r200 array
            for (int i=0;i<N9_16;i++){
                for (int j=0;j<=i;j++){
                    double rsho=(r200[j]>r200[i])?r200[j]:r200[i];
                    double omga=k*rsho;
                    rsho=log10(rsho);
                    Qk1[(jk*N9_16+i)*N9_16+j]=-kfac2*( sin(omga) - omga*cos(omga)  );
                    Qk2[(jk*N9_16+i)*N9_16+j]=kfac*gsl_spline_eval(qk_spline,rsho,qk_acc);
                    Qk1[(jk*N9_16+j)*N9_16+i]=Qk1[(jk*N9_16+i)*N9_16+j];
                    Qk2[(jk*N9_16+j)*N9_16+i]=Qk2[(jk*N9_16+i)*N9_16+j];
                }
            }

            /// Free the spline
            gsl_interp_accel_free(qk_acc);
            gsl_spline_free(qk_spline);

        }else{
            /// Fill Qk1 and Qk2 using the r200 array
            double kfac=k*k*k/(2.0*M_PI*M_PI);
            double kfac2=(4.0*M_PI)/k/k/k;
            for (int i=0;i<N9_16;i++){
                for (int j=0;j<=i;j++){
                    double rsho=(r200[j]>r200[i])?r200[j]:r200[i];
                    double omga=k*rsho;
                    Qk1[(jk*N9_16+i)*N9_16+j]=-kfac2*( sin(omga) - omga*cos(omga)  );
                    Qk2[(jk*N9_16+i)*N9_16+j]=Delta2_L_num(k,z)/kfac;
                    Qk1[(jk*N9_16+j)*N9_16+i]=Qk1[(jk*N9_16+i)*N9_16+j];
                    Qk2[(jk*N9_16+j)*N9_16+i]=Qk2[(jk*N9_16+i)*N9_16+j];
                }
            }
        }
    }

    bool_initQk=true;

}

/// First non-linear part for the redshift space distortions
double cosmology::xiNL_bar_num(double r,double z){
    double result=0.0;
    if(z!=z_glob){
	setnew_z(z);
    }
    if(!bool_init_xiNL_bar){
        init_xiNL_bar(z);
    }
    double logr=log10(r);
    result=gsl_spline_eval (xiNL_bar_spline,logr, xiNL_bar_acc);
    return pow(10.,result)-1.;
}

void cosmology::init_xiNL_bar(double z)
{
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for numerical calculation of xiNLbar was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance
    double rrmin=-3.0;
    double rrmax=3.0;

    double xx[Nxibar],yy[Nxibar];
    for (int i=0;i<Nxibar;i++)
    {
        double r=rrmin+i/(Nxibar-1.)*(rrmax-rrmin);
        xx[i]=r;
        r=pow(10.,r);
        yy[i]=log10(1.+xiNL_bar(r,z));
    }
    xiNL_bar_acc = gsl_interp_accel_alloc ();
    xiNL_bar_spline = gsl_spline_alloc (gsl_interp_cspline, Nxibar);
    gsl_spline_init (xiNL_bar_spline,xx,yy,Nxibar);

    bool_init_xiNL_bar=true;
    if(verbose){
    std::cout<<"# "<<"The spline for numerical calculation of xiNLbar is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
}

double dxinlbar(double x, void* params)
{
    ///Int_{0}^{1.0} dy y^2 \xi_gg( y*r )
    xi_params c1 = *(xi_params *)params;

    cosmology *c2;
    double *r;
    double *z;

    c2=c1.cptr;
    r=c1.r;
    z=c1.z;

    double arg=x*(*r);
    double res=(*c2).xi_NL_num(arg,*z);
    res=res*x*x;

    return res;
}

double cosmology::xiNL_bar(double r,double z){
	/// QAGW integral
	///Int_{0}^{1.0} dy y^2 \xi_gg( y*r )

	double result, error;
	gsl_function F;
	F.function = &(dxinlbar);

	xi_params p;
	p.cptr = this;
	p.r = &(r);
	p.z = &z;

	F.params = &p;

	int stat;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	//stat=gsl_integration_qag (&F, -15.0, 0.0, 0, 1.e-4, 1000, 5, w, &result, &error);
	stat=gsl_integration_qags (&F, 0.0, 1.0, 0, 1.e-3, 1000, w, &result, &error);
	gsl_integration_workspace_free (w);

        return result*3.;
}

double cosmology::xiNL_barbar_num(double r,double z){
    double result=0.0;
    if(z!=z_glob){
	setnew_z(z);
    }
    if(!bool_init_xiNL_barbar){
        init_xiNL_barbar(z);
    }
    double logr=log10(r);
    result=gsl_spline_eval (xiNL_barbar_spline,logr, xiNL_barbar_acc);
    
    return pow(10.,result)-1.;
}

void cosmology::init_xiNL_barbar(double z)
{
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for numerical calculation of xiNLbarbar was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance
    double rrmin=-3.0;
    double rrmax=3.0;

    double xx[Nxibar],yy[Nxibar];
    for (int i=0;i<Nxibar;i++)
    {
        double r=rrmin+i/(Nxibar-1.)*(rrmax-rrmin);
        xx[i]=r;
        r=pow(10.,r);
        yy[i]=log10(1.+xiNL_barbar(r,z));
    }
    xiNL_barbar_acc = gsl_interp_accel_alloc ();
    xiNL_barbar_spline = gsl_spline_alloc (gsl_interp_cspline, Nxibar);
    gsl_spline_init (xiNL_barbar_spline,xx,yy,Nxibar);

    bool_init_xiNL_barbar=true;
    if(verbose){
    std::cout<<"# "<<"The spline for numerical calculation of xiNLbarbar is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
}

double dxinlbarbar(double x, void* params)
{
    ///Int_{0}^{1.0} dy y^4 \xi_gg( y*r )
    xi_params c1 = *(xi_params *)params;

    cosmology *c2;
    double *r;
    double *z;

    c2=c1.cptr;
    r=c1.r;
    z=c1.z;

    double arg=x*(*r);
    double res=(*c2).xi_NL_num(arg,*z);
    res=res*pow(x,4.);

    return res;
}

double cosmology::xiNL_barbar(double r,double z){
	/// QAGW integral
	///Int_{0}^{1.0} dy y^2 \xi_gg( y*r )
	///Int_{-15.0}^{0.0} d ln y y^3 \xi_gg( y*r )

	double result, error;
	gsl_function F;
	F.function = &(dxinlbarbar);

	xi_params p;
	p.cptr = this;
	p.r = &(r);
	p.z = &z;

	F.params = &p;

	int stat;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	//stat=gsl_integration_qag (&F, -15.0, 0.0, 0, 1.e-4, 1000, 5, w, &result, &error);
	stat=gsl_integration_qags (&F, 0.0, 1.0, 0, 1.e-3, 1000, w, &result, &error);
	gsl_integration_workspace_free (w);

        return result*5.;
}

double cosmology::xi_NL_kaiser(double r, double z, double mu, double fkai){
    double xinl=xi_NL_num(r,z);
    double xinlbar=xiNL_bar_num(r,z);
    double xinlbarbar=xiNL_barbar_num(r,z);
    double musq=mu*mu;
    double p2mu=(3.*musq-1.)/2.;
    double p4mu=(35.*musq*musq-30.*musq+3.)/8.;

    double res=(xinl*(1.+2./3.*fkai+1./5.*fkai*fkai)+p2mu*(4./3.*fkai+4./7.*fkai*fkai)*(xinl-xinlbar)+p4mu*(8./35*fkai*fkai)*(xinl+5./2.*xinlbar-7./2.*xinlbarbar));
    //res=res*bbiassq;
    return res;
}

/*
double dwpnl_kaiser(double x, void* params)
{
    /// Wp(R) = 2.* R Int_{0}^{xmax} xi( \sqrt(x^2+1) R ) dx
    /// R is multiplied outside the integral
    projwpk_params c1 = *(projwpk_params *)params;

    cosmology *c2;
    double *R;
    double *z;
    double fkai;

    c2=c1.cptr;
    R=c1.R;
    z=c1.z;
    fkai=*(c1.fkai);

    double arg=sqrt(x*x+1.)*(*R);
    double mu=x/sqrt(1.+x*x);
    double res=(*c2).xi_NL_kaiser(arg,*z,mu,fkai);

    return res;
}

double cosmology::wpnl_kaiser(double z, double rad, double projmax,double fkai){

	/// QAGW integral
	///Int_{0}^{ymax} dy 2  R \xi_gg( \sqrt(y^2+1) R)
	/// ymax=pimax/R
	double result, error;
	double ymax=(projmax/rad);

	gsl_function F;
	F.function = &(dwpnl_kaiser);

	projwpk_params p;
	p.cptr = this;
	p.R = &(rad);
	p.z = &z;
        p.fkai=&fkai;

	F.params = &p;

	int stat;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	stat=gsl_integration_qag (&F, 0.0, ymax, 0, 1.e-3, 1000, 5, w, &result, &error);
	gsl_integration_workspace_free (w);

        result=2.0*result*rad;
        return result;
}
*/

double dwpnl_kaiser(double x, void* params)
{
    /// Wp(R) = 2.*Int_{0}^{xmax} xi(rp/sqrt(1-mu*mu)) d mu/(1-mu*mu)^{3/2}
    /// R is multiplied outside the integral
    projwpk_params c1 = *(projwpk_params *)params;

    cosmology *c2;
    double *R;
    double *z;
    double fkai;

    c2=c1.cptr;
    R=c1.R;
    z=c1.z;
    fkai=*(c1.fkai);

    double sinm=sqrt(1-x*x);
    double arg=(*R)/sinm;
    double res=(*c2).xi_NL_kaiser(arg,*z,x,fkai);

    return res/pow(sinm,3.0);
}

double cosmology::wpnl_kaiser(double z, double rad, double projmax,double fkai){

	/// QAGW integral
	///Int_{0}^{ymax} dy 2  R \xi_gg( \sqrt(y^2+1) R)
	/// ymax=pimax/R
	double result, error;
	double ymax=(1./sqrt(rad/projmax*rad/projmax+1.0));

	gsl_function F;
	F.function = &(dwpnl_kaiser);

	projwpk_params p;
	p.cptr = this;
	p.R = &(rad);
	p.z = &z;
        p.fkai=&fkai;

	F.params = &p;

	int stat;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	stat=gsl_integration_qag (&F, 0.0, ymax, 0, 1.e-3, 1000, 5, w, &result, &error);
	gsl_integration_workspace_free (w);

        result=2.0*rad*result;
        return result;
}

/*
void cosmology::setbias(double b){
    bbiassq=b*b;
}
*/

/// ============================ Linear part ======================
double cosmology::xiL_bar_num(double r,double z){
    double result=0.0;
    if(z!=z_glob){
	setnew_z(z);
    }
    if(!bool_init_xiL_bar){
        init_xiL_bar(z);
    }
    double logr=log10(r);
    result=gsl_spline_eval (xiL_bar_spline,logr, xiL_bar_acc);

    return pow(10.,result)-1.;
}

void cosmology::init_xiL_bar(double z)
{
    gsl_set_error_handler_off();
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for numerical calculation of xiLbar was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance
    double rrmin=-3.0;
    double rrmax=3.0;

    double xx[Nxibar],yy[Nxibar];
    for (int i=0;i<Nxibar;i++)
    {
        double r=rrmin+i/(Nxibar-1.)*(rrmax-rrmin);
        xx[i]=r;
        r=pow(10.,r);
        yy[i]=log10(1.+xiL_bar(r,z));
        //printf("%lg %lg\n",xx[i],yy[i]);
    }
    xiL_bar_acc = gsl_interp_accel_alloc ();
    xiL_bar_spline = gsl_spline_alloc (gsl_interp_cspline, Nxibar);
    gsl_spline_init (xiL_bar_spline,xx,yy,Nxibar);

    bool_init_xiL_bar=true;
    if(verbose){
    std::cout<<"# "<<"The spline for numerical calculation of xiLbar is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
}

double dxilbar(double x, void* params)
{
    ///Int_{0}^{1.0} dy y^2 \xi_gg( y*r )
    xi_params c1 = *(xi_params *)params;

    cosmology *c2;
    double *r;
    double *z;

    c2=c1.cptr;
    r=c1.r;
    z=c1.z;

    double arg=x*(*r);
    double res=(*c2).xi_L_num(arg,*z);
    res=res*x*x;

    return res;
}

double cosmology::xiL_bar(double r,double z){
	/// QAGW integral
	///Int_{0}^{1.0} dy y^2 \xi_gg( y*r )
	///Int_{-15.0}^{0.0} d ln y y^3 \xi_gg( y*r )

	double result, error;
	gsl_function F;
	F.function = &(dxilbar);

	xi_params p;
	p.cptr = this;
	p.r = &(r);
	p.z = &z;

	F.params = &p;

	int stat;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	//stat=gsl_integration_qag (&F, -15.0, 0.0, 0, 1.e-4, 1000, 5, w, &result, &error);
	stat=gsl_integration_qags (&F, 0., 1.0, 0, 1.e-3, 1000, w, &result, &error);
	gsl_integration_workspace_free (w);

        return result*3.;
}

double cosmology::xiL_barbar_num(double r,double z){
    double result=0.0;
    if(z!=z_glob){
	setnew_z(z);
    }
    if(!bool_init_xiL_barbar){
        init_xiL_barbar(z);
    }
    double logr=log10(r);
    result=gsl_spline_eval (xiL_barbar_spline,logr, xiL_barbar_acc);
    
    return pow(10.,result)-1.;
}

void cosmology::init_xiL_barbar(double z)
{
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for numerical calculation of xiLbarbar was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance
    double rrmin=-3.0;
    double rrmax=3.0;

    double xx[Nxibar],yy[Nxibar];
    for (int i=0;i<Nxibar;i++)
    {
        double r=rrmin+i/(Nxibar-1.)*(rrmax-rrmin);
        xx[i]=r;
        r=pow(10.,r);
        yy[i]=log10(1.+xiL_barbar(r,z));
    }
    xiL_barbar_acc = gsl_interp_accel_alloc ();
    xiL_barbar_spline = gsl_spline_alloc (gsl_interp_cspline, Nxibar);
    gsl_spline_init (xiL_barbar_spline,xx,yy,Nxibar);

    bool_init_xiL_barbar=true;
    if(verbose){
    std::cout<<"# "<<"The spline for numerical calculation of xiLbarbar is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
}

double dxilbarbar(double x, void* params)
{
    ///Int_{0}^{1.0} dy y^4 \xi_gg( y*r )
    xi_params c1 = *(xi_params *)params;

    cosmology *c2;
    double *r;
    double *z;

    c2=c1.cptr;
    r=c1.r;
    z=c1.z;

    double arg=x*(*r);
    double res=(*c2).xi_L_num(arg,*z);
    res=res*pow(x,4.);

    return res;
}

double cosmology::xiL_barbar(double r,double z){
	/// QAGW integral
	///Int_{0}^{1.0} dy y^4 \xi_gg( y*r )

	double result, error;
	gsl_function F;
	F.function = &(dxilbarbar);

	xi_params p;
	p.cptr = this;
	p.r = &(r);
	p.z = &z;

	F.params = &p;

	int stat;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	//stat=gsl_integration_qag (&F, -15.0, 0.0, 0, 1.e-4, 1000, 5, w, &result, &error);
	stat=gsl_integration_qags (&F, 0., 1.0, 0, 1.e-3, 1000, w, &result, &error);
	gsl_integration_workspace_free (w);

        return result*5.;
}

double cosmology::xi_L_kaiser(double r, double z, double mu, double fkai){
    double xil=xi_L_num(r,z);
    double xilbar=xiL_bar_num(r,z);
    double xilbarbar=xiL_barbar_num(r,z);
    double musq=mu*mu;
    double p2mu=(3.*musq-1.)/2.;
    double p4mu=(35.*musq*musq-30.*musq+3.)/8.;

    /*
    double bb0=1.0 + (2./3.) * fkai + (1./5.) * fkai*fkai;
    double bb2= (4./3.) * fkai + (4./7.) * fkai*fkai;
    double bb4= (8./35.) * fkai*fkai;

    double fac1=bb0+bb2*p2mu+bb4*p4mu;
    double fac2=-bb2*p2mu+5./2.*bb4*p4mu;
    double fac3=-7./2.*bb4*p4mu;
    */

    double res=(xil*(1.+2./3.*fkai+1./5.*fkai*fkai)+p2mu*(4./3.*fkai+4./7.*fkai*fkai)*(xil-xilbar)+p4mu*(8./35.*fkai*fkai)*(xil+5./2.*xilbar-7./2.*xilbarbar));

    //double res=(fac1*xil+fac2*xilbar+fac3*xilbarbar);
    //fprintf(stdout," %e %e %e ",xil,xilbar,xilbarbar);
    //res=res*bbiassq;
    return res;
}

/*
double dwpl_kaiser(double x, void* params)
{
    /// Wp(R) = 2.* R Int_{0}^{xmax} xi( \sqrt(x^2+1) R ) dx
    /// R is multiplied outside the integral
    projwpk_params c1 = *(projwpk_params *)params;

    cosmology *c2;
    double *R;
    double *z;
    double fkai;

    c2=c1.cptr;
    R=c1.R;
    z=c1.z;
    fkai=*(c1.fkai);

    double arg=sqrt(x*x+1.)*(*R);
    double mu=x/sqrt(1.+x*x);
    double res=(*c2).xi_L_kaiser(arg,*z,mu,fkai);

    return res;
}

double cosmology::wpl_kaiser(double z, double rad, double projmax,double fkai){

	/// QAGW integral
	///Int_{0}^{ymax} dy 2  R \xi_gg( \sqrt(y^2+1) R)
	/// ymax=pimax/R
	double result, error;
	double ymax=(projmax/rad);

	gsl_function F;
	F.function = &(dwpl_kaiser);

	projwpk_params p;
	p.cptr = this;
	p.R = &(rad);
	p.z = &z;
        p.fkai=&fkai;

	F.params = &p;

	int stat;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	stat=gsl_integration_qag (&F, 0.0, ymax, 0, 1.e-3, 1000, 5, w, &result, &error);
	gsl_integration_workspace_free (w);

        result=2.0*result*rad;
        return result;
}
*/

double dwpl_kaiser(double x, void* params)
{
    /// Wp(R) = 2.*Int_{0}^{xmax} xi(rp/sqrt(1-mu*mu)) d mu/(1-mu*mu)^{3/2}
    /// R is multiplied outside the integral
    projwpk_params c1 = *(projwpk_params *)params;

    cosmology *c2;
    double *R;
    double *z;
    double fkai;

    c2=c1.cptr;
    R=c1.R;
    z=c1.z;
    fkai=*(c1.fkai);

    double sinm=sqrt(1-x*x);
    double arg=(*R)/sinm;
    //fprintf(stdout,"%e %e %e ",*R,x,arg);
    double res=(*c2).xi_L_kaiser(arg,*z,x,fkai);

    res=res/pow(sinm,3.0);
    //fprintf(stdout," %e \n",res);
    return res;
}

double cosmology::wpl_kaiser(double z, double rad, double projmax,double fkai){

	/// QAGW integral
	///Int_{0}^{ymax} dy 2  R \xi_gg( \sqrt(y^2+1) R)
	/// ymax=pimax/R
	double result, error;
	double ymax=(1./sqrt(rad/projmax*rad/projmax+1.0));

	gsl_function F;
	F.function = &(dwpl_kaiser);

	projwpk_params p;
	p.cptr = this;
	p.R = &(rad);
	p.z = &z;
        p.fkai=&fkai;

	F.params = &p;

	int stat;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	stat=gsl_integration_qag (&F, 0.0, ymax, 0, 1.e-3, 1000, 5, w, &result, &error);
	gsl_integration_workspace_free (w);

        result=2.0*rad*result;
        //exit(1111);
        return result;
}

