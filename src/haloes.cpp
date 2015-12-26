/// Power spectrum routines part of cosmology class
#include "cosmology.h"
#include <gsl/gsl_sf_trig.h>


/// Mass function wrapper
double cosmology::nofm(double M, double z)
{
    double result;
    if(opt_mf==1)
    {
        result=MF_TI10(M,z);
    }else if(opt_mf==2)
    {
	result=MF_ST(M,z);
    }else if(opt_mf==3)
    {
	result=MF_BH(M,z);
    }else
    {
        std::cout<<"Mass function option not supported yet"<<std::endl;
        exit(0);
    }
    return result;
}

/// Tinker et al. 2009 SO (350) mass function The mass function is in comoving coordinates.
/// n(M) dM has the units of h^3 Mpc^{-3}
double cosmology::MF_TI09_350(double M, double z)
{
    double sig=sqrt(varM_TH_num(M,z));
    double dlogsigdlogm= 0.5*(log(varM_TH_num(0.99*M,z))-log(varM_TH_num(1.01*M,z)))/(log(0.99*M)-log(1.01*M));
    double A0=0.206;
    double sa0=1.54,b0=2.15,c0=1.31;
    double Aofz=A0*pow(1.+z,-0.14);
    double saofz=sa0*pow(1.+z,-0.06);
    double alp350=pow(10.,-pow(0.75/log10(350.0/75.0),1.2));
    double bofz=b0*pow(1.+z,-alp350);
    double cofz=c0;
    double fsig=Aofz*(pow(sig/bofz,-saofz)+1.)*exp(-cofz/pow(sig,2.));

    return rho_crit_0*Omega0*fsig/pow(M,2.)*fabs(dlogsigdlogm);
}

/// Tinker et al. 2010 SO (200) mass function The mass function is in comoving coordinates.
/// n(M) dM has the units of h^3 Mpc^{-3}
double cosmology::MF_TI10(double M, double z)
{
    if(z!=z_glob){
	setnew_z(z);
    }
    if(!init_Tink){
	initTinker(z);
    }

    double dc=1.686;
    double sig=sqrt(varM_TH_num(M,z));
    double xnu=dc/sig;

    // Note that the variance can be calced either at z=0.0 or z=z.
    double dlogsigdlogm= 0.5*(log(varM_TH_num(0.99*M,z))-log(varM_TH_num(1.01*M,z)))/(log(0.99*M)-log(1.01*M));

    double beta= 0.589*pow(1.+z, 0.20);
    double phi =-0.729*pow(1.+z,-0.08);
    double eta =-0.243*pow(1.+z, 0.27);
    double gama= 0.864*pow(1.+z,-0.01);

    double result=alpTink*(1.+pow(beta*xnu,-2.*phi))*pow(xnu,2.0*eta);
    result=result*exp(-gama*pow(xnu,2.0)/2.0);
    result=result*xnu;
    //printf("#DEBUG: %e %e %e %e ",M,result,fabs(dlogsigdlogm),xnu);
    result=rho_crit_0*Omega0*result/pow(M,2.)*fabs(dlogsigdlogm);

    //std::cout<<M<<" "<<result<<" "<<sig<<" "<<dlogsigdlogm<<" "<<alpTink/0.368<<std::endl;
    
    return result;

}

/// Tinker et al. 2010 mass function normalization. This is defines so that:
//  \int_{0}^{\infty} f(\nu) d\nu = 1
void cosmology::initTinker(double z)
{
    
    double beta= 0.589*pow(1.+z, 0.20);
    double phi =-0.729*pow(1.+z,-0.08);
    double eta =-0.243*pow(1.+z, 0.27);
    double gama= 0.864*pow(1.+z,-0.01);

    double fac1=pow(2.0/gama,eta+0.5);
    double fac2=pow(2.0/gama,eta-phi+0.5);
    double fac3=pow(beta,-2.0*phi);
    double fac4=exp(gsl_sf_lngamma(eta+0.5));
    double fac5=exp(gsl_sf_lngamma(eta-phi+0.5));

    alpTink=(2.0/(fac1*fac4+fac3*fac2*fac5));

    init_Tink=true;
}

/// Tinker et al. 2009 SO (200) mass function The mass function is in comoving coordinates.
/// n(M) dM has the units of h^3 Mpc^{-3}
double cosmology::MF_TI09(double M, double z)
{
    double sig=sqrt(varM_TH_num(M,z));
    double dlogsigdlogm= 0.5*(log(varM_TH_num(0.99*M,z))-log(varM_TH_num(1.01*M,z)))/(log(0.99*M)-log(1.01*M));
    double A0=0.186;
    double sa0=1.47,b0=2.57,c0=1.19;
    double Aofz=A0*pow(1.+z,-0.14);
    double saofz=sa0*pow(1.+z,-0.06);
    double alp200=pow(10.,-pow(0.75/log10(200.0/75.0),1.2));
    double bofz=b0*pow(1.+z,-alp200);
    double cofz=c0;
    double fsig=Aofz*(pow(sig/bofz,-saofz)+1.)*exp(-cofz/pow(sig,2.));

    return rho_crit_0*Omega0*fsig/pow(M,2.)*fabs(dlogsigdlogm);
}

/// Sheth Tormen mass function The mass function is in comoving coordinates.
/// n(M) dM has the units of h^3 Mpc^{-3}
double cosmology::MF_ST(double M, double z)
{
    double sig=sqrt(varM_TH_num(M,z)); //_num(M,z);
    double dlogsigdlogm= 0.5*(log(varM_TH_num(0.99*M,z))-log(varM_TH_num(1.01*M,z)))/(log(0.99*M)-log(1.01*M));

    double nu=1.686/sig;
    /// One may want to use Delta_c which depends weakly on Omega(z) here.

    double nu2=0.84083292 * nu;

    double dodl=0.3222*1.0/(sqrt(2.0*M_PI))*nu2*exp(-0.5*pow(nu2,2.0))*(1.+1.0/pow(nu2,0.6));

    return 2.0*rho_crit_0*Omega0/pow(M,2.)*fabs(dlogsigdlogm)*dodl;
}

/// Bhattacharya et al. mass function The mass function is in comoving coordinates.
/// n(M) dM has the units of h^3 Mpc^{-3}
double cosmology::MF_BH(double M, double z)
{
    double sig=sqrt(varM_TH_num(M,z));
    double dlogsigdlogm= 0.5*(log(varM_TH_num(0.99*M,z))-log(varM_TH_num(1.01*M,z)))/(log(0.99*M)-log(1.01*M));
    double Atilde=0.333/pow(1+z,0.11);
    double atilde=0.788/pow(1+z,0.01);
    double ptilde=0.807;
    double qtilde=1.795;
    double nu=1.686/sig;
    double fsig=Atilde*sqrt(2/M_PI)*exp(-atilde/2.*nu*nu)*(1+pow(nu*nu*atilde,-ptilde))*pow(nu*sqrt(atilde),qtilde);

    /// Eq. 5 from Warren et al. 2006
    return rho_crit_0*Omega0*fsig/pow(M,2.)*fabs(dlogsigdlogm);
}

/// Warren et al. 2006 mass function The mass function is in comoving coordinates.
/// n(M) dM has the units of h^3 Mpc^{-3}
double cosmology::MF_WA(double M, double z)
{
    double sig=sqrt(varM_TH_num(M,z)); //_num(M,z);
    double dlogsigdlogm= 0.5*(log(varM_TH_num(0.99*M,z))-log(varM_TH_num(1.01*M,z)))/(log(0.99*M)-log(1.01*M));
    double fsig=0.7234*(pow(sig,-1.625)+0.2538)*exp(-1.1982/pow(sig,2.));

    /// Eq. 5 from Warren et al. 2006
    return rho_crit_0*Omega0*fsig/pow(M,2.)*fabs(dlogsigdlogm);
}

// Bias wrapper
double cosmology::bias(double M, double z)
{
    double result;
    if(opt_b==1)
    {
        result=bias_TI10(M,z);
    }else
    {
        std::cout<<"Bias function option not supported yet"<<std::endl;
        exit(0);
    }
    return result;

}

/// Tinker et al. 2010 bias for SO(200)
/// This should be used alongwith MF_TI10
double cosmology::bias_TI10(double M, double z)
{
    double sig=sqrt(varM_TH_num(M,0.0));
    double dc=1.686;
    double dcz=dc;
    if(z>0.0)
    {
        dcz=dc/growthfactor_num(z);
    }
    double xnu=dcz/sig;

    double Delta=200.0;
    double y=log10(Delta);

    double expfac=exp(-pow(4./y,4.0));
    double A=1.0+0.24*y*expfac;
    double a=0.44*y-0.88;
    double B=0.183;
    double b=1.5;
    double C=0.019+0.107*y+0.19*expfac;
    double c=2.4;

    double xnua=pow(xnu,a);
    double xnub=pow(xnu,b);
    double xnuc=pow(xnu,c);

    return 1.0 - A*xnua/(xnua+pow(dc,a)) + B*xnub + C*xnuc;
}

/// Tinker et al. 2006 bias 
/// There is a small Omega(z) dependence on dc which has been ignored here.
double cosmology::bias_TWZZ(double M, double z)
{
    double sig=sqrt(varM_TH_num(M,z));
    double dc=1.68647;
    /* Was a bug here. dc should not be changed.
    if(z>0.0)
    {
        dc=dc/growthfactor_num(z);
    }*/
    double xnu=dc/sig;

    double aa=0.707,bb=0.35,cc=0.80;

    double f1= aa*xnu*xnu;

    double fac1=sqrt(aa)*(f1+bb*pow(f1,1.-cc));
    double fac2=pow(f1,cc)+bb*(1.-cc)*(1.-cc/2.);
    fac2=pow(f1,cc)/fac2; //Potential for improvement
    return 1.0 + (fac1-fac2)/(sqrt(aa)*dc);
}

/// Number of haloes with mass greater than or equal to M
double cosmology::Nplus(double M, double z)
{
    //Int_{10^9}^{10^{16}} n(M) dM
    double result=0;
    const int N=100;
    double x[N],w[N];

    gauleg(log10(M),17.0,x,w,N);
    result=0;

    for (int i=0;i<N;i++)
    {
        double m=pow(10.,x[i]);
        result=result+w[i]*nofm(m,z)*m*log(10.);
    }
    return result;

}

/// Root of this function yields the mass which has a given N+
double findM(double xm, void *params)
{
    np_params c1 = *(np_params *) params;
//    std::cout<<1<<std::endl;
    cosmology *c2;
    double *z;
    double *Np;
    c2=c1.cptr;
    z=c1.z;
    Np=c1.Np;
    double m=pow(10.0,xm);
    double Npl=(*c2).Nplus(m,*z);

    return Npl-(*Np);

}

/// Get the M which has N_{+} number of haloes above mass M.
double cosmology::getM(double Np,double z)
{
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    // m is in log_{10}(m)
    double xm_lo = 9.0, xm_hi = 15.5;

    gsl_function F;
    np_params p;
    p.cptr = this;
    p.z = &z;
    p.Np = &Np;

    F.function = &findM;
    F.params = &p;
   
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, xm_lo, xm_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        xm_lo = gsl_root_fsolver_x_lower (s);
        xm_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (xm_lo, xm_hi,0, 1e-6);

        if (status == GSL_SUCCESS)
        {
            std::cout<<"# "<<"getM:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return pow(10.,res);

}


/// The root of this function gives the collapse redshift
double findzcoll(double z, void *params)
{
    coll_params c1 = *(coll_params *) params;
    cosmology *c2;
    double *sig;
    c2=c1.cptr;
    sig=c1.sig;

    return (*sig) - 1.686/(*c2).growthfactor_num(z);

}

/// Get the collapse redshift for a given variance
double cosmology::getzcoll(double sigma)
{
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double z_lo = 0.0, z_hi = 10.0;

    gsl_function F;
    coll_params p;
    p.cptr = this;
    p.sig = &sigma;

    F.function = &findzcoll;
    F.params = &p;
   
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
            //std::cout<<"# "<<"zcollapse:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return res;

}

/// Concentration parameter wrapper
double cosmology::conc(double M, double z)
{
    double result;

    if(opt_c==1)
    {
        result=c_MAC(M, z);
    } else
    {
        std::cout<<"Concentration parameter option not supported yet."<<std::endl;
        exit(0);
    }

    return result*cfactor;
}

/// The concentration parameter ala Maccio et al. 2007. Note that the mass is Mvir in h^{-1} Msun
double cosmology::c_MAC(double M, double z)
{
    double concen=0;
    double F=0.001;
    double K=2.6;

    double sigma=sqrt(varM_TH_num(F*M,0.0));

    if(sigma<1.686/growthfactor_num(z))
    {
        concen=K;
    }else
    {
        concen=K*((1.0 + getzcoll(sigma)) / (1.0 + z));
    }

    return concen;
}

/// The root of this function gives the Mvir for a particular M200
double findmvir(double mvir, void *params)
{
    mvir_params c1 = *(mvir_params *) params;
    cosmology *c2;
    double *m200;
    double *z;
    c2=c1.cptr;
    m200=c1.m200;
    z=c1.z;

    double cvir=(*c2).conc(mvir,*z);
    double c200=(*c2).getc200(cvir,*z);

    double fcvir=log(1.+cvir) - cvir/(1.+cvir);
    double fc200=log(1.+c200) - c200/(1.+c200);

    return *m200/mvir-fc200/fcvir;

}

/// Get the Mvir for a given M200
double cosmology::getMvir(double M, double z)
{
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double m_lo = 0.1*M, m_hi = 5.0*M;

    gsl_function F;
    mvir_params p;
    p.cptr = this;
    p.m200 = &M;
    p.z=&z;

    F.function = &findmvir;
    F.params = &p;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, m_lo, m_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        m_lo = gsl_root_fsolver_x_lower (s);
        m_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (m_lo, m_hi,0, 1.e-6);

        if (status == GSL_SUCCESS)
        {
            //std::cout<<"# "<<"zcollapse:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return res;

}

/// The root of this function gives the c200 for a particular cvir
double findc200(double c200, void *params)
{
    c200_params c1 = *(c200_params *) params;
//    std::cout<<1<<std::endl;
    cosmology *c2;
    double *cvir;
    double *omegaz;
    double *dcz;
    c2=c1.cptr;
    cvir=c1.cvir;
    omegaz=c1.omegaz;
    dcz=c1.dcz;

    double fcvir=log(1.+*cvir) - *cvir/(1.+*cvir);
    double fc200=log(1.+c200) - c200/(1.+c200);

    return (*dcz)/(200.0*(*omegaz))*pow((*cvir)/c200,3.0) - fcvir/fc200;

}

/// Get the c200 for a given cvir
double cosmology::getc200(double cvir, double z)
{
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double c_lo = 0.01*cvir, c_hi = 10.0*cvir;

    gsl_function F;
    c200_params p;
    p.cptr = this;
    p.cvir = &cvir;
    double omz=Omega(z);
    double dcz=Delta_crit(z);
    p.omegaz=&omz; 
    p.dcz=&dcz;

    F.function = &findc200;
    F.params = &p;
   
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, c_lo, c_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        c_lo = gsl_root_fsolver_x_lower (s);
        c_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (c_lo, c_hi,0, 1e-6);

        if (status == GSL_SUCCESS)
        {
            //std::cout<<"# "<<"zcollapse:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return res;

}

/// The root of this function gives Mstar such that sigma(Mstar)=1.68
double findmstar(double m, void *params)
{
    c_params c1 = *(c_params *) params;
    cosmology *c2;
    c2=c1.cptr;

    return 1.68 - sqrt((*c2).varM_TH_num(pow(10.,m),0.0));

}

/// Get Mstar such that sigma(Mstar)=1.68
double cosmology::getmstar()
{
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double m_lo = 10.0, m_hi = 15.0;

    gsl_function F;
    c_params p;
    p.cptr = this;

    F.function = &findmstar;
    F.params = &p;
   
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, m_lo, m_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        m_lo = gsl_root_fsolver_x_lower (s);
        m_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (m_lo, m_hi,0, 1e-6);

        if (status == GSL_SUCCESS)
        {
            //std::cout<<"# "<<"zcollapse:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return res;

}

/*
/// Radii are in comoving units here
void cosmology::modelNFWhalo_com(double m200,double z,double &Mvir, double &Rvir, double &cvir)
{
    Mvir=getMvir(m200,z);

    Rvir=0.169*pow(Mvir/1.e12,1./3.);
    Rvir*=pow(Delta_crit(z)/178.0,-1.0/3.0);
    Rvir*=pow(Eofz(z),-2./3.);

    cvir=conc(Mvir,z);

    Rvir=Rvir*(1.+z);

}
*/

/// Get comoving virial radius from virial mass
double cosmology::getRvirfromMvir(double Mvir, double z){

    double Rvir=0.169*pow(Mvir/1.e12,1./3.);
    Rvir*=pow(Delta_crit(z)/178.0,-1.0/3.0);
    Rvir*=pow(Eofz(z),-2./3.);
    Rvir=Rvir*(1.+z);

    return Rvir;
}

/// Radii are in comoving units here
void cosmology::modelNFWhalo_com(double m200,double z,double &Mvir, double &Rvir, double &cvir, double &R200, double &c200)
{
    Mvir=getMvir(m200,z);
    Rvir=0.169*pow(Mvir/1.e12,1./3.);
    Rvir*=pow(Delta_crit(z)/178.0,-1.0/3.0);
    Rvir*=pow(Eofz(z),-2./3.);

    double Delta=200.;
    R200=pow( 1.0e12/(4./3.*Delta*3e4*0.3/(8*gee)),1./3.);
    R200*=pow(m200/1.e12,1./3.);
    R200*=pow(Omega(z)/0.3,-1.0/3.0);
    R200*=pow(Eofz(z),-2./3.);

    cvir=conc(Mvir,z);

    c200=getc200(cvir,z);
    R200=R200*(1.+z);
    Rvir=Rvir*(1.+z);

}

/// Radii are in physical units here, so be careful
void cosmology::modelNFWhalo(double m200,double z,double &Mvir, double &Rvir, double &cvir, double &R200, double &c200)
{
    Mvir=getMvir(m200,z);

    Rvir=0.169*pow(Mvir/1.e12,1./3.);
    Rvir*=pow(Delta_crit(z)/178.0,-1.0/3.0);
    Rvir*=pow(Eofz(z),-2./3.);

    double Delta=200.;
    R200=pow( 1.0e12/(4./3.*Delta*3e4*0.3/(8*gee)),1./3.);
    R200*=pow(m200/1.e12,1./3.);
    R200*=pow(Omega(z)/0.3,-1.0/3.0);
    R200*=pow(Eofz(z),-2./3.);

    cvir=conc(Mvir,z);

    c200=getc200(cvir,z);

}

/// The fourier transform of the NFW profile on a grid of k*rs and c
/// krsmin=10^-6 and krsmax=10.0^8.0
/// log(cmin)=0.0 and logcmax=2.5
/// The arrays are calced from a Fortran numerical recipes code
void cosmology::ukinit2()
{
    double krsdiff=(krsmax-krsmin)/(Nuk-1.);
    double cdiff=(cmax-cmin)/(Nuk-1.);

    FILE *fp;
    fp=fopen("DATA/utildearrays","r");

    for(int i=0;i<Nuk;i++)
    {
	uk_krs[i]=krsmin+i*krsdiff;
	uk_c  [i]=cmin  +i*cdiff;
    }

    /// Bilinear interpolation is not defined in GSL. So code this up yourself
    //std::cout<<"# Start check here"<<std::endl;
    for(int i=0;i<Nuk;i++)
    {
	for(int j=0;j<Nuk;j++)
	{
	    double jk1,jk2;
	    if(fscanf(fp,"%lg %lg %lg",&jk1,&jk2,&(ukrsc[i][j])) != 3){
		std::cout<<"Error reading DATA/utildearrays. Aborting\n"<<std::endl;
		exit(0);
	    }
	    std::cout<<std::scientific;
	    //std::cout<<uk_krs[j]<<" "<<uk_c[i]<<" "<<ukrsc[j][i]<<std::endl;
	    //std::cout<<pow(10.,uk_krs[j])<<" "<<pow(10.,uk_c[i])<<" "<<pow(10.,ukrsc[j][i])<<std::endl;
	}
	//std::cout<<std::endl;
    }

    fclose(fp);
    
    uk_c_acc=gsl_interp_accel_alloc();
    uk_krs_acc=gsl_interp_accel_alloc();
    bool_inituk=true;
}

/// The fourier transform of the NFW profile on a grid of k*rs and c
/// krsmin=10^-6 and krsmax=10.0^8.0
/// log(cmin)=0.0 and logcmax=2.5
void cosmology::ukinit()
{
    double krsdiff=(krsmax-krsmin)/(Nuk-1.);
    double cdiff=(cmax-cmin)/(Nuk-1.);

    for(int i=0;i<Nuk;i++)
    {
	uk_krs[i]=krsmin+i*krsdiff;
	uk_c  [i]=cmin  +i*cdiff;
    }

    /// Bilinear interpolation is not defined in GSL. So code this up yourself
    //std::cout<<"# Start check here"<<std::endl;
    for(int i=0;i<Nuk;i++)
    {
	double con=pow(10.,uk_c[i]);
	double fac=log(1.+con) - con/(1.+con);

	bool extrap=false;
	double x_extrap,y_extrap;
	for(int j=0;j<Nuk;j++)
	{
	    double krs=pow(10.,uk_krs[j]);

	    // Do the cisi calculation only when u(k)>10^-6
	    // A log-linear interpolation for values below
	    if(!extrap)
	    {
		double arg1=krs;
		double arg2=(1.0+con)*krs;
		double arg3=con*krs;

		double ci1,ci2;
		double si1,si2;

		ci1=gsl_sf_Ci(arg1);
		ci2=gsl_sf_Ci(arg2);
		si1=gsl_sf_Si(arg1);
		si2=gsl_sf_Si(arg2);
		
		double res1=sin(arg1)*(si2-si1);
		res1-=sin(arg3)/arg2;
		res1+=cos(arg1)*(ci2-ci1);

		if(res1/fac<=0.0){
		    std::cout<<"Some serious trouble here "<<res1<<" "<<fac<<" "<<krs<<" "<<con<<std::endl;
		    exit(0);
		}
		ukrsc[i][j]=log10(res1/fac);
	    //std::cout<<std::scientific;
	    //std::cout<<krs<<" "<<con<<" "<<res1/fac<<" "<<std::endl;
		if (ukrsc[i][j] < -6.0) {
		    extrap=true;
		    x_extrap=uk_krs[j];
		    y_extrap=ukrsc[i][j];
		}
	    }else{
		ukrsc[i][j]=y_extrap-1.99*(uk_krs[j]-x_extrap);
	    //std::cout<<std::scientific;
	    //std::cout<<krs<<" "<<con<<" "<<pow(10.,ukrsc[i][j])<<" "<<std::endl;
	    }
	}
	//std::cout<<std::endl;
    }
	//exit(0);
    
    uk_c_acc=gsl_interp_accel_alloc();
    uk_krs_acc=gsl_interp_accel_alloc();
    bool_inituk=true;
    
}

/// Bilinear interpolation for u given k*rs and c
double cosmology::ukinterp(double krs, double c)
{
    // Bilinear interpolation is not defined in GSL. So code this up yourself
    if(!bool_inituk)
    {
	// Init ukm
	ukinit();
    }
    double xc=log10(c);
    double xkrs=log10(krs);
    double result;

    int i,j;
    i=j=-99;
    // Complain for out of bound krs and c values
    if(xkrs>uk_krs[Nuk-1]){
        if(verbose){
	std::cout<<"krs out of bound in ukinterp, log10[krs]="<<xkrs<<std::endl;
	std::cout<<"Extrapolating to "<<xkrs<<std::endl;
        }
        j=Nuk-2;
	//exit(0);
    }
    if(xc>uk_c[Nuk-1]||xc<uk_c[0]){
        if(verbose){
	std::cout<<"c out of bound in ukinterp, log10[c]="<<xc<<std::endl;
	std::cout<<"Extrapolating to "<<xc<<std::endl;
        }
        i=(xc>uk_c[Nuk-1])?(Nuk-2):0;
	//exit(0);
    }
    if(xkrs<uk_krs[0]){
	// For small k u(k|M)=1.0
	result=1.0;
	return result;
    }

    // Look up con and krs index
    if(i==-99)
	i=gsl_interp_accel_find (uk_c_acc,uk_c,Nuk,xc);
    if(j==-99)
        j=gsl_interp_accel_find (uk_krs_acc,uk_krs,Nuk,xkrs);

    double y11=ukrsc[ i ][ j ];
    double y12=ukrsc[ i ][j+1];
    double y21=ukrsc[i+1][ j ];
    double y22=ukrsc[i+1][j+1];

    double kdis=(xkrs-uk_krs[j])/(uk_krs[j+1]-uk_krs[j]);
    double cdis=(xc-uk_c[i])/(uk_c[i+1]-uk_c[i]);

    double r1=y11+cdis*(y21-y11);
    double r2=y12+cdis*(y22-y12);

    result=r1+kdis*(r2-r1);
    return pow(10.,result);

}

/// Fourier transform of NFW profile
double cosmology::ukofm(double k, double m, double z)
{
    double mvir,rvir, cvir, r200, c200;
    modelNFWhalo(m,z,mvir,rvir,cvir,r200,c200);

    double rs=r200/c200;

    double fac=log(1.+c200) - c200/(1.+c200);

    double ci1,ci2;
    double si1,si2;

    double arg1=k*rs;
    double arg2=(1.+c200)*k*rs;
    double arg3=c200*k*rs;

    ci1=gsl_sf_Ci(arg1);
    ci2=gsl_sf_Ci(arg2);
    si1=gsl_sf_Si(arg1);
    si2=gsl_sf_Si(arg2);

    double res1=sin(arg1)*(si2-si1);
    res1-=sin(arg3)/arg2;
    res1+=cos(arg1)*(ci2-ci1);

    return res1/fac;

}

/// Fourier transform of NFW profile for satellites
double cosmology::uskofm(double k, double m, double z, double csbycdm)
{
    double mvir,rvir, cvir, r200, c200;
    modelNFWhalo(m,z,mvir,rvir,cvir,r200,c200);

    c200=c200*csbycdm;

    double rs=r200/c200;

    double fac=log(1.+c200) - c200/(1.+c200);

    double ci1,ci2;
    double si1,si2;

    double arg1=k*rs;
    double arg2=(1.+c200)*k*rs;
    double arg3=c200*k*rs;

    ci1=gsl_sf_Ci(arg1);
    ci2=gsl_sf_Ci(arg2);
    si1=gsl_sf_Si(arg1);
    si2=gsl_sf_Si(arg2);

    double res1=sin(arg1)*(si2-si1);
    res1-=sin(arg3)/arg2;
    res1+=cos(arg1)*(ci2-ci1);

    return res1/fac;

}

double cosmology::munfw(double x){
    return log(1.+x)-x/(1.+x);
}


/// The root of this function gives the concentration for a halo
/// evolving in isolation
double findcmarch(double cdelz, void *params)
{
    march_params c1 = *(march_params *) params;
    double fac;
    fac=c1.fac;

    double mucdelz = log(1.+cdelz)-cdelz/(1.+cdelz);
    double res= fac - pow(cdelz,3)/mucdelz;
    //std::cerr<<cdelz<<" "<<res<<std::endl;
    return res;

}

void cosmology::pevolve_fixed(double cdel, int opt, double z, double zstart, double
&cdelz, double &fdelz){

    double fac;
    if(opt==1){
        fac=pow(cdel*(1.+zstart)/(1.+z),3)/munfw(cdel);
    }else if(opt==2){
        fac=pow(cdel,3)/munfw(cdel);
        fac=fac*pow(Eofz(zstart)/Eofz(z),2);
        fac=fac*Delta_crit(zstart)/Delta_crit(z);
    }else if(opt==3){
        fac=pow(cdel,3)/munfw(cdel);
        fac=fac*pow(Eofz(zstart)/Eofz(z),2);
    }else{
        fprintf(stderr,"Option %d not supported yet, bailing out...",opt);
        exit(100);
    }

    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double c_lo = 0.01*cdel, c_hi = 100000.0*cdel;

    gsl_function F;
    march_params p;
    p.fac = fac;

    F.function = &findcmarch;
    F.params = &p;
   
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, c_lo, c_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        c_lo = gsl_root_fsolver_x_lower (s);
        c_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (c_lo, c_hi,0, 1e-6);

        if (status == GSL_SUCCESS)
        {
            //std::cout<<"# "<<"zcollapse:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    cdelz=res;
    fdelz=munfw(cdelz)/munfw(cdel);

}

/// The root of this function gives the cDel for a particular cvir
double findcDel(double cDel, void *params)
{
    cDel_params c1 = *(cDel_params *) params;
    cosmology *c2;
    double *cvir;
    double *omegaz;
    double *dcz;
    double *Delta;
    c2=c1.cptr;
    cvir=c1.cvir;
    omegaz=c1.omegaz;
    dcz=c1.dcz;
    Delta=c1.Delta;

    double fcvir=log(1.+*cvir) - *cvir/(1.+*cvir);
    double fcDel=log(1.+cDel) - cDel/(1.+cDel);

    double res= (*dcz)/(*Delta*(*omegaz))*pow((*cvir)/cDel,3.0) - fcvir/fcDel;
    //printf("%lg %lg %lg %lg %lg \n",fcvir,fcDel,*cvir,cDel,res);
    //fflush(stdout);
    return res;

}

/// Get the cDel for a given cvir
double cosmology::getcDel(double cvir, double z, double Delta)
{
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double c_lo = 0.01*cvir, c_hi = 10.0*cvir;

    gsl_function F;
    cDel_params p;
    p.cptr = this;
    p.cvir = &cvir;
    double omz=Omega(z);
    double dcz=Delta_crit(z);
    p.omegaz=&omz; 
    p.dcz=&dcz;
    p.Delta=&Delta;

    F.function = &findcDel;
    F.params = &p;
   
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, c_lo, c_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        c_lo = gsl_root_fsolver_x_lower (s);
        c_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (c_lo, c_hi,0, 1e-6);

        if (status == GSL_SUCCESS)
        {
            //std::cout<<"# "<<"zcollapse:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return res;

}

