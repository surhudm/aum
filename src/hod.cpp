#include "hod.h"
///Halo occupation distribution function
///Derived from the class cosmology.

int hod::initialize_const(){
    // Initialize some constants
    hod_rmin=-4.6;
    hod_rmax= 4.0;
    hod_rmax_u= 2.1;

    return 1;
}

///Constructor
hod::hod(cosmo p, hodpars h)
:cosmology(p)
{
    initialize_const();
    hodp=h; 
    if(verbose){
        std::cout<<"# HOD inited to: "<<"\n"
            <<"# Mmin "<<hodp.Mmin<<"\n" \
            <<"# siglogM "<<hodp.siglogM<<"\n" \
            <<"# Msat "<<hodp.Msat<<"\n" \
            <<"# alpsat "<<hodp.alpsat<<"\n" \
            <<"# Mcut "<<hodp.Mcut<<"\n" \
            <<"# fac "<<hodp.fac<<"\n" \
            <<"# csbycdm "<<hodp.csbycdm \
            <<std::endl;
    }
    halo_exc=true;
    whichopt=0;
    //std::cerr<<"HOD constructor initialized"<<std::endl;
    //exit(101);
    bool_init_D2gg=false;
    bool_init_D2gd=false;
    bool_init_xigg=false;
    bool_init_xigd=false;
    bool_init_xigg_bar=false;
    bool_init_xigg_barbar=false;
    miyatake21switch=false;
    fKaiser=-99.0;
    off_rbyrs=0.0;
    fcen_off=0.0;
    inc_alp=0.0;
    inc_xM=12.0;
#if TINK==2
    bool_init_nc=false;
    bool_init_ns=false;
#endif
}

hod::hod()
:cosmology()
{
    initialize_const();
    hodp.Mmin=13.0;
    hodp.siglogM=0.25;
    hodp.Msat=14.0;
    hodp.alpsat=1.0;
    hodp.Mcut=13.5;
    hodp.csbycdm=1.0;
    hodp.fac=1.0;
    if(verbose){
        std::cout<<"# HOD inited to: "<<"\n"
            <<"# Mmin "<<hodp.Mmin<<"\n" \
            <<"# siglogM "<<hodp.siglogM<<"\n" \
            <<"# Msat "<<hodp.Msat<<"\n" \
            <<"# alpsat "<<hodp.alpsat<<"\n" \
            <<"# Mcut "<<hodp.Mcut<<"\n" \
            <<"# fac "<<hodp.fac<<"\n" \
            <<"# csbycdm "<<hodp.csbycdm \
            <<std::endl;
    }
    halo_exc=true;
    whichopt=0;

    bool_init_D2gg=false;
    bool_init_D2gd=false;
    bool_init_xigg=false;
    bool_init_xigd=false;
    bool_init_xigg_bar=false;
    bool_init_xigg_barbar=false;
    miyatake21switch=false;
    fKaiser=-99.0;
    off_rbyrs=0.0;
    fcen_off=0.0;
    inc_alp=0.0;
    inc_xM=12.0;
#if TINK==2
    bool_init_nc=false;
    bool_init_ns=false;
#endif
}

void hod::sethalo_exc(bool mark){
    halo_exc=mark;
}

void hod::setMiyatake21_switch(bool mark){
    miyatake21switch=mark;
}

#if TINK==1
///And now all the different functions 
///Average number of central galaxies in a halo of mass m
double hod::ncen(double xm)
{
    double arg=(xm-hodp.Mmin)/hodp.siglogM;
    double res=0.5*(1+gsl_sf_erf(arg));
    return res;
}

///Average number of satellite galaxies in a halo of mass m
double hod::nsat(double xm)
{
    double arg=(xm-hodp.Mmin)/hodp.siglogM;
    double res=0.5*(1+gsl_sf_erf(arg));
    res=res*pow(10.,hodp.alpsat*(xm-hodp.Msat))*exp(-pow(10.,hodp.Mcut-xm));
    return res;
}
#elif TINK==2
///And now all the different functions 
///Average number of central galaxies in a halo of mass m
double hod::ncen(double xm)
{
    if(!bool_init_nc){
        std::cout<<"Spline for Nc not initialized\n";
    }
    double res=0.0;
    if(xm<nc_mmin || xm>nc_mmax){
        res=0.0;
    }else{
        res=pow(10.,gsl_spline_eval(nc_spline,xm,nc_acc));
    }
    if(res>1.0) res=1.0;
    return res;
}

///Average number of satellite galaxies in a halo of mass m
double hod::nsat(double xm)
{
    if(!bool_init_ns){
        std::cout<<"Spline for Ns not initialized\n";
    }
    double res=0.0;
    if(xm<ns_mmin || xm>ns_mmax){
        res=0.0;
    }else{
        res=pow(10.,gsl_spline_eval(ns_spline,xm,ns_acc));
    }
    return res;
}

void hod::init_Nc_spl(double xx[],double yy[],int Ncspl){
    if(!bool_init_nc){
        gsl_interp_accel_free(nc_acc);
        gsl_spline_free(nc_spline);
        bool_init_nc=false;
    }
    nc_acc=gsl_interp_accel_alloc();
    nc_spline=gsl_spline_alloc(gsl_interp_cspline,Ncspl);
    gsl_spline_init(nc_spline,xx,yy,Ncspl);
    bool_init_nc=true;
    nc_mmin=xx[0];
    nc_mmax=xx[Ncspl-1];
}

void hod::init_Ns_spl(double xx[],double yy[],int Nsspl){
    if(!bool_init_ns){
        gsl_interp_accel_free(ns_acc);
        gsl_spline_free(ns_spline);
        bool_init_ns=false;
    }
    ns_acc=gsl_interp_accel_alloc();
    ns_spline=gsl_spline_alloc(gsl_interp_cspline,Nsspl);
    gsl_spline_init(ns_spline,xx,yy,Nsspl);
    bool_init_ns=true;
    ns_mmin=xx[0];
    ns_mmax=xx[Nsspl-1];
}
#elif TINK==3
// This uses a modified prescription of Martin White's 2012 model
// White -> This model
// Mcut -> Mmin
// kMcut -> Mcut
// M1 -> Msat
// alp -> alpsat
// fac -> central occupation maximum
///And now all the different functions 
///Average number of central galaxies in a halo of mass m
double hod::ncen(double xm)
{
    double arg=(xm-hodp.Mmin)/hodp.siglogM;
    double res=hodp.fac*0.5*(1+gsl_sf_erf(arg));
    return res;
}

///Average number of satellite galaxies in a halo of mass m
double hod::nsat(double xm)
{
    if(xm<hodp.Mcut){
        return 0;
    }else{
        double res=pow(pow(10.,xm-hodp.Msat)-pow(10.,hodp.Mcut-hodp.Msat),hodp.alpsat);
        return res;
    }
}
#elif TINK==4
// This uses a modified prescription of Martin White's 2012 model
// White -> This model
// Mcut -> Mmin
// kMcut -> Mcut -> Mass at which the central occupation becomes one. Use this and fac to calculate the slope
// M1 -> Msat
// alp -> alpsat
// fac -> central occupation maximum at Mmin
///And now all the different functions 
///Average number of central galaxies in a halo of mass m
double hod::ncen(double xm)
{
    double arg=(xm-hodp.Mmin)/hodp.siglogM;
    double slp=(1.-hodp.fac)/(hodp.Mcut-hodp.Mmin);
    double factor=hodp.fac+(xm-hodp.Mmin)*slp;
    if(factor<0.0) factor=0.0;
    if(factor>1.0) factor=1.0;
    double res=factor*0.5*(1+gsl_sf_erf(arg));
    return res;
}

///Average number of satellite galaxies in a halo of mass m
double hod::nsat(double xm)
{
    double res=pow(pow(10.,xm-hodp.Msat),hodp.alpsat);
    return res;
}
#else
// This can use the Martin White prescription from White et al. 2012
// White -> This model
// Mcut -> Mmin
// kMcut -> Mcut
// M1 -> Msat
// alp -> alpsat
///And now all the different functions 
///Average number of central galaxies in a halo of mass m
double hod::ncen(double xm)
{
    double arg=(xm-hodp.Mmin)/hodp.siglogM;
    double res=0.5*(1+gsl_sf_erf(arg));
    double inc=1.0;
    if(xm<inc_xM){
        inc=1.0+inc_alp*(xm-inc_xM);
        if(inc>1.0)inc=1.0; // Guarantees the code works even if inc_alp is negative, effectively ignoring inc_alp in that case
        if(inc<0.0)inc=0.0; // Guarantees the code works even if inc_alp is negative, effectively ignoring inc_alp in that case
    }
    return res*inc;
}

///Average number of satellite galaxies in a halo of mass m
double hod::nsat(double xm)
{
    //if(xm<hodp.Mcut||xm<hodp.Mmin){
    if(xm<hodp.Mcut){
        return 0;
    }else{
        double arg=(xm-hodp.Mmin)/hodp.siglogM;
        double res=0.5*(1+gsl_sf_erf(arg));
        double inc=1.0;
        if(xm<inc_xM){
            inc=1.0+inc_alp*(xm-inc_xM);
            if(inc>1.0)inc=1.0; // Guarantees the code works even if inc_alp is negative, effectively ignoring inc_alp in that case
        }
        res=res*inc;
        res=res*pow(pow(10.,xm-hodp.Msat)-pow(10.,hodp.Mcut-hodp.Msat),hodp.alpsat);
        return res;
    }
}
#endif

///Average number of satellite galaxies 
double hod::nsatz(double z)
{

    double result=0.0;
    for(int i=0;i<N9_16;i++)
    {
        double m=pow(10.,x9_16[i]);
        result=result+log(10.0)*nsat(x9_16[i])*m*nofm(m,z)*w9_16[i];
    }
    return result;
}

///Average number of central galaxies at redshift z
double hod::ncenz(double z)
{

    double result=0.0;
    for(int i=0;i<N9_16;i++)
    {
        double m=pow(10.,x9_16[i]);
        result=result+log(10.0)*ncen(x9_16[i])*m*nofm(m,z)*w9_16[i];
    }
    return result;
}

///h value
double hod::gets8()
{
    return sigma8;
}
double hod::geth()
{
    return h;
}

///Omb value
double hod::getOmb()
{
    return Omegab;
}

///Omk value
double hod::getOmk()
{
    return 1.-Omegal-Omega0;
}

/// Set new redshift if not same as the earlier one
/// Set power spectra to zero.
void hod::resetz(double z)
{
    if(z_glob!=z) setnew_z(z);
    if(bool_init_D2gg){
	fprintf(stderr,"Attempting to free 1\n");
        gsl_interp_accel_free(D2gg_acc);
        gsl_spline_free(D2gg_spline);
        bool_init_D2gg=false;
    }

    if(bool_init_D2gd){
	fprintf(stderr,"Attempting to free 2\n");
        gsl_interp_accel_free(D2gd_acc);
        gsl_spline_free(D2gd_spline);
        bool_init_D2gd=false;
    }

    if(bool_init_xigg){
	fprintf(stderr,"Attempting to free 3\n");
        gsl_interp_accel_free(xigg_acc);
        gsl_spline_free(xigg_spline);
        bool_init_xigg=false;
    }

    if(bool_init_xigd){
	fprintf(stderr,"Attempting to free 4\n");
        gsl_interp_accel_free(xigd_acc);
        gsl_spline_free(xigd_spline);
        bool_init_xigd=false;
    }

    if(bool_init_xigg_bar){
	fprintf(stderr,"Attempting to free 5\n");
        gsl_interp_accel_free(xigg_bar_acc);
        gsl_spline_free(xigg_bar_spline);
        bool_init_xigg_bar=false;
    }

    if(bool_init_xigg_barbar){
	fprintf(stderr,"Attempting to free 6\n");
        gsl_interp_accel_free(xigg_barbar_acc);
        gsl_spline_free(xigg_barbar_spline);
        bool_init_xigg_barbar=false;
    }
    fKaiser=-99.0;

}

/// Average galaxy bias of the sample
/// exclusion and radial dependence of the bias.
double hod::galaxy_bias(double z)
{
    double btot=0.0;
    double totNc=0.0;
    double totNs=0.0;
    for(int i=0;i<N9_16;i++)
    {
        // Numerator HOD factors
        double mass=pow(10.,x9_16[i]);
        double avNcen=ncen(x9_16[i]);
        double avNsat=nsat(x9_16[i]);
        double nofmdm=nofm(mass,z)*w9_16[i]*mass*log(10.0);
        double bmnmdm=bias(mass,z)*nofm(mass,z)*w9_16[i]*mass*log(10.0);

        totNc=totNc+avNcen*nofmdm;
        totNs=totNs+avNsat*nofmdm;
        btot=btot+(avNcen+avNsat)*bmnmdm;

    }
    btot=btot/(totNc+totNs);
    return btot;
}

/// Halo model for gg power spectrum and gd power spectrum, with halo
/// exclusion and radial dependence of the bias.
double hod::Pk_gg_gd_he(double z)
{

    if(verbose){
        std::cout<<"# Calling Pk_gg_gd_he\n";
    }

    // Set the k-array
    double karr[kbins], Pk_gg[kbins], Pk_gd[kbins];
    double xx[kbins],yy[kbins],zz[kbins];
    double hod_kdiff=(hod_kmax-hod_kmin)/(kbins-1.);

    // First the various denominators. 
    //  The rho is in comoving units, so no z dependence.
    double rho_aver=Omega0*rho_crit_0;

    // Set up the mass dependent factors in the integral
    double mdep_1hcs[N9_16],mdep_1hss[N9_16],mdep_1hcd[N9_16],mdep_1hsd[N9_16];
    double rvir[N9_16],mvir[N9_16],cvir[N9_16];
    double r200[N9_16],m200[N9_16],c200[N9_16];
    double fofM[N9_16];

    double mdep_2hs_1[N9_16],mdep_2hd_1[N9_16],mdep_2hc_1[N9_16];
    double mdep_2hs_2[N9_16],mdep_2hd_2[N9_16],mdep_2hc_2[N9_16];
    double mdep_2h_bcorr=0.0;
    double mdep_2h_mcorr=0.0;

    double totNc=0.0;
    double totNs=0.0;
    double btot=0.0;
    for(int i=0;i<N9_16;i++)
    {
        // Numerator HOD factors
        double mass=pow(10.,x9_16[i]);
        double avNcen=ncen(x9_16[i]);
        double avNsat=nsat(x9_16[i]);
        double nofmdm=nofm(mass,z)*w9_16[i]*mass*log(10.0);
        double bmnmdm=bias(mass,z)*nofm(mass,z)*w9_16[i]*mass*log(10.0);

        m200[i]=x9_16[i];
        modelNFWhalo_com(mass, z, mvir[i], rvir[i], cvir[i], r200[i],c200[i]);

        // 1 halo M dependent terms that need a kterm multiplication
        mdep_1hcs[i]=avNcen*avNsat*nofmdm;
        mdep_1hss[i]=avNsat*avNsat*nofmdm;
        mdep_1hcd[i]=avNcen* mass *nofmdm/rho_aver;
        mdep_1hsd[i]=avNsat* mass *nofmdm/rho_aver;

        mdep_2hc_1[i] =avNcen*nofmdm;
        mdep_2hs_1[i] =avNsat*nofmdm;
        mdep_2hd_1[i] = mass *nofmdm/rho_aver;

        mdep_2hc_2[i] =avNcen*bmnmdm;
        mdep_2hs_2[i] =avNsat*bmnmdm;
        mdep_2hd_2[i] = mass *bmnmdm/rho_aver;

        /// 2 halo integral correction items
        mdep_2h_bcorr=mdep_2h_bcorr+mass*bmnmdm/rho_aver;
        mdep_2h_mcorr=mdep_2h_mcorr+mass*nofmdm/rho_aver;

        totNc=totNc+mdep_2hc_1[i];
        totNs=totNs+mdep_2hs_1[i];
        btot=btot+(avNcen+avNsat)*bmnmdm;
        //printf("DEBUG: %e %e %e %e \n",totNc,totNs,avNcen,nofmdm);

        if(fgm_m0>0.0){
            if(fgm_slp!=-99.0){
                fofM[i]=0.5+fgm_slp*(x9_16[i]-fgm_m0);
                if(fofM[i]>1.0) fofM[i]=1.0;
                if(fofM[i]<0.0) fofM[i]=0.0;
            }else{
                fofM[i]=(x9_16[i]<fgm_m0)?0.0:1.0;
            }
        }

    }
    mdep_2h_bcorr=1.0-mdep_2h_bcorr;
    mdep_2h_mcorr=1.0-mdep_2h_mcorr;

    if(fgm_m0>0.0){
        if(fgm_slp>0.0 || fgm_slp==-99.0){
            mdep_2h_bcorr=0.0;
            mdep_2h_mcorr=0.0;
            //printf("# I am assuming no correction due to low mass haloes as f(M) slope is greater than zero\n");
        }
    }
    btot=btot/(totNc+totNs);

    /// Set the Kaiser factor
    fKaiser=-dlnDdln1pz(z);
    fKaiser=fKaiser/btot;

    // Central satellite fractions
    double fcen=totNc/(totNc+totNs);
    double fsat=1.-fcen;
    if(fcen>1.0){
        printf("Error in fcen: %e %e %e %e\n",fcen,fsat,totNc,totNs);
        exit(101);
    }

    for(int i=0;i<N9_16;i++)
    {
        mdep_1hcs[i]= mdep_1hcs[i]/totNc/totNs;
        mdep_1hss[i]= mdep_1hss[i]/totNs/totNs;
        mdep_1hcd[i]= mdep_1hcd[i]/totNc;
        mdep_1hsd[i]= mdep_1hsd[i]/totNs;
        mdep_2hc_1[i]=mdep_2hc_1[i]/totNc;
        mdep_2hs_1[i]=mdep_2hs_1[i]/totNs;
        mdep_2hc_2[i]=mdep_2hc_2[i]/totNc;
        mdep_2hs_2[i]=mdep_2hs_2[i]/totNs;
    }
    //exit(101);

    if(!bool_initQk){
        initQk(z,r200);
    }

    // k and M dependent parts
    int kk=0; /// -> Keeps account of how many were filled
    for(int k=0;k<kbins;k++)
    {

        karr[kk]=hod_kmin+k*hod_kdiff;
        xx[kk]=karr[kk];
        karr[kk]=pow(10.,karr[kk]);

        /// 1 halo terms, single integration
        double int_1hcs,int_1hss,int_1hcd,int_1hsd;
        int_1hcs=int_1hss=int_1hcd=int_1hsd=0.0;
        double uk_s[N9_16],uk_d[N9_16];
        double uk_cen[N9_16];
        for(int i=0;i<N9_16;i++)
        {
#if OFFNEW==1
/// Note fcen_off in this new parameterization is not used, off_rbyrs is in units of rs or 100 hinv kpc
            //uk_cen[i]=exp(-pow(karr[kk]*off_rbyrs*r200[i]/c200[i],2.)/2.);
            //uk_cen[i]=exp(-pow(karr[kk]*off_rbyrs*0.1,2.)/2.);
            uk_cen[i]=1.0-fcen_off+fcen_off*exp(-pow(karr[kk]*off_rbyrs*0.10,2.)/2.);
#else
            if(off_rbyrs==0.0 || fcen_off==0.0){
                uk_cen[i]=1.0;
            }else{
                uk_cen[i]=1.0-fcen_off+fcen_off*exp(-pow(karr[kk]*off_rbyrs*r200[i]/c200[i],2.)/2.);
            }
#endif

            //printf("%f %f %f:\n",karr[k],r200[i],c200[i]);
            uk_d[i]=ukinterp(karr[kk]*r200[i]/c200[i], c200[i]);
            //fprintf(stderr, "uk_d: %le %le %le \n", karr[kk]*r200[i]/c200[i], c200[i], uk_d[i]);
            //printf("DEBUG: %e %e %e \n",karr[k]*r200[i]/c200[i], c200[i],uk_d[i]);

            if(hodp.csbycdm==1.0){
                uk_s[i]=uk_d[i];
            }else{
                //printf("DEBUG: I should not be here, %e \n",hodp.csbycdm);
                uk_s[i]=ukinterp(karr[kk]*r200[i]/(c200[i]*hodp.csbycdm),hodp.csbycdm*c200[i]);
            }

            //double rfac=0.080*pow(10.,1./3.*(x9_16[i]-12.0))/r200[i];
            double rfac=1.0;
            if(rfac!=1.0){
                //printf("DEBUG: I should not be here 2 \n");
                uk_d[i]=ukinterp(karr[kk]*r200[i]/c200[i], c200[i]*rfac);
            }
            //exit(101);

            // Now the integrals for the 1 halo term
            if (miyatake21switch){
                double avNcen = ncen(x9_16[i]);
                if (avNcen>0){
		    int_1hcs+=mdep_1hcs[i]*uk_s[i]*uk_cen[i]/ncen(x9_16[i]); /// --> Change here Divide by ncen(M)
		    int_1hss+=mdep_1hss[i]*pow(uk_s[i],2.)/ncen(x9_16[i]);  /// --> Change here Divide by ncen(M)
                }
            }else{
		int_1hcs+=mdep_1hcs[i]*uk_s[i]*uk_cen[i];
		int_1hss+=mdep_1hss[i]*pow(uk_s[i],2.);  
            }
            if(fgm_m0<0.0){
                int_1hcd+=mdep_1hcd[i]*uk_d[i]*uk_cen[i];
                int_1hsd+=mdep_1hsd[i]*uk_s[i]*uk_d[i];
            }else{
                int_1hcd+=mdep_1hcd[i]*uk_d[i]*fofM[i]*uk_cen[i];
                int_1hsd+=mdep_1hsd[i]*uk_s[i]*uk_d[i]*fofM[i];
            }

        }

        /// 2 halo terms, double integration
        double int_2hcs,int_2hss,int_2hcc;
        double int_2hcd,int_2hsd;
        int_2hcc=int_2hcs=int_2hss=int_2hcd=int_2hsd=0.0;
        for(int i=0;i<N9_16;i++)
        {
            for(int j=0;j<N9_16;j++)
            {
                // Now the double integrals, assuming that both the Qk
                // factors have been initialized
                int_2hcc +=(mdep_2hc_1[i]*mdep_2hc_1[j]*uk_cen[i]*uk_cen[j]*Qk1[(k*N9_16+i)*N9_16+j]);
                int_2hcc +=(mdep_2hc_2[i]*mdep_2hc_2[j]*uk_cen[i]*uk_cen[j]*Qk2[(k*N9_16+i)*N9_16+j]);

                int_2hcs +=(mdep_2hc_1[i]*mdep_2hs_1[j]*uk_cen[i]*uk_s[j]*Qk1[(k*N9_16+i)*N9_16+j]);
                int_2hcs +=(mdep_2hc_2[i]*mdep_2hs_2[j]*uk_cen[i]*uk_s[j]*Qk2[(k*N9_16+i)*N9_16+j]);

                int_2hss +=(mdep_2hs_1[i]*uk_s[i]*mdep_2hs_1[j]*uk_s[j]*Qk1[(k*N9_16+i)*N9_16+j]);
                int_2hss +=(mdep_2hs_2[i]*uk_s[i]*mdep_2hs_2[j]*uk_s[j]*Qk2[(k*N9_16+i)*N9_16+j]);

                // Now the double integrals for gd, assuming that both the Qk
                // factors have been initialized
                if(fgm_m0<0.0){
                    int_2hcd +=(mdep_2hc_1[i]*mdep_2hd_1[j]*uk_cen[i]*uk_d[j]*Qk1[(k*N9_16+i)*N9_16+j]);
                    int_2hcd +=(mdep_2hc_2[i]*mdep_2hd_2[j]*uk_cen[i]*uk_d[j]*Qk2[(k*N9_16+i)*N9_16+j]);

                    int_2hsd +=(mdep_2hs_1[i]*uk_s[i]*mdep_2hd_1[j]*uk_d[j]*Qk1[(k*N9_16+i)*N9_16+j]);
                    int_2hsd +=(mdep_2hs_2[i]*uk_s[i]*mdep_2hd_2[j]*uk_d[j]*Qk2[(k*N9_16+i)*N9_16+j]);
                }else{
                    int_2hcd +=(mdep_2hc_1[i]*mdep_2hd_1[j]*uk_cen[i]*uk_d[j]*Qk1[(k*N9_16+i)*N9_16+j])*fofM[i];
                    int_2hcd +=(mdep_2hc_2[i]*mdep_2hd_2[j]*uk_cen[i]*uk_d[j]*Qk2[(k*N9_16+i)*N9_16+j])*fofM[i];

                    int_2hsd +=(mdep_2hs_1[i]*uk_s[i]*mdep_2hd_1[j]*uk_d[j]*Qk1[(k*N9_16+i)*N9_16+j])*fofM[i];
                    int_2hsd +=(mdep_2hs_2[i]*uk_s[i]*mdep_2hd_2[j]*uk_d[j]*Qk2[(k*N9_16+i)*N9_16+j])*fofM[i];
                }

            }

            /// Correction for the finite range of integration
            int_2hcd +=(mdep_2hc_1[i]*uk_cen[i]*Qk1[(k*N9_16+i)*N9_16+i]*mdep_2h_mcorr);
            int_2hcd +=(mdep_2hc_2[i]*uk_cen[i]*Qk2[(k*N9_16+i)*N9_16+i]*mdep_2h_bcorr);
            int_2hsd +=(mdep_2hs_1[i]*uk_s[i]*Qk1[(k*N9_16+i)*N9_16+i]*mdep_2h_mcorr);
            int_2hsd +=(mdep_2hs_2[i]*uk_s[i]*Qk2[(k*N9_16+i)*N9_16+i]*mdep_2h_bcorr);

        }

        // Power spectra various terms from all the summations

        double kfac=pow(karr[kk],3.)/(2.0*M_PI*M_PI);

        //printf("k: %e 1hcd: %e 1hsd: %e \n",karr[kk],int_1hcd,int_1hsd);

        double P_gg_1hcs=int_1hcs;
        double P_gg_1hss=int_1hss;
        double P_gg_2hcc=int_2hcc;
        double P_gg_2hcs=int_2hcs;
        double P_gg_2hss=int_2hss;

        double P_gd_1hc =int_1hcd;
        double P_gd_1hs =int_1hsd;
        double P_gd_2hc =int_2hcd; 
        double P_gd_2hs =int_2hsd; 

        /// Everything below here stays the same!!!
        // Total gg and gd power spectra
        //std::cout<<whichopt<<" whichopt is";
        if(whichopt==0){
            Pk_gg[kk]=2.0*fcen*fsat*P_gg_1hcs+fsat*fsat*P_gg_1hss+fcen*fcen*P_gg_2hcc+2.0*fcen*fsat*P_gg_2hcs+fsat*fsat*P_gg_2hss;
            //fprintf(stderr, "%le %le %le %le %le %le %le %le\n", karr[kk], fcen, fsat, P_gg_1hcs, P_gg_1hss, P_gg_2hcc, P_gg_2hcs, P_gg_2hss);
        }else{
            Pk_gg[kk]=0.0;
            if(whichopt==1){
                Pk_gg[kk]=2.0*fcen*fsat*P_gg_1hcs;
                //Pk_gg[kk]=fcen*fcen*P_gg_2hcc;
            }else if(whichopt==2){
                Pk_gg[kk]=fsat*fsat*P_gg_1hss;
                //Pk_gg[kk]=2.0*fcen*fsat*P_gg_2hcs;
            }else if(whichopt==3){
                Pk_gg[kk]=fcen*fcen*P_gg_2hcc+2.0*fcen*fsat*P_gg_2hcs+fsat*fsat*P_gg_2hss;
                //Pk_gg[kk]=fsat*fsat*P_gg_2hss;
            }
        }

        if(whichopt==0){
            Pk_gd[kk]=fcen*P_gd_1hc+fsat*P_gd_1hs+fcen*P_gd_2hc+fsat*P_gd_2hs;
        }else{
            Pk_gd[kk]=0.0;
            if(whichopt==1){
                Pk_gd[kk]=fcen*P_gd_1hc;
            }else if(whichopt==2){
                Pk_gd[kk]=fsat*P_gd_1hs;
            }else if(whichopt==3){
                Pk_gd[kk]=fcen*P_gd_2hc+fsat*P_gd_2hs;
            }
        }

        if(Pk_gg[kk]<=0||Pk_gd[kk]<=0){
            std::cout<<"# Power spectra negative, will skip this k-point in the interpolation "<<std::endl;
            printf("# xx:%e Pk_gg:%e Pk_gd:%e Pgg1hcs:%e Pgg1hss:%e Pgg2hcc:%e Pgg2hcs:%e Pgg2hss:%e fcen:%e fsat:%e \n",xx[kk],Pk_gg[kk],Pk_gd[kk],P_gg_1hcs,P_gg_1hss,P_gg_2hcc,P_gg_2hcs,P_gg_2hss,fcen,fsat);
            print();
	    continue;
            if(Pk_gg[kk]<0){
                yy[kk]=-15.0;
            }else{
                yy[kk]=log10(Pk_gg[kk]*kfac);
            }
            if(Pk_gd[kk]<0){
                zz[kk]=-15.0;
            }else{
                zz[kk]=log10(Pk_gd[kk]*kfac);
            }
        }else{
            yy[kk]=log10(Pk_gg[kk]*kfac);
            zz[kk]=log10(Pk_gd[kk]*kfac);
	    if(verbose)
		printf("%e %e %e\n",xx[kk],Pk_gg[kk],Pk_gd[kk]);
	    kk++;
            //printf("%d %le %le\n",k,Pk_gg[k],Pk_gd[k]);
            //printf("%e %e %e\n",xx[k],yy[k],zz[k]);
        }


    } // Loop over k
    //exit(101);

    if(kk!=kbins && verbose)
	printf("# Only %d valid kk values found to be used for spline, proceeding with crossed fingers \n",kk);

    // Setup the splines for D2_gg
    D2gg_acc = gsl_interp_accel_alloc ();
    D2gg_spline = gsl_spline_alloc (gsl_interp_cspline, kk);
    gsl_spline_init (D2gg_spline,xx,yy,kk);
    D2gg_dlow=(yy[2]-yy[0])/(xx[2]-xx[0]);
    D2gg_xlow=xx[1]; D2gg_ylow=yy[1];
    D2gg_dhigh=(yy[kk-3]-yy[kk-1])/(xx[kk-3]-xx[kk-1]);
    D2gg_xhigh=xx[kk-2]; D2gg_yhigh=yy[kk-2];
    bool_init_D2gg=true;

    // Setup the splines for D2_gd
    D2gd_acc = gsl_interp_accel_alloc ();
    D2gd_spline = gsl_spline_alloc (gsl_interp_cspline, kk);
    gsl_spline_init (D2gd_spline,xx,zz,kk);
    D2gd_dlow=(zz[2]-zz[0])/(xx[2]-xx[0]);
    D2gd_xlow=xx[1]; D2gd_ylow=zz[1];
    D2gd_dhigh=(zz[kk-3]-zz[kk-1])/(xx[kk-3]-xx[kk-1]);
    D2gd_xhigh=xx[kk-2]; D2gd_yhigh=zz[kk-2];
    bool_init_D2gd=true;

    // Find the place where k P(k) has a maximum
    //int imax=0;
    double xkPk;
    xkPk=-1.0E30;
    for (int i=kk-1;i>0;i--){
        double comp=yy[i]-2.0*log10(karr[i]);
        if(xkPk < comp ){
            xkPk=comp;
            kPkmax=karr[i];
        }
    }
    kPkmax=kPkmax+3.0*M_PI;
    xkPk=-1.0E30;
    for (int i=kk-1;i>0;i--){
        double comp=zz[i]-2.0*log10(karr[i]);
        if(xkPk < comp ){
            xkPk=comp;
            kPkmax2=karr[i];
        }
    }
    kPkmax2=kPkmax2+2.0*M_PI;

    if(verbose){
        std::cout<<"# Call to Pk_gg_gd_he finished\n";
        std::cout<<"# kPkmax kPkmax2 "<<kPkmax<<" "<<kPkmax2<<std::endl;
    }
    return 0;
}

void hod::set_cen_offset_params(double x, double y){
    fcen_off=x;
    off_rbyrs=y;
}

void hod::set_inc_params(double x, double y){
    inc_alp=x;
    inc_xM=y;
}


/// Halo model for gg power spectrum and gd power spectrum
double hod::Pk_gg_gd(double z)
{
    if(verbose){
    std::cout<<"# Calling Pk_gg_gd\n";
    }

    // Sanity check 
    if(fcen_off!=0.0||off_rbyrs!=0.0){
        std::cout<<"The offset central case is not handled while the halo exclusion is not on. Please enable halo exclusion"<<std::endl;
        std::cout<<"I am going to abort execution inside Pk_gg_gd"<<std::endl;
        exit(110);
    }

    // Set the k-array
    double karr[kbins], Pk_gg[kbins], Pk_gd[kbins];
    double xx[kbins],yy[kbins],zz[kbins];
    double hod_kdiff=(hod_kmax-hod_kmin)/(kbins-1.);
    for(int i=0;i<kbins;i++){
        karr[i]=hod_kmin+i*hod_kdiff;
        xx[i]=karr[i];
        karr[i]=pow(10.,karr[i]);
    }

    // First the various denominators. 
    //  The rho is in comoving units, so no z dependence.
    double rho_aver=Omega0*rho_crit_0;
    
    // Set up the mass dependent factors in the integral
    double mdep_1hcs[N9_16],mdep_1hss[N9_16],mdep_1hcd[N9_16],mdep_1hsd[N9_16],mdep_2hs[N9_16],mdep_2hd[N9_16];
    double rvir[N9_16],mvir[N9_16],cvir[N9_16];
    double r200[N9_16],c200[N9_16];

    // No u(k|M) term in the central contribution to the 2halo term. So integrate it directly
    double mdep_2hc=0.0;

    // Correction term for the halo mass correlations
    double mdep_2hcorr=0.0;

    double totNc=0.0;
    double totNs=0.0;
    double btot=0.0;
    for(int i=0;i<N9_16;i++)
    {
	// Numerator HOD factors
	double mass=pow(10.,x9_16[i]);
	double avNcen=ncen(x9_16[i]);
	double avNsat=nsat(x9_16[i]);
	double nofmdm=nofm(mass,z)*w9_16[i]*mass*log(10.0);
	double bmnmdm=bias(mass,z)*nofm(mass,z)*w9_16[i]*mass*log(10.0);
	modelNFWhalo_com(mass, z, mvir[i], rvir[i], cvir[i], r200[i],c200[i]);

	// M dependent terms that need a kterm multiplication
	mdep_1hcs[i]=avNcen*avNsat*nofmdm;
	mdep_1hss[i]=avNsat*avNsat*nofmdm;
	mdep_1hcd[i]=avNcen* mass *nofmdm/rho_aver;
	mdep_1hsd[i]=avNsat* mass *nofmdm/rho_aver;

	mdep_2hs[i] =avNsat*bmnmdm;
	mdep_2hd[i] = mass *bmnmdm/rho_aver;

	// M dependent terms that are k-independent
	mdep_2hc=mdep_2hc+avNcen*bmnmdm;
	mdep_2hcorr=mdep_2hcorr+mass*bmnmdm/rho_aver;

        totNc=totNc+avNcen*nofmdm;
        totNs=totNs+avNsat*nofmdm;

        btot=btot+(avNcen+avNsat)*bmnmdm;

    }

    mdep_2hcorr=1.0-mdep_2hcorr;
    btot=btot/(totNc+totNs);

    //std::cout<<"Nc Ns:"<<totNc<<" "<<totNs<<std::endl;

    // Central satellite fractions
    double fcen=totNc/(totNc+totNs);
    double fsat=1.-fcen;

    /// Set the Kaiser factor
    fKaiser=-dlnDdln1pz(z)/btot;

    for(int i=0;i<N9_16;i++)
    {
	mdep_1hcs[i]=mdep_1hcs[i]/totNc/totNs;
	mdep_1hss[i]=mdep_1hss[i]/totNs/totNs;
	mdep_1hcd[i]=mdep_1hcd[i]/totNc;
	mdep_1hsd[i]=mdep_1hsd[i]/totNs;
                                 
	mdep_2hs[i] =mdep_2hs[i]/totNs;

    }
    mdep_2hc=mdep_2hc/totNc;

    // k and M dependent parts
    for(int k=0;k<kbins;k++)
    {
	double int_1hcs,int_1hss,int_1hcd,int_1hsd,int_2hs,int_2hd;
	int_1hcs=int_1hss=int_1hcd=int_1hsd=int_2hs=int_2hd=0.0;
	for(int i=0;i<N9_16;i++)
	{
	    double uk_s,uk_d;
	    uk_d=ukinterp(karr[k]*r200[i]/c200[i], c200[i]);
	    if(hodp.csbycdm==1.0){
		uk_s=uk_d;
	    }else{
		uk_s=ukinterp(karr[k]*r200[i]/(c200[i]*hodp.csbycdm), hodp.csbycdm*c200[i]);
	    }

	    // Now the integrals
	    int_1hcs+=mdep_1hcs[i]*uk_s;
	    int_1hss+=mdep_1hss[i]*pow(uk_s,2.);
	    int_1hcd+=mdep_1hcd[i]*uk_d;
	    int_1hsd+=mdep_1hsd[i]*uk_s*uk_d;
	    int_2hs +=mdep_2hs[i] *uk_s;
	    int_2hd +=mdep_2hd[i] *uk_d;

	    //std::cout<<std::scientific;
	    //std::cout<<karr[k]<<" "<<c200[i]<<" "<<uk_s<<std::endl;
	}

	// Power spectra various terms from all the summations

	double kfac=pow(karr[k],3.)/(2.0*M_PI*M_PI);
	double pofk=Delta2_L_num(karr[k],z)/kfac;

	double P_gg_1hcs=int_1hcs;
	double P_gg_1hss=int_1hss;
	double P_gg_2hcc=mdep_2hc*mdep_2hc*pofk;
	double P_gg_2hcs=mdep_2hc* int_2hs*pofk;
	double P_gg_2hss=int_2hs * int_2hs*pofk;

	double P_gd_1hc =int_1hcd;
	double P_gd_1hs =int_1hsd;
	double P_gd_2hc =mdep_2hc*(int_2hd+mdep_2hcorr)*pofk; 
	double P_gd_2hs =int_2hs *(int_2hd+mdep_2hcorr)*pofk; 

	// Total gg and gd power spectra
	Pk_gg[k]=2.0*fcen*fsat*P_gg_1hcs+fsat*fsat*P_gg_1hss+fcen*fcen*P_gg_2hcc+2.0*fcen*fsat*P_gg_2hcs+fsat*fsat*P_gg_2hss;
	Pk_gd[k]=fcen*P_gd_1hc+fsat*P_gd_1hs+fcen*P_gd_2hc+fsat*P_gd_2hs;

	if(Pk_gg[k]<=0||Pk_gd[k]<=0){
	    std::cout<<"Power spectra negative, please debug"<<std::endl;
	    exit(0);
	}else{
	    yy[k]=log10(Pk_gg[k]*kfac);
	    zz[k]=log10(Pk_gd[k]*kfac);
	}

    }

    // Setup the splines for D2_gg
    D2gg_acc = gsl_interp_accel_alloc ();
    D2gg_spline = gsl_spline_alloc (gsl_interp_cspline, kbins);
    gsl_spline_init (D2gg_spline,xx,yy,kbins);
    D2gg_dlow=(yy[2]-yy[0])/(xx[2]-xx[0]);
    D2gg_xlow=xx[1]; D2gg_ylow=yy[1];
    D2gg_dhigh=(yy[kbins-3]-yy[kbins-1])/(xx[kbins-3]-xx[kbins-1]);
    D2gg_xhigh=xx[kbins-2]; D2gg_yhigh=yy[kbins-2];
    bool_init_D2gg=true;

    // Setup the splines for D2_gd
    D2gd_acc = gsl_interp_accel_alloc ();
    D2gd_spline = gsl_spline_alloc (gsl_interp_cspline, kbins);
    gsl_spline_init (D2gd_spline,xx,zz,kbins);
    D2gd_dlow=(zz[2]-zz[0])/(xx[2]-xx[0]);
    D2gd_xlow=xx[1]; D2gd_ylow=zz[1];
    D2gd_dhigh=(zz[kbins-3]-zz[kbins-1])/(xx[kbins-3]-xx[kbins-1]);
    D2gd_xhigh=xx[kbins-2]; D2gd_yhigh=zz[kbins-2];
    bool_init_D2gd=true;

    // Find the place where k P(k) has a maximum
    for (int i=kbins-1;i>0;i--){
	if(yy[i]-2.0*log10(karr[i]) > yy[i-1]-2.0*log10(karr[i-1]) )
	{
	    kPkmax=karr[i];
	    break;
	}
    }
    kPkmax=kPkmax+3.0*M_PI;
    for (int i=kbins-1;i>0;i--){
	if(zz[i]-2.0*log10(karr[i]) > zz[i-1]-2.0*log10(karr[i-1]) )
	{
	   kPkmax2=karr[i];
	    break;
	}
    }
    kPkmax2=kPkmax2+2.0*M_PI;

    if(verbose){
    std::cout<<"# Call to Pk_gg_gd finished\n";
    std::cout<<"# "<<kPkmax<<" "<<kPkmax2<<std::endl;
    }
    return 0;
}

/// Numerical calculation of Delta_gg
double hod::D2gg_num(double k, double z)
{
    if(!bool_init_D2gg)
    {
        if(halo_exc){
            Pk_gg_gd_he(z);
        }else{
            Pk_gg_gd(z);
        }
    }
    double logk=log10(k);
    double result;
    if(logk>hod_kmin&&logk<hod_kmax)
    {
        //result=gsl_interp_eval (D2gg_interp, D2gg_x, D2gg_y, logk, D2gg_acc);
	//std::cout<<logk<<" "<<result<<" ";
        result=gsl_spline_eval (D2gg_spline,logk, D2gg_acc);
	//std::cout<<result<<" "<<std::endl;
    } else if (logk<=hod_kmin)
    {
	//Extrapolate now
	result=D2gg_ylow+D2gg_dlow*(logk-D2gg_xlow);
    } else if (logk>=hod_kmax)
    {
	result=D2gg_yhigh+D2gg_dhigh*(logk-D2gg_xhigh);
    }

    result = pow(10.,result);
    return result;

}

/// Numerical calculation of Delta_gd
double hod::D2gd_num(double k, double z)
{
    if(!bool_init_D2gd)
    {
        if(halo_exc){
            Pk_gg_gd_he(z);
        }else{
            Pk_gg_gd(z);
        }
    }
    double logk=log10(k);
    double result;
    if(logk>hod_kmin&&logk<hod_kmax)
    {
        result=gsl_spline_eval (D2gd_spline,logk, D2gd_acc);
    } else if (logk<=hod_kmin)
    {
        //Extrapolate now
        result=D2gd_ylow+D2gd_dlow*(logk-D2gd_xlow);
    } else if (logk>=hod_kmax)
    {
        result=D2gd_yhigh+D2gd_dhigh*(logk-D2gd_xhigh);
    }
    result = pow(10.,result);
    return result;

}

/// Correlation functions and their interpolations
/// Friend function for xi_gd
double dxi_gd(double x, void* params)
{
    /// Int_{0}^{infty} Delta(x/r)/(x*x)*sin(x) dx
    xi_gd_params c1 = *(xi_gd_params *)params;

    hod *c2;
    double *r;
    double *z;

    c2=c1.hptr;
    r=c1.r;
    z=c1.z;

    double arg=x/(*r);
    double res=(*c2).D2gd_num(arg,*z)/(x*x);

    //std::cout<<"Check "<<x<<" "<<res<<std::endl;

    return res;
}

/// Friend function for xi_gg
double dxi_gg(double x, void* params)
{
    /// Int_{0}^{infty} Delta(x/r)/(x*x)*sin(x) dx
    xi_gg_params c1 = *(xi_gg_params *)params;

    hod *c2;
    double *r;
    double *z;

    c2=c1.hptr;
    r=c1.r;
    z=c1.z;

    double arg=x/(*r);
    double res=(*c2).D2gg_num(arg,*z)/(x*x);

    return res;
}

/// g-g and g-d correlation functions as FT(P_gg) and FT(P_gd)
double hod::xi_gg_gd(double z)
{
    if(verbose){
        std::cout<<"# Calling xi_gg_gd\n";
    }
    // Define rbins 
    // 240 req accuracy 1e-3 in projection integral
    // 80 req accuracy 1e-3 in projection integral
    static const int rbins=40;
    double rarr[rbins];
    double rdiff=(hod_rmax-hod_rmin)/(rbins-1.);
    double xx[rbins],yy[rbins],zz[rbins];
    double xigg[rbins],xigd[rbins];

    for(int i=0;i<rbins;i++)
    {
        rarr[i]=hod_rmin+rdiff*i;
        xx[i]=rarr[i];
        rarr[i]=pow(10.,rarr[i]);
    }

    if(!bool_init_D2gg||!bool_init_D2gd)
    {
        if(halo_exc){
            Pk_gg_gd_he(z);
        }else{
            Pk_gg_gd(z);
        }
    }

    //std::cout<<kPkmax<<" "<<kPkmax2<<std::endl;
    gsl_integration_qawo_table * qt = gsl_integration_qawo_table_alloc (1.0,3.0*M_PI,GSL_INTEG_SINE,1000);
    /// Fourier transform D2gg
    for(int i=0;i<rbins;i++)
    {
        double result0,result,error;
        result=result0=0.0;
        gsl_function F;
        F.function=&(dxi_gg);

        xi_gg_params p;
        p.hptr =this;
        p.r=&(rarr[i]);
        p.z=&z;

        F.params=&p;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_integration_workspace * cw = gsl_integration_workspace_alloc (1000);
        gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
        int status=0;

        if(!mock){
            status=gsl_integration_qawf (&F, 1.E-30, 1.e-3, 1000, w, cw, qt, &result, &error); 
        }else{
            //    status=gsl_integration_qawf (&F, 4.35E-2*rarr[i], 1.e-3, 1000, w, cw, qt, &result, &error); 
            status=gsl_integration_qawf (&F, 4.35E-2*rarr[i]/4., 1.e-3, 1000, w, cw, qt, &result, &error); 
        }
        gsl_set_error_handler(oldhand);
        if(status!=0 && verbose){
            printf("# Problems in convergence, xigg: z:%e r:%e\n",z,rarr[i]);
        }

        gsl_integration_workspace_free (w);
        gsl_integration_workspace_free (cw);

        yy[i]=log10(1.+result+result0);
        //std::cout<<rarr[i]<<" "<<result<<std::endl;
        xigg[i]=result+result0;

    }
    //gsl_integration_qawo_table_free (qt0);
    //exit(0);

    // Initialise the splines
    xigg_acc = gsl_interp_accel_alloc ();
    xigg_spline = gsl_spline_alloc (gsl_interp_cspline, rbins);
    gsl_spline_init (xigg_spline,xx,yy,rbins);
    xigg_dlow=(yy[2]-yy[0])/(xx[2]-xx[0]);
    xigg_xlow=xx[1]; xigg_ylow=yy[1];
    xigg_dhigh=(yy[rbins-3]-yy[rbins-1])/(xx[rbins-3]-xx[rbins-1]);
    xigg_xhigh=xx[rbins-2]; xigg_yhigh=yy[rbins-2];
    bool_init_xigg=true;
    if(verbose){
        std::cout<<"# xi_gg inited \n";
    }

    /// Fourier transform D2gd
    for(int i=0;i<rbins;i++)
    {
        double result0,result,error;
        result=result0=0.0;
        gsl_function F;
        F.function=&(dxi_gd);

        xi_gd_params p;
        p.hptr =this;
        p.r=&(rarr[i]);
        p.z=&z;

        F.params=&p;

        //std::cout<<"Done 1 "<<kPkmax2<<std::endl;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_integration_workspace * cw = gsl_integration_workspace_alloc (1000);
        gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
        int status=0;

        if(!mock){
            gsl_integration_qawf (&F, 1.E-30, 1.e-3, 1000, w, cw, qt, &result, &error); 
        }else{
            //gsl_integration_qawf (&F, 4.35E-2*rarr[i], 1.e-3, 1000, w, cw, qt, &result, &error); 
            gsl_integration_qawf (&F, 4.35E-2*rarr[i]/4., 1.e-3, 1000, w, cw, qt, &result, &error); 
        }
        gsl_set_error_handler(oldhand);
        if(status!=0 && verbose){
            printf("# Problems in convergence, xigd: z:%e r:%e\n",z,rarr[i]);
        }

        gsl_integration_workspace_free (w);
        gsl_integration_workspace_free (cw);
        //std::cout<<"Done 2\n";

        zz[i]=log10(1.+result+result0);
        xigd[i]=result+result0;
        //std::cout<<rarr[i]<<" "<<xigg[i]<<" "<<xigd[i]<<std::endl;

    }
    //gsl_integration_qawo_table_free (qt0);
    gsl_integration_qawo_table_free (qt);
    //exit(0);

    // Initialise the splines
    xigd_acc = gsl_interp_accel_alloc ();
    xigd_spline = gsl_spline_alloc (gsl_interp_cspline, rbins);
    gsl_spline_init (xigd_spline,xx,zz,rbins);
    xigd_dlow=(zz[2]-zz[0])/(xx[2]-xx[0]);
    xigd_xlow=xx[1]; xigd_ylow=zz[1];
    xigd_dhigh=(zz[rbins-3]-zz[rbins-1])/(xx[rbins-3]-xx[rbins-1]);
    xigd_xhigh=xx[rbins-2]; xigd_yhigh=zz[rbins-2];
    bool_init_xigd=true;
    if(verbose){
        std::cout<<"# xi_gd inited \n";
        std::cout<<"# Call to xi_gg_gd finished \n";
    }

    return 0;
}

/// g-d correlation functions as FT(P_gd)
double hod::xi_gd(double z)
{
    if(verbose){
        std::cout<<"# Calling xi_gd\n";
    }
    // Define rbins 
    // 240 req accuracy 1e-3 in projection integral
    // 80 req accuracy 1e-3 in projection integral
    static const int rbins=40;
    double rarr[rbins];
    double rdiff=(hod_rmax-hod_rmin)/(rbins-1.);
    double xx[rbins],zz[rbins];
    double xigd[rbins];

    for(int i=0;i<rbins;i++)
    {
        rarr[i]=hod_rmin+rdiff*i;
        xx[i]=rarr[i];
        rarr[i]=pow(10.,rarr[i]);
    }

    if(!bool_init_D2gd)
    {
        if(halo_exc){
            Pk_gg_gd_he(z);
        }else{
            Pk_gg_gd(z);
        }
    }

    //std::cout<<kPkmax<<" "<<kPkmax2<<std::endl;

    gsl_integration_qawo_table * qt = gsl_integration_qawo_table_alloc (1.0,3.0*M_PI,GSL_INTEG_SINE,1000);

    /// Fourier transform D2gd
    for(int i=0;i<rbins;i++)
    {
        double result0,result,error;
        result=result0=0.0;
        gsl_function F;
        F.function=&(dxi_gd);

        xi_gd_params p;
        p.hptr =this;
        p.r=&(rarr[i]);
        p.z=&z;

        F.params=&p;

        //std::cout<<"Done 1 "<<kPkmax2<<std::endl;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_integration_workspace * cw = gsl_integration_workspace_alloc (1000);
        gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
        int status=0;

        if(!mock){
            status=gsl_integration_qawf (&F, 1.E-30, 1.e-3, 1000, w, cw, qt, &result, &error); 
        }else{
            //status=gsl_integration_qawf (&F, 4.35E-2*rarr[i], 1.e-3, 1000, w, cw, qt, &result, &error); 
            status=gsl_integration_qawf (&F, 4.35E-2*rarr[i]/4., 1.e-3, 1000, w, cw, qt, &result, &error); 
        }
        gsl_set_error_handler(oldhand);
        if(status!=0 && verbose){
            printf("# Problems in convergence, xigd: z:%e r:%e\n",z,rarr[i]);
        }

        gsl_integration_workspace_free (w);
        gsl_integration_workspace_free (cw);
        //std::cout<<"Done 2\n";

        zz[i]=log10(1.+result+result0);
        xigd[i]=result+result0;
        //std::cout<<rarr[i]<<" "<<xigd[i]<<std::endl;

    }
    //gsl_integration_qawo_table_free (qt0);
    gsl_integration_qawo_table_free (qt);
    //exit(0);

    // Initialise the splines
    xigd_acc = gsl_interp_accel_alloc ();
    xigd_spline = gsl_spline_alloc (gsl_interp_cspline, rbins);
    gsl_spline_init (xigd_spline,xx,zz,rbins);
    xigd_dlow=(zz[2]-zz[0])/(xx[2]-xx[0]);
    xigd_xlow=xx[1]; xigd_ylow=zz[1];
    xigd_dhigh=(zz[rbins-3]-zz[rbins-1])/(xx[rbins-3]-xx[rbins-1]);
    xigd_xhigh=xx[rbins-2]; xigd_yhigh=zz[rbins-2];
    bool_init_xigd=true;
    if(verbose){
        std::cout<<"# xi_gd inited \n";
        std::cout<<"# Call to xi_gg_gd finished \n";
    }

    return 0;

}

/// g-g and g-d correlation functions as FT(P_gg) and FT(P_gd)
double hod::xi_gg(double z)
{
    if(verbose){
        std::cout<<"# Calling xi_gg\n";
    }
    // Define rbins 
    // 240 req accuracy 1e-3 in projection integral
    // 80 req accuracy 1e-3 in projection integral
    static const int rbins=240;
    double rarr[rbins];
    double rdiff=(hod_rmax-hod_rmin)/(rbins-1.);
    double xx[rbins],yy[rbins];
    double xigg[rbins];

    for(int i=0;i<rbins;i++)
    {
        rarr[i]=hod_rmin+rdiff*i;
        xx[i]=rarr[i];
        rarr[i]=pow(10.,rarr[i]);
    }

    if(!bool_init_D2gg)
    {
        if(halo_exc){
            Pk_gg_gd_he(z);
        }else{
            Pk_gg_gd(z);
        }
    }

    //std::cout<<kPkmax<<" "<<kPkmax2<<std::endl;

    gsl_integration_qawo_table * qt = gsl_integration_qawo_table_alloc (1.0,3.0*M_PI,GSL_INTEG_SINE,1000);
    /// Fourier transform D2gg
    for(int i=0;i<rbins;i++)
    {
        //std::cerr<<rarr[i]<<std::endl;
        double result0,result,error;
        result0=result=0.0;
        gsl_function F;
        F.function=&(dxi_gg);

        xi_gg_params p;
        p.hptr =this;
        p.r=&(rarr[i]);
        p.z=&z;

        F.params=&p;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_integration_workspace * cw = gsl_integration_workspace_alloc (1000);
        gsl_error_handler_t *oldhand=gsl_set_error_handler_off();

        int status=0;
        if(!mock){
            status=gsl_integration_qawf (&F, 1.0E-30, 1.e-3, 1000, w, cw, qt, &result, &error); 
        }else{
            //status=gsl_integration_qawf (&F, 4.35E-2*rarr[i], 1.e-3, 1000, w, cw, qt, &result, &error); 
            status=gsl_integration_qawf (&F, 4.35E-2*rarr[i]/4., 1.e-3, 1000, w, cw, qt, &result, &error); 
        }
        gsl_set_error_handler(oldhand);
        if(status!=0 && verbose){
            printf("# Problems in convergence, xigg: z:%e r:%e\n",z,rarr[i]);
        }

        gsl_integration_workspace_free (w);
        gsl_integration_workspace_free (cw);

        yy[i]=log10(1.+result+result0);
        xigg[i]=result+result0;
        //std::cin>>rarr[i]>>xigg[i];
        yy[i]=log10(1.+xigg[i]);
        //std::cout<<rarr[i]<<" "<<xigg[i]<<std::endl;

    }
    //gsl_integration_qawo_table_free (qt0);
    gsl_integration_qawo_table_free (qt);
    //exit(0);

    // Initialise the splines
    xigg_acc = gsl_interp_accel_alloc ();
    xigg_spline = gsl_spline_alloc (gsl_interp_cspline, rbins);
    gsl_spline_init (xigg_spline,xx,yy,rbins);
    xigg_dlow=(yy[2]-yy[0])/(xx[2]-xx[0]);
    xigg_xlow=xx[1]; xigg_ylow=yy[1];
    xigg_dhigh=(yy[rbins-3]-yy[rbins-1])/(xx[rbins-3]-xx[rbins-1]);
    xigg_xhigh=xx[rbins-2]; xigg_yhigh=yy[rbins-2];
    bool_init_xigg=true;
    if(verbose){
        std::cout<<"# xi_gg inited \n";
    }
    return 0;
}

/// Numerical interpolation for xi_gg
double hod::xigg_num(double r, double z)
{
    if(!bool_init_xigg)
    {
        xi_gg(z);
    }
    double logr=log10(r);
    double result;
    if(logr>hod_rmin&&logr<hod_rmax)
    {
        result=gsl_spline_eval (xigg_spline,logr, xigg_acc);
    } else if (logr<=hod_rmin)
    {
        //Extrapolate now
        result=xigg_ylow+xigg_dlow*(logr-xigg_xlow);
    } else if (logr>=hod_rmax_u)
    {
        result=xigg_yhigh+xigg_dhigh*(logr-xigg_xhigh);
        //result=0.0;
        //return result;
    }

    result = pow(10.,result);
    return result-1.;

}

/// Numerical interpolation for xi_gd
double hod::xigd_num(double r, double z)
{
    if(!bool_init_xigd)
    {
        xi_gd(z);
    }
    double logr=log10(r);
    double result;

    if(logr>hod_rmin&&logr<hod_rmax)
    {
        result=gsl_spline_eval (xigd_spline,logr, xigd_acc);
    } else if (logr<=hod_rmin)
    {
        //Extrapolate now
        result=xigd_ylow+xigd_dlow*(logr-xigd_xlow);
    } else if (logr>=hod_rmax_u)
    {
        //result=xigd_yhigh+xigd_dhigh*(logr-xigd_xhigh);
        //std::cout<<"Set to zero\n";
        result=0.0;
        return result;
    }

    result = pow(10.,result);
    //std::cout<<r<<" "<<result-1.<<std::endl;
    return result-1.;

}

/// Friend function for Wp projection integral
double dwp2(double x, void* params)
{
    /// Wp(R) = 2.* R Int_{0}^{xmax} xi( x R ) x/sqrt(x^2-1.) dx
    /// R is multiplied outside the integral
    proj_params c1 = *(proj_params *)params;

    hod *c2;
    double *R;
    double *z;

    c2=c1.hptr;
    R=c1.R;
    z=c1.z;

    double arg=x*(*R);
    double res=(*c2).xigg_num(arg,*z)*x/sqrt(x*x-1.);

    return res;
}

double dwp(double x, void* params)
{
    /// Wp(R) = 2.* R Int_{0}^{xmax} xi( \sqrt(x^2+1) R ) dx
    /// R is multiplied outside the integral
    proj_params c1 = *(proj_params *)params;

    hod *c2;
    double *R;
    double *z;

    c2=c1.hptr;
    R=c1.R;
    z=c1.z;

    double arg=sqrt(x*x+1.)*(*R);
    double res=(*c2).xigg_num(arg,*z);

    //std::cout<<res<<" "<<x<<std::endl;

    return res;
}

/// Friend function for surface density projection integral
double dsigma(double x, void* params)
{
    /// Sigma(R) = 2. * rhobar* R Int_{0}^{\infty} xi_{gd}( \sqrt(x^2+1) R ) dx
    /// 2. rhobar R is multiplied outside the integral
    proj_params c1 = *(proj_params *)params;

    hod *c2;
    double *R;
    double *z;

    c2=c1.hptr;
    R=c1.R;
    z=c1.z;

    double arg=sqrt(x*x+1.)*(*R);
    //if(arg>10.0) return 0.0;
    double res=(*c2).xigd_num(arg,*z);

    //std::cout<<arg<<" "<<(*c2).xigd_num(arg,*z)<<std::endl;

    return res;
}

/// Projected correlation functions. See the DOCS/NOTES/doc.pdf for
/// documentation in this relation
/// The inputs should be: rpbins_wp, rpbins_esd, wp_Rp, esd_Rp
/// The output should be: wp and esd
double hod::Wp_ESD(double z, int rpbins_wp, int rpbins_esd, double wp_Rp[], double esd_Rp[], double wp[], double esd[], int rpbins_esd2,double pimax,bool reset)
{
    if(reset){
        // Reset the redshift
        resetz(z);
    }

    // We are going to use qagw integration.
    // I need to calculate these at the N projected radius bins from
    // the data. For now just choose some bins for convenience

    double rho_aver=Omega0*rho_crit_0;

    double *xx;
    xx=new double[rpbins_wp];

    for(int i=0;i<rpbins_wp;i++)
    {
        xx[i]=wp_Rp[i];
    }

    for(int i=0;i<rpbins_wp;i++)
    {
        /// QAGW integral
        ///Int_{0}^{ymax} dy 2  R \xi_gg( \sqrt(y^2+1) R)
        /// ymax=pimax/R
        double result, error;
        double ymax=(pimax/xx[i]);

        gsl_function F;
        F.function = &(dwp);

        proj_params p;
        p.hptr = this;
        p.R = &(xx[i]);
        p.z = &z;

        F.params = &p;
        //std::cout<<i<<" YMAX "<<xx[i]<<" "<<pimax/xx[i]<<std::endl;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
        int status=gsl_integration_qag (&F, 0.0, ymax, 0, 1.e-3, 1000, 5, w, &result, &error);
        gsl_set_error_handler(oldhand);
        if(status!=0 && verbose){
            printf("# Problems in convergence, wp: z:%e r:%e\n",z,xx[i]);
        }
        gsl_integration_workspace_free (w);

        wp[i]=2.0*result*xx[i];
        //std::cout<<xx[i]<<" "<<wp[i]<<std::endl;

    }
    //exit(0);
    if(verbose){
        std::cout<<"# Wp calced\n";
    }

    /// Surface density integral
    double *sd,*ww,*xx_esd;
    sd=new double[rpbins_esd2];
    ww=new double[rpbins_esd2];
    xx_esd=new double[rpbins_esd2];

    int first=rpbins_esd2-rpbins_esd;

    for(int i=first;i<rpbins_esd2;i++)
    {
        xx_esd[i]=esd_Rp[i-first];
    }
    double minrp=log10(esd_Rp[0]);
    for(int i=first-1;i>-1;i--)
    {
        xx_esd[i]=pow(10.0,minrp-0.1*(first-i));
    }


    for(int i=0;i<rpbins_esd2;i++)
    {
        //std::cout<<"SD integral"<<xx_esd[i]<<std::endl;
        /// QAGIU integral
        ///Int_{0}^{\infty} dy R \xi_gm( \sqrt(z^2+1) R)
        double result1, result2, error;
        result1=result2=error=0.0;

        gsl_function F;
        F.function = &(dsigma);

        proj_params p;
        p.hptr = this;
        p.R = &(xx_esd[i]);
        p.z = &z;

        F.params = &p;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
        int status0=gsl_integration_qags (&F, 0.0, 1./xx_esd[i], 0, 1.e-3, 1000, w, &result1, &error);
        gsl_integration_workspace_free (w);
        if(status0!=0 && verbose){
            printf("# Problems in convergence, esd1: z:%e r:%e\n",z,xx_esd[i]);
        }

        gsl_integration_workspace * v = gsl_integration_workspace_alloc (1000);
        int status=gsl_integration_qagiu (&F, 1./xx_esd[i], 1.e-3, 0, 1000, v, &result2, &error);
        gsl_set_error_handler(oldhand);
        if(status!=0 && verbose){
            printf("# Problems in convergence, esd2: z:%e r:%e\n",z,xx_esd[i]);
        }
        //stat=gsl_integration_qagiu (&F, 1.e-30, 1.e-3, 0, 1000, v, &result2, &error);
        gsl_integration_workspace_free (v);

        sd[i]=2./1.E12*xx_esd[i]*rho_aver*(result1+result2);
        ww[i]=sd[i]*xx_esd[i];

    }
    if(verbose){
        std::cout<<"# ESD calced\n";
    }

    // Initialize spline for ESD
    sd_acc = gsl_interp_accel_alloc ();
    sd_spline = gsl_spline_alloc (gsl_interp_cspline, rpbins_esd2);
    gsl_spline_init (sd_spline,xx_esd,ww,rpbins_esd2);
    bool_init_sd=true; // This variable is not used anywhere

    for(int i=first;i<rpbins_esd2;i++)
    {
        esd[i-first]=gsl_spline_eval_integ(sd_spline,0.0,xx_esd[i],sd_acc);
        esd[i-first]=2.0/pow(xx_esd[i],2.0)*esd[i-first];
        esd[i-first]=esd[i-first]-sd[i];
        //std::cout<<esd_Rp[i-first]<<" "<<esd[i]<<" "<<wp_Rp[i-first]<<" "<<wp[i-first]<<std::endl;
    }
    //exit(101);

    gsl_spline_free(sd_spline);
    gsl_interp_accel_free(sd_acc);

    delete []xx;
    delete []sd;
    delete []ww;
    delete []xx_esd;

    /// End of story!!!
    hod_free();
    return 0;
}

/// The inputs should be: rpbins_wp, rpbins_esd
/// The output should be: wp
double hod::Wp(double z, int rpbins_wp, double wp_Rp[], double wp[],double pimax,bool reset)
{
    if(reset){
        // Reset the redshift
        resetz(z);
    }


    // We are going to use qagw integration.
    // I need to calculate these at the N projected radius bins from
    // the data. For now just choose some bins for convenience

    //double rho_aver=Omega0*rho_crit_0;

    double *xx;
    xx=new double[rpbins_wp];

    for(int i=0;i<rpbins_wp;i++)
    {
        xx[i]=wp_Rp[i];
        //printf("%e \n",wp_Rp[i]);
    }

    for(int i=0;i<rpbins_wp;i++)
    {
        /// QAGW integral
        ///Int_{0}^{ymax} dy 2  R \xi_gg( \sqrt(y^2+1) R)
        /// ymax=pimax/R
        double result, error;
        double ymax=(pimax/xx[i]);

        gsl_function F;
        F.function = &(dwp);

        proj_params p;
        p.hptr = this;
        p.R = &(xx[i]);
        p.z = &z;

        F.params = &p;
        //std::cout<<i<<" YMAX "<<xx[i]<<" "<<pimax/xx[i]<<std::endl;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
        int status=gsl_integration_qag (&F, 0.0, ymax, 0, 1.e-3, 1000, 5, w, &result, &error);
        gsl_set_error_handler(oldhand);
        if(status!=0 && verbose){
            printf("# Problems in convergence, wp: z:%e r:%e\n",z,xx[i]);
        }
        gsl_integration_workspace_free (w);

        wp[i]=2.0*result*xx[i];
        //std::cout<<xx[i]<<" "<<wp[i]<<std::endl;

    }
    //exit(0);
    if(verbose){
        std::cout<<"# Wp calced\n";
    }

    delete []xx;
    //gsl_set_error_handler();

    /// End of story!!!
    hod_free();
    return 0;
}

/// The inputs should be: rpbins_wp, rpbins_esd, wp_Rp, esd_Rp
/// The output should be: wp and esd
double hod::ESD(double z, int rpbins_esd, double esd_Rp[], double esd[], int rpbins_esd2,bool reset)
{
    if(reset){
        // Reset the redshift
        resetz(z);
    }

    if(verbose) printf("# Bools %d %d:\n",bool_init_D2gg,bool_init_D2gd);

    // We are going to use qagw integration.
    // I need to calculate these at the N projected radius bins from
    // the data. For now just choose some bins for convenience

    double rho_aver=Omega0*rho_crit_0;

    /// Surface density integral
    double *sd,*ww,*xx_esd;
    sd=new double[rpbins_esd2];
    ww=new double[rpbins_esd2];
    xx_esd=new double[rpbins_esd2];

    int first=rpbins_esd2-rpbins_esd;

    for(int i=first;i<rpbins_esd2;i++)
    {
        xx_esd[i]=esd_Rp[i-first];
    }
    double minrp=log10(esd_Rp[0]);
    for(int i=first-1;i>-1;i--)
    {
        xx_esd[i]=pow(10.0,minrp-0.1*(first-i));
    }


    for(int i=0;i<rpbins_esd2;i++)
    {
        //std::cout<<"SD integral"<<xx_esd[i]<<std::endl;
        /// QAGIU integral
        ///Int_{0}^{\infty} dy R \xi_gm( \sqrt(z^2+1) R)
        double result1, result2, error;
        result1=result2=error=0.0;

        gsl_function F;
        F.function = &(dsigma);

        proj_params p;
        p.hptr = this;
        p.R = &(xx_esd[i]);
        p.z = &z;

        F.params = &p;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
        int status0=gsl_integration_qags (&F, 0.0, 1./xx_esd[i], 0, 1.e-3, 1000, w, &result1, &error);
        gsl_integration_workspace_free (w);
        if(status0!=0 && verbose){
            printf("# Problems in convergence, esd: z:%e r:%e\n",z,xx_esd[i]);
        }

        gsl_integration_workspace * v = gsl_integration_workspace_alloc (1000);
        int status=gsl_integration_qagiu (&F, 1./xx_esd[i], 1.e-3, 0, 1000, v, &result2, &error);
        gsl_set_error_handler(oldhand);
        if(status!=0 && verbose){
            printf("# Problems in convergence, esd: z:%e r:%e\n",z,xx_esd[i]);
        }
        //stat=gsl_integration_qagiu (&F, 1.e-30, 1.e-3, 0, 1000, v, &result2, &error);
        gsl_integration_workspace_free (v);

        sd[i]=2./1.E12*xx_esd[i]*rho_aver*(result1+result2);
        ww[i]=sd[i]*xx_esd[i];

        //std::cout<<xx_esd[i]<<" "<<ww[i]<<std::endl;

    }
    if(verbose){
        std::cout<<"# ESD calced\n";
    }
    //exit(101);
    xx_esd[0]=0;
    ww[0]=0;

    // Initialize spline for ESD
    sd_acc = gsl_interp_accel_alloc ();
    sd_spline = gsl_spline_alloc (gsl_interp_cspline, rpbins_esd2);
    gsl_spline_init (sd_spline,xx_esd,ww,rpbins_esd2);
    bool_init_sd=true; // This variable is not used anywhere

    for(int i=first;i<rpbins_esd2;i++)
    {
        //std::cerr<<"Calculating ESD for: "<<esd_Rp[i-first]<<std::endl;
        esd[i-first]=gsl_spline_eval_integ(sd_spline,0.0,xx_esd[i],sd_acc);
        esd[i-first]=2.0/pow(xx_esd[i],2.0)*esd[i-first];
        esd[i-first]=esd[i-first]-sd[i];
        //std::cerr<<esd_Rp[i-first]<<" "<<esd[i]<<std::endl;
    }
    //exit(101);

    gsl_spline_free(sd_spline);
    gsl_interp_accel_free(sd_acc);

    delete []sd;
    delete []ww;
    delete []xx_esd;

    /// End of story!!!
    hod_free();
    return 0;
}

/// The inputs should be: rpbins_wp, rpbins_esd, wp_Rp, esd_Rp
/// The output should be: wp and esd
double hod::Sigma(double z, int rpbins_sd, double sd_Rp[], double sigma[], bool reset, double xfgm_m0, double xfgm_slp, double pimax)
{
    if(reset){
        // Reset the redshift
        resetz(z);
    }

    fgm_m0=xfgm_m0;
    fgm_slp=xfgm_slp;

    if(verbose) printf("# Bools %d %d:\n",bool_init_D2gg,bool_init_D2gd);

    // We are going to use qagw integration.
    // I need to calculate these at the N projected radius bins from
    // the data. For now just choose some bins for convenience

    double rho_aver=Omega0*rho_crit_0;

    /// Surface density integral
    for(int i=0;i<rpbins_sd;i++)
    {
        //std::cout<<"SD integral"<<xx_esd[i]<<std::endl;
        /// QAGIU integral
        ///Int_{0}^{\infty} dy R \xi_gm( \sqrt(z^2+1) R)
        double result1, result2, error;
        result1=result2=error=0.0;

        gsl_function F;
        F.function = &(dsigma);

        proj_params p;
        p.hptr = this;
        p.R = &(sd_Rp[i]);
        p.z = &z;

        F.params = &p;

        gsl_integration_workspace * v = gsl_integration_workspace_alloc (1000);
        gsl_error_handler_t *oldhand=gsl_set_error_handler_off();

        int status;
        if(pimax<0.0){
            status=gsl_integration_qagiu (&F, 1E-3, 1.e-1, 0, 1000, v, &result2, &error);
        }else{
            double ymax=pimax/sd_Rp[i];
            status=gsl_integration_qag (&F, 0.0, ymax, 0, 1.e-3, 1000, 5, v, &result2, &error);
        }
        gsl_set_error_handler(oldhand);
        if(status!=0 && verbose){
            printf("# Problems in convergence, esd: z:%e r:%e\n",z,sd_Rp[i]);
        }
        //stat=gsl_integration_qagiu (&F, 1.e-30, 1.e-3, 0, 1000, v, &result2, &error);
        gsl_integration_workspace_free (v);

        if(pimax<0.0){
            sigma[i]=2./1E12*rho_aver*(result2)*sd_Rp[i];
        }else{
            // This one is DD(rp)/RR(rp) both num and denom integrated along pi upto pimax?
            //sigma[i]=(result2)*sd_Rp[i]/pimax;
            // This one is 2 \int xigm dz
            sigma[i]=2*(result2)*sd_Rp[i];
        }
        //std::cout<<sd_Rp[i]<<" "<<sigma[i]<<std::endl;

    }
    if(verbose){
        std::cout<<"# ESD calced\n";
    }

    /// End of story!!!
    hod_free();
    return 0;
}

/// Free unrequired memory
void hod::hod_free()
{
    if(bool_init_D2gg){
        gsl_interp_accel_free(D2gg_acc);
        gsl_spline_free(D2gg_spline);
        bool_init_D2gg=false;
    }

    if(bool_init_D2gd){
        gsl_interp_accel_free(D2gd_acc);
        gsl_spline_free(D2gd_spline);
        bool_init_D2gd=false;
    }

    if(bool_init_xigg){
        gsl_interp_accel_free(xigg_acc);
        gsl_spline_free(xigg_spline);
        bool_init_xigg=false;
    }

    if(bool_init_xigd){
        gsl_interp_accel_free(xigd_acc);
        gsl_spline_free(xigd_spline);
        bool_init_xigd=false;
    }

    if(bool_init_xigg_bar){
        gsl_interp_accel_free(xigg_bar_acc);
        gsl_spline_free(xigg_bar_spline);
        bool_init_xigg_bar=false;
    }

    if(bool_init_xigg_barbar){
        gsl_interp_accel_free(xigg_barbar_acc);
        gsl_spline_free(xigg_barbar_spline);
        bool_init_xigg_barbar=false;
    }
}

void hod::print(cosmo p, hodpars h)
{
    printf(" %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg \n", p.Om0, p.Omk, p.w0, p.wa, p.s8, p.nspec, p.Omb, p.hval, p.ximax, p.cfac, h.Mmin, h.siglogM, h.Msat, h.alpsat, h.Mcut, h.csbycdm, h.fac );
}

void hod::print(FILE *&fp, cosmo p, hodpars h)
{
    fprintf(fp,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg \n", p.Om0, p.s8, p.nspec, p.Omb, p.hval, p.ximax, p.cfac, h.Mmin, h.siglogM, h.Msat, h.alpsat, h.Mcut, h.csbycdm, h.fac);
}

void hod::print(FILE *&fp, cosmo p, hodpars h, double facc, double c2)
{
    fprintf(fp," %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg \n", p.Om0, p.Omk, p.w0, p.wa, p.s8, p.nspec, p.Omb, p.hval, p.ximax, p.cfac, h.Mmin, h.siglogM, h.Msat, h.alpsat, h.Mcut, h.csbycdm, h.fac, facc, c2);
}

void hod::print()
{
    //printf("# %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg \n", Omega0, Omegal, w0, wa, sigma8, ns, Omegab, h, xiNLzetamax, cfactor, hodp.Mmin, hodp.siglogM, hodp.Msat, hodp.alpsat, hodp.Mcut, hodp.csbycdm, hodp.fac );
    printf("# %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg %.8lg \n", Omega0, Omegal, w0, wa, sigma8, ns, Omegab, h, xiNLzetamax, cfactor, hodp.Mmin, hodp.siglogM, hodp.Msat, hodp.alpsat, hodp.Mcut, hodp.csbycdm, hodp.fac,off_rbyrs,fcen_off,inc_alp,inc_xM);
}

/*
   Kaiser correction part, only for the nonlinear power spectrum
 */

double hod::xigg_bar_num(double r,double z){
    double result=0.0;
    if(!bool_init_xigg_bar){
        init_xigg_bar(z);
    }
    double logr=log10(r);
    result=gsl_spline_eval (xigg_bar_spline,logr, xigg_bar_acc);
    return pow(10.,result)-1.;
}

void hod::init_xigg_bar(double z)
{
    if(verbose){
        std::cout<<"# "<<"============================================================="<<std::endl;
        std::cout<<"# "<<"The spline for numerical calculation of xiggbar was not initialized. Initializing..."<<std::endl;
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
        yy[i]=log10(1.+xigg_bar(r,z));
    }
    xigg_bar_acc = gsl_interp_accel_alloc ();
    xigg_bar_spline = gsl_spline_alloc (gsl_interp_cspline, Nxibar);
    gsl_spline_init (xigg_bar_spline,xx,yy,Nxibar);

    bool_init_xigg_bar=true;
    if(verbose){
        std::cout<<"# "<<"The spline for numerical calculation of xiggbar is now initialized"<<std::endl;
        std::cout<<"# "<<"============================================================="<<std::endl;
    }
}

double dxiggbar(double x, void* params)
{
    ///Int_{0}^{1.0} dy y^2 \xi_gg( y*r )
    xi_gg_params c1 = *(xi_gg_params *)params;

    hod *c2;
    double *r;
    double *z;

    c2=c1.hptr;
    r=c1.r;
    z=c1.z;

    double arg=x*(*r);
    double res=(*c2).xigg_num(arg,*z);
    res=res*x*x;

    return res;
}

double hod::xigg_bar(double r,double z){
    /// QAGW integral
    ///Int_{0}^{1.0} dy y^2 \xi_gg( y*r )

    double result, error;
    gsl_function F;
    F.function = &(dxiggbar);

    xi_gg_params p;
    p.hptr = this;
    p.r = &(r);
    p.z = &z;

    F.params = &p;

    int stat;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
    stat=gsl_integration_qags (&F, 0.0, 1.0, 0, 1.e-3, 1000, w, &result, &error);
    gsl_set_error_handler(oldhand);
    if(stat!=0){
        printf("# Problems in convergence, xigg_bar: z:%e r:%e\n",z,r);
        print();
    }
    gsl_integration_workspace_free (w);

    return result*3.;
}

double hod::xigg_barbar_num(double r,double z){
    double result=0.0;
    if(!bool_init_xigg_barbar){
        init_xigg_barbar(z);
    }
    double logr=log10(r);
    result=gsl_spline_eval (xigg_barbar_spline,logr, xigg_barbar_acc);

    return pow(10.,result)-1.;
}

void hod::init_xigg_barbar(double z)
{
    if(verbose){
        std::cout<<"# "<<"============================================================="<<std::endl;
        std::cout<<"# "<<"The spline for numerical calculation of xiggbarbar was not initialized. Initializing..."<<std::endl;
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
        yy[i]=log10(1.+xigg_barbar(r,z));
    }
    xigg_barbar_acc = gsl_interp_accel_alloc ();
    xigg_barbar_spline = gsl_spline_alloc (gsl_interp_cspline, Nxibar);
    gsl_spline_init (xigg_barbar_spline,xx,yy,Nxibar);

    bool_init_xigg_barbar=true;
    if(verbose){
        std::cout<<"# "<<"The spline for numerical calculation of xiggbarbar is now initialized"<<std::endl;
        std::cout<<"# "<<"============================================================="<<std::endl;
    }
}

double dxiggbarbar(double x, void* params)
{
    ///Int_{0}^{1.0} dy y^4 \xi_gg( y*r )
    xi_gg_params c1 = *(xi_gg_params *)params;

    hod *c2;
    double *r;
    double *z;

    c2=c1.hptr;
    r=c1.r;
    z=c1.z;

    double arg=x*(*r);
    double res=(*c2).xigg_num(arg,*z);
    res=res*pow(x,4.);

    return res;
}

double hod::xigg_barbar(double r,double z){
    /// QAGW integral
    ///Int_{0}^{1.0} dy y^2 \xi_gg( y*r )
    ///Int_{-15.0}^{0.0} d ln y y^3 \xi_gg( y*r )

    double result, error;
    gsl_function F;
    F.function = &(dxiggbarbar);

    xi_gg_params p;
    p.hptr = this;
    p.r = &(r);
    p.z = &z;

    F.params = &p;

    int stat;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    //stat=gsl_integration_qag (&F, -15.0, 0.0, 0, 1.e-4, 1000, 5, w, &result, &error);
    gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
    stat=gsl_integration_qags (&F, 0.0, 1.0, 0, 1.e-3, 1000, w, &result, &error);
    gsl_set_error_handler(oldhand);
    if(stat!=0){
        printf("# Problems in convergence, xigg_barbar: z:%e r:%e\n",z,r);
    }
    gsl_integration_workspace_free (w);

    return result*5.;
}

double dwp_kaiser(double x, void* params)
{
    /// Wp(R) = 2.*Int_{0}^{xmax} xi(rp/sqrt(1-mu*mu)) d mu/(1-mu*mu)^{3/2}
    /// R is multiplied outside the integral
    proj_params c1 = *(proj_params *)params;

    hod *c2;
    double *R;
    double z;

    c2=c1.hptr;
    R=c1.R;
    z=*(c1.z);

    double sinm=sqrt(1-x*x);
    double arg=(*R)/sinm;

    double xi=(*c2).xigg_num(arg,z);
    double xibar=(*c2).xigg_bar_num(arg,z);
    double xibarbar=(*c2).xigg_barbar_num(arg,z);
    //printf("%e %e %e %e\n",arg,xi,xibar,xibarbar);
    double musq=x*x;
    double p2mu=(3.*musq-1.)/2.;
    double p4mu=(35.*musq*musq-30.*musq+3.)/8.;

    /// Kaiser factor will be set up by now!
    double fkai=(*c2).fKaiser;

    double res=(xi*(1.+2./3.*fkai+1./5.*fkai*fkai)+p2mu*(4./3.*fkai+4./7.*fkai*fkai)*(xi-xibar)+p4mu*(8./35*fkai*fkai)*(xi+5./2.*xibar-7./2.*xibarbar));

    return res/pow(sinm,3.0);
}

/// The inputs should be: rpbins_wp, rpbins_esd
/// The output should be: wp
double hod::Wp_Kaiser(double z, int rpbins_wp, double wp_Rp[], double wp[], double pimax,bool reset)
{
    if(reset){
        // Reset the redshift
        resetz(z);
    }

    // We are going to use qagw integration.
    // I need to calculate these at the N projected radius bins from
    // the data. For now just choose some bins for convenience

    double *xx;
    xx=new double[rpbins_wp];

    for(int i=0;i<rpbins_wp;i++)
    {
        xx[i]=wp_Rp[i];
    }

    for(int i=0;i<rpbins_wp;i++)
    {
        /// QAGW integral
        ///Int_{0}^{ymax} dy 2  R \xi_gg( \sqrt(y^2+1) R)
        /// ymax=pimax/R
        double result, error;
        double ymax=(1./sqrt(xx[i]/pimax*xx[i]/pimax+1.0));

        gsl_function F;
        F.function = &(dwp_kaiser);

        proj_params p;
        p.hptr = this;
        p.R = &(xx[i]);
        p.z = &z;

        F.params = &p;

        int stat;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_error_handler_t *oldhand=gsl_set_error_handler_off();
        stat=gsl_integration_qag (&F, 0.0, ymax, 0, 1.e-3, 1000, 5, w, &result, &error);
        gsl_set_error_handler(oldhand);
        if(stat!=0){
            printf("# Problems in convergence, Wp_kaiser: z:%e r:%e\n",z,xx[i]);
        }
        gsl_integration_workspace_free (w);

        wp[i]=2.0*result*xx[i];
    }

    if(verbose){
        std::cout<<"# Wp calced\n";
    }

    delete []xx;

    /// End of story!!!
    hod_free();
    return 0;
}

/// Test the matter power spectrum as a Fourier transform of the simple halo model components
double hod::pspec_halomodel(double z){
    // Reset the redshift
    resetz(z);

    // Set the k-array
    double karr[kbins];
    double xx[kbins];
    double hod_kdiff=(hod_kmax-hod_kmin)/(kbins-1.);
    for(int i=0;i<kbins;i++){
        karr[i]=hod_kmin+i*hod_kdiff;
        xx[i]=karr[i];
        karr[i]=pow(10.,karr[i]);
    }

    // First the various denominators. 
    //  The rho is in comoving units, so no z dependence.
    double rho_aver=Omega0*rho_crit_0;

    double mdep_1h[N9_16],mdep_2h[N9_16];
    double mdep_2h_bcorr=0.0;
    double rvir[N9_16],mvir[N9_16],cvir[N9_16];
    double r200[N9_16],m200[N9_16],c200[N9_16];
    double avmsq=0.0;

    double remainder;
    for(int i=0;i<N9_16;i++)
    {
        double mass=pow(10.,x9_16[i]);
        double nofmdm=nofm(mass,z)*w9_16[i]*mass*log(10.0);
        double bmnmdm=bias(mass,z)*nofm(mass,z)*w9_16[i]*mass*log(10.0);

        m200[i]=x9_16[i];
        modelNFWhalo_com(mass, z, mvir[i], rvir[i], cvir[i], r200[i],c200[i]);
        mdep_1h[i]=pow(mass/rho_aver,2)*nofmdm;
        mdep_2h[i]=(mass/rho_aver)*bmnmdm;
        mdep_2h_bcorr=mdep_2h_bcorr+mass*bmnmdm/rho_aver;
        avmsq=avmsq+mdep_1h[i];
        if(i==N9_16-1){
            remainder=nofm(mass,z)*mass*mass*mass/1.1/pow(rho_aver,2);
            //printf("%e %f\n",mass,remainder);
        }
    }
    mdep_2h_bcorr=1.0-mdep_2h_bcorr;
    std::cout<<"# "<<avmsq+remainder<<std::endl;
    //exit(101);

    /// Initialize Qk here
    if(!bool_initQk){
        initQk(z,r200);
    }
    /// Now the k-dependent portion
    for(int k=0;k<kbins;k++)
    {
        double kfac=karr[k]*karr[k]*karr[k]/(2.0*M_PI*M_PI);
        double int1h=0.0;
        double int2h=0.0;

        double uk_d[N9_16];

        for(int i=0;i<N9_16;i++)
        {
            uk_d[i]=ukinterp(karr[k]*r200[i]/c200[i], c200[i]);
            int1h=int1h+mdep_1h[i]*pow(uk_d[i],2.0);
            int2h=int2h+(mdep_2h[i]*uk_d[i]);
        }
        int2h=int2h+mdep_2h_bcorr;
        double pkl=Delta2_L_num(karr[k],z)/kfac;
        //double pknl=Delta2_NL_num(karr[k],z)/kfac;
        double pknl;
        if(karr[k]<0.01){
            pknl=Delta2_NL_num(karr[k],z)/kfac;
        }else{
            pknl=Pktest_zetaNL(karr[k],z);
        }

        double int2hnl=int2h*int2h*pknl;
        int2h=int2h*int2h*pkl;

        double int2h_he=0.0;
        for(int i=0;i<N9_16;i++)
        {
            for(int j=0;j<N9_16;j++)
            {
                int2h_he+=(mdep_2h[i]*uk_d[i]*mdep_2h[j]*uk_d[j]*Qk1[(k*N9_16+i)*N9_16+j]);
                int2h_he+=(mdep_2h[i]*uk_d[i]*mdep_2h[j]*uk_d[j]*Qk2[(k*N9_16+i)*N9_16+j]);
            }
            int2h_he+=(mdep_2h[i]*uk_d[i]*Qk1[(k*N9_16+i)*N9_16+i]*mdep_2h_bcorr*2);
            int2h_he+=(mdep_2h[i]*uk_d[i]*Qk2[(k*N9_16+i)*N9_16+i]*mdep_2h_bcorr*2);
        }
        int2h_he+=(Qk1[(k*N9_16+0)*N9_16+0]*mdep_2h_bcorr*mdep_2h_bcorr);
        int2h_he+=(Qk2[(k*N9_16+0)*N9_16+0]*mdep_2h_bcorr*mdep_2h_bcorr);

        //if(karr[k]<pow(10.,-1.5)) int1h=0.0;

        std::cout<<karr[k]<<" "<<pknl<<" "<<int1h+int2h<<" "<<int1h+int2h_he<<" "<<int1h<<" "<<pkl<<" "<<int1h+int2hnl<<std::endl;
    }

    return 0;
}

void hod::hod_renew(cosmo p, hodpars h){
    // First free up all the memory associated with all splines
    hod_free();

    // Set up new hod parameters
    hodp=h;

    if(verbose){
        std::cout<<"# New HOD parameters set, setting up the new cosmology now\n";
    }
    // Next call the cosmology renew method
    renew(p);

}

double hod::set_cfactor(double cfac){
    return set_cfac(cfac);
}

double hod::avmass_tot(double z){
    // Get the average mass
    double avMasst_n=0.0;
    double avMass_d=0.0;

    for(int i=0;i<N9_16;i++)
    {
        double m=pow(10.,x9_16[i]);
        double tmp_m=pow(10.,x9_16[i]-12.0);
        avMass_d=avMass_d+log(10.0)*(ncen(x9_16[i])+nsat(x9_16[i]))*m*nofm(m,z)*w9_16[i];
        avMasst_n=avMasst_n+tmp_m*log(10.0)*(ncen(x9_16[i])+nsat(x9_16[i]))*m*nofm(m,z)*w9_16[i];
    }
    return avMasst_n/avMass_d;
}

double hod::avmass_cen(double z){
    // Get the average mass
    double avMassc_n=0.0;
    double avMass_d=0.0;

    for(int i=0;i<N9_16;i++)
    {
        double m=pow(10.,x9_16[i]);
        double tmp_m=pow(10.,x9_16[i]-12.0);
        avMass_d=avMass_d+log(10.0)*ncen(x9_16[i])*m*nofm(m,z)*w9_16[i];
        avMassc_n=avMassc_n+tmp_m*log(10.0)*ncen(x9_16[i])*m*nofm(m,z)*w9_16[i];
    }
    return avMassc_n/avMass_d;
}

double hod::scale_dep_bias_crossr(double z, int rbins, double rr[], double bias[], double crossr[], bool reset){
    if(reset){
        // Reset the redshift
        resetz(z);
    }

    /// Calculating the bias, and cross correlation coefficient
    for(int i=0;i<rbins;i++){
        double xigg=xigg_num(rr[i],z);
        double xigd=xigd_num(rr[i],z);
        double ximm=xi_NL_num(rr[i],z);
        bias[i]=pow(xigg/ximm,0.5);
        crossr[i]=xigd/pow(ximm*xigg,0.5);
    }
    
    return 0;
}
