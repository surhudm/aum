#include "hod.h"
#include <iostream>

int main(){

    cosmo p;
    hodpars q;

    // Mmin:13.15012 siglogM:0.50470 Msat:13.98671 alpsat:1.00000 Mcut:1.39500 Mstel:5.93047 cfactor:0.70408 Omegam:0.29330 sigma8:0.83242

    p.Om0  = 0.2930;
    p.Omk  = 0.0;
    p.w0=-1.;
    p.wa=0.0;
    p.hval = 0.697;
    p.Omb  = 0.0461;
    p.th   = 1.0093*2.7;
    p.s8   = 0.83242;
    p.nspec= 0.9646;
    p.ximax=log10(8.0);
    p.cfac=0.70408;

    q.Mmin   =13.15012;
    q.siglogM=0.50470;
    q.Msat   =13.98671;
    q.alpsat =1.;
    q.Mcut   =log10(1.39500)+q.Mmin;
    q.csbycdm=1;
    q.fac    =-0.25;

    hod h(p,q);

    //h.Pk_gg_gd_he(0.55);
    double *rarr;
    double *wparr;
    int Nbins=12;
    rarr=new double [Nbins];
    wparr=new double [Nbins];
    for(int i=0;i<Nbins;i++){
	rarr[i]=pow(10.,-2.0+i*3./(Nbins-1.));
    }
    double zcalc=0.55;
    double pimax=200.0;
    h.Wp_Kaiser(zcalc,Nbins,rarr,wparr,pimax);
    for(int i=0;i<Nbins;i++){
	rarr[i]=pow(10.,-2.0+i*3./(Nbins-1.));
	fprintf(stderr,"%.4f %.4f\n",rarr[i],wparr[i]);
    }

    fprintf(stderr,"=======================\n");
    h.ESD(zcalc,Nbins,rarr,wparr,Nbins+6);
    for(int i=0;i<Nbins;i++){
	rarr[i]=pow(10.,-2.0+i*3./(Nbins-1.));
	fprintf(stderr,"%.4f %.4f\n",rarr[i],wparr[i]);
    }

}
