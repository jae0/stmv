#ifndef M_PI_2
#define M_PI_2 (6.283185307179518647692) //  2*(3.14159265358979323846)
#endif

#include <R.h>
#include <Rmath.h>
#include <stdlib.h>
#include <fftw3.h>
#include <cmath>
#include <stdint.h> // not allowed to use cstdint as it needs C++11
#include <iostream>
#include <vector>
#include <Ziggurat.h>
// #include <Rcpp.h>  .. not needed due to RcppArmadillo.h included below 
#include <RcppArmadillo.h>


using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

// copied from spate library with some minor optimzations 
// .. moved here to make stm operate more flexibly 

static Ziggurat::Ziggurat::Ziggurat zigg;



// [[Rcpp::export]]
arma::mat inverse_crossproduct (arma::mat A) {
  arma::mat invcp ;
  invcp = (A.t() * A ).i()  ;
  return(invcp) ;
}                                



extern "C" void real_fft_stm(int *n, double yh[], int *inverse, int indCos[], int indW[], int indWCon[], int *NFc){
  fftw_complex *in, *out;
  fftw_plan p;
  int i;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * *n * *n);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * *n * *n);
  if(*inverse==1){
    for(i=0; i<(*n * *n);i++){
      in[i][0]=yh[i];
      in[i][1]=0.0;
    }
    p = fftw_plan_dft_2d(*n,*n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    for(i=0; i<4;i++){
      yh[i]=out[indW[i]-1][0]/ *n;
    }
    for(i=0; i< *NFc;i++){
      yh[indCos[i]-1]=sqrt(2) * out[indW[indCos[i]-1]-1][0]/ *n;
      yh[indCos[i]]=-sqrt(2) * out[indW[indCos[i]]-1][1]/ *n;
    }
  }else{
    for(i=0; i<4;i++){
       in[indW[i]-1][0]=yh[i];
       in[indW[i]-1][1]=0;
    }
    for(i=0; i< *NFc;i++){
      in[indW[indCos[i]-1]-1][0]=yh[indCos[i]-1]/sqrt(2);
      in[indW[indCos[i]-1]-1][1]=-yh[indCos[i]]/sqrt(2);
      in[indWCon[i]-1][0]=yh[indCos[i]-1]/sqrt(2);
      in[indWCon[i]-1][1]=yh[indCos[i]]/sqrt(2);
    }
    p = fftw_plan_dft_2d(*n,*n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    for(i=0; i<(*n * *n);i++){
      yh[i]=out[i][0]/ *n;
    }
  }
  if (NULL != in) fftw_free(in);
  if (NULL != out) fftw_free(out);
  if (NULL != p) fftw_destroy_plan(p);
}



extern "C" void propagate_spectral_stm(double xtp1[], double xt[], double G11C[], double G11[], double G12[], int *NFc, int *ns){
  int j, k;
  for(j=0; j<*ns;j++){
    xtp1[j]=G11C[j] * xt[j];
  }
  for(j=0; j<  *NFc;j++){
    k = 2 * j+*ns; 
    xtp1[k]=G11[j] * xt[k] + G12[j] * xt[k+1];
    xtp1[k+1]=G11[j] * xt[k+1] - G12[j] * xt[k];
  }
}


extern "C" void kf_spectral_stm(double wFT[], double mtt1[], double mtt[], double rtt1[], double rtt[], double specCosOnly[], double G11C[], double specCosSine[], double G11[], double G12[], double *tau2, int *T, int *NFc, int *ns){
  int t,j,k,l,NF;
  NF=(2 * *NFc + *ns);
  double rtt1Temp,rttTemp;
  for(j=0; j<*ns;j++){
    rttTemp=specCosOnly[j];
    rtt1Temp=0;
    rtt[j]=rttTemp;
    for(t=0; t<*T;t++){
      rtt1Temp=rttTemp*(pow(G11C[j],2))+specCosOnly[j];
      rtt1[j + t*NF]=rtt1Temp;
      rttTemp=rtt1Temp*(1-rtt1Temp/(*tau2+rtt1Temp));
      rtt[j + (t+1)*NF]=rttTemp;
    }
  }
  for(j=0; j < *NFc;j++){
    k = *ns + 2*j;
    rttTemp=specCosSine[j];
    rtt1Temp=0;
    rtt[k]=rttTemp;
    rtt[k+1]=rttTemp;
    for(t=0; t<*T;t++){
      rtt1Temp=rttTemp*(pow(G11[j],2)+pow(G12[j],2))+specCosSine[j];
      rtt1[k + t*NF]=rtt1Temp;
      rtt1[k + 1 +t*NF]=rtt1Temp;
      rttTemp=rtt1Temp*(1-rtt1Temp/(*tau2+rtt1Temp));
      rtt[k + (t+1)*NF]=rttTemp; 
      rtt[k + 1 +(t+1)*NF]=rttTemp; 
    }
  }
  for(t=0; t<*T;t++){
    l = t * NF; 
    propagate_spectral_stm(&mtt1[l], &mtt[l], G11C, G11, G12, NFc, ns);
    for(j=0; j<NF;j++){
      k = l +j ;
      mtt[(t+1) * NF+j] = mtt1[k]+(rtt1[k] * (wFT[k]-mtt1[k])/(*tau2+rtt1[k]));
    }
  }
}


extern "C" void bs_spectral_stm(double simAlpha[], double mtt[], double mtt1[], double rtt[], double rtt1[], double spec[], double G11C[], double G11[], double G12[], int *T, int *NFc, int *ns){
  int NF=(2 * *NFc + *ns);
  int j,k,l,t;
  double *AlMinusMtt1 = (double *) malloc(NF * sizeof(*AlMinusMtt1));
  double *Prop = (double *) malloc(NF * sizeof(*Prop));
  double simTemp, mt, RttBar;
  double *G12t = (double *) malloc(NF * sizeof(*G12t));
  for(j=0; j < *NFc; j++){
    G12t[j]=-G12[j];
  }
  for(j=0; j<NF;j++){
    simTemp=mtt[*T * NF +j]+sqrt(rtt[*T * NF +j])*zigg.norm(); // Ziggurat's rnorm(0, 1)
    simAlpha[(*T-1) * NF +j] = simTemp;
    AlMinusMtt1[j]=simTemp-mtt1[(*T-1) * NF +j];
  }
  for(t=(*T-1); t>0;t--){
    propagate_spectral_stm(Prop, AlMinusMtt1, G11C, G11, G12t, NFc, ns);
    for(j=0; j<NF;j++){
      k = t * NF +j ;
      l = (t-1) * NF +j ;
      mt = mtt[k]+rtt[k]/rtt1[k]*Prop[j];
      RttBar = rtt[k] * (1-(rtt1[k]-spec[j])/rtt1[k]);
      simTemp=mt+sqrt(RttBar)*zigg.norm();  
      simAlpha[l] = simTemp;
      AlMinusMtt1[j]=simTemp-mtt1[l];
    }
  }
}



extern "C" void ll_spectral_stm(double *ll, double wFT[], double mtt1[], double rtt1[],  int *T, int *NF, double *tau2){
  *ll=0;
  int t,j,k,l;
  for(t=0; t<*T;t++){
    l = t * *NF;
    for(j=0; j<*NF; j++){
      k = l + j ;
      *ll = *ll - log(*tau2+rtt1[k]) - (wFT[k]-mtt1[k])*(wFT[k]-mtt1[k])/(*tau2+rtt1[k]);
    }
  }
  *ll = *ll/2-*T * *NF * log(M_PI_2)/2;
}





// [[Rcpp::export]]
extern "C" void TSreal_fft_stm(int *n, int *T, double yh[], int *inverse, int indCos[], int indW[], int indWCon[], int *NFc){
  int t;
  for(t=0; t<*T;t++){
    real_fft_stm(n, &yh[t * *n * *n], inverse, indCos, indW, indWCon, NFc);
  }
}


extern "C" void ffbs_spectral_stm(double wFT[], double *bw, double *ll, double specCosOnly[], double G11C[], double specCosSine[], double G11[], double G12[], double specAll[], double *tau2, int *T, int *NFc, int *ns){
/* , double simAlpha[] */
  int NF=(2 * *NFc + *ns);
  int j;
  double *rtt = (double *) malloc(NF * (*T +1)*sizeof(*rtt));
  double *rtt1 = (double *) malloc(NF * (*T)*sizeof(*rtt1));
  double *mtt = (double *) malloc(NF * (*T+1)*sizeof(*mtt));
  double *mtt1 = (double *) malloc(NF * (*T)*sizeof(*mtt1));

  for(j=0; j<NF;j++){
    mtt[j]=0.0;
  }
  
  kf_spectral_stm(wFT, mtt1, mtt, rtt1, rtt, specCosOnly, G11C, specCosSine, G11, G12, tau2, T, NFc, ns);

  if(*ll==1){
    ll_spectral_stm(ll, wFT, mtt1, rtt1, T, &NF, tau2);
   // printf("%f",*ll); 
  }
  if(*bw==1){
    bs_spectral_stm(wFT, mtt, mtt1, rtt, rtt1, specAll, G11C, G11, G12, T, NFc, ns);
  }

  free(rtt);
  free(rtt1);
  free(mtt);
  free(mtt1);
  rtt=NULL;
  rtt1=NULL;
  mtt=NULL;
  mtt1=NULL;
}


extern "C" void ffbs_spectral_stm_oneshot(  
    double yh[], 
    double *ll, 
    double *tau2,
    int indCosOnly[],
    int indCos[],
    int indW[],
    int indWCon[],
    double spec[],
    double DiffDamp[],
    double Adv[],
    int *T, 
    int *n, 
    int *NF, 
    int *NFc, 
    int *ns ) {
  
  int i, j, k, l, m, q, t;

  int nn;
  double sqrt2, invsqrt2, sqrt2_n;
  
  nn = *n * *n;
  sqrt2 = sqrt(2);
  invsqrt2 = 1/sqrt(2);
  sqrt2_n = sqrt(2) / *n;
 
  fftw_complex *in, *out;
  fftw_plan p;
  
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nn);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nn);

  double *rtt = (double *) malloc(*NF * (*T +1) * sizeof(*rtt));
  double *rtt1 = (double *) malloc(*NF * (*T) * sizeof(*rtt1));
  double *mtt = (double *) malloc(*NF * (*T+1) * sizeof(*mtt));
  double *mtt1 = (double *) malloc(*NF * (*T) * sizeof(*mtt1));
  double *wFT = (double *) malloc(*NF * (*T) * sizeof(*wFT)) ;

  double *G11C = (double *) malloc( *ns * sizeof(*G11C)) ;
  double *G11  = (double *) malloc( *NFc * sizeof(*G11)) ;
  double *G12  = (double *) malloc( *NFc * sizeof(*G12)) ;

  double *specCosOnly  = (double *) malloc( *ns * sizeof(*specCosOnly)) ;
  double *specCosSine  = (double *) malloc( *NFc * sizeof(*specCosSine)) ;

  int *iC  = (int *) malloc( *ns * sizeof(*iC)) ;
  int *iCos  = (int *) malloc( *NFc * sizeof(*iCos)) ;
  int *iWCon  = (int *) malloc( *NFc * sizeof(*iWCon)) ;
  int *iW  = (int *) malloc( *NF * sizeof(*iW)) ;

  for (i=0; i < *ns; i++) {
    iC[i] = indCosOnly[i]-1 ; // -1 for C indexing
    j = iC[i];
    specCosOnly[i] = spec[j]; 
    G11C[i] = exp( DiffDamp[j] ) ;    
  }

  for (i=0; i < *NFc; i++) {
    iWCon[i] = indWCon[i] - 1;
    iCos[i]  = indCos[i]-1 ; // -1 for C indexing
    j = iCos[i];
    specCosSine[i] = spec[j]; 
    G11[i]  =  exp( DiffDamp[j] ) * cos( Adv[j] ) ;   
    G12[i]  = -exp( DiffDamp[j] ) * sin( Adv[j] ) ;
  }

  for(i=0; i < *NF; i++){
    mtt[i]=0.0;
    iW[i] = indW[i] - 1;
  }


  // convert to spectral space
  for(t=0; t<*T; t++) {
    j = t * nn ;

    for(i=0; i<nn; i++){
      k = j +i;
      in[i][0]=yh[k];
      in[i][1]=0.0;
    }
    p = fftw_plan_dft_2d(*n, *n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    for(i=0; i<4; i++){
      k = j +i;
      wFT[k]=out[iW[i]][0]/ *n;
    }
    for(i=0; i< *NFc; i++){
      l = iCos[i] ;
      m = j + l;
      wFT[m]  = sqrt2_n * out[iW[l]][0] ;
      wFT[m+1]=-sqrt2_n * out[iW[l+1]][1] ;
    }
  }

  kf_spectral_stm(wFT, mtt1, mtt, rtt1, rtt, specCosOnly, G11C, specCosSine, G11, G12, tau2, T, NFc, ns);
  ll_spectral_stm(ll, wFT, mtt1, rtt1, T, NF, tau2);
  bs_spectral_stm(wFT, mtt, mtt1, rtt, rtt1, spec, G11C, G11, G12, T, NFc, ns);

  
  // return to normal space
  for(t=0; t<*T; t++){
    j = t * nn ;
    for(i=0; i<4;i++){
       k = iW[i]; 
       l = j+i; 
       in[k][0]=wFT[l];
       in[k][1]=0;
    }
    for(i=0; i< *NFc; i++){
      k = iW[iCos[i]];
      l = iWCon[i];
      m = j+k;
      in[k][0] = wFT[m]   * invsqrt2;
      in[k][1] =-wFT[m+1] * invsqrt2;
      in[l][0] = wFT[m]   * invsqrt2;
      in[l][1] = wFT[m+1] * invsqrt2;
    }
    p = fftw_plan_dft_2d(*n,*n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    for(i=0; i<nn; i++){
      k = j+i; 
      yh[k] = out[i][0] / *n;
    }
  }


  if (NULL != G11C) free(G11C);
  if (NULL != G11) free(G11);
  if (NULL != G12) free(G12);
  if (NULL != wFT) free(wFT);
  if (NULL != rtt) free(rtt);
  if (NULL != rtt1) free(rtt1);
  if (NULL != mtt) free(mtt);
  if (NULL != mtt1) free(mtt1);
  if (NULL != specCosOnly) free(specCosOnly);
  if (NULL != specCosSine) free(specCosSine);

  if (NULL != in) fftw_free(in);
  if (NULL != out) fftw_free(out);
  if (NULL != p) fftw_destroy_plan(p);
  

}


/*

Testing:

  git commit -am"debug" && R 

  install_git( "/home/jae/bio/stm", branch="clibs", force=TRUE) 
  
    require(stm)
    n = 12
    nn=n*n
    T=10
    
    parh_proposal = c(rho0 = 0.25, 
        sigma2 = 0.2, zeta = 0.25, rho1 = 0.25, gamma = 1, alpha = 0.1, muX = 0, muY = 0, tau2 = 0.2 )
    
    y = spate::spate.sim(par=parh_proposal, n=n, T=T)$xi

    NF=n*n
    Z = spate::spate.init(n=n, T=T, NF=NF)
    ns = Z$ns
    NFc = (NF-Z$ns)/2
    indCosOnly = 1:4
 
    dt = 1
    nu = 1
    w2 = apply(Z$wave^2, 2, sum)
    d = 2

    Tr = cbind( c(cos(parh_proposal["alpha"]), - parh_proposal["gamma"] * sin(parh_proposal["alpha"])),
                c(sin(parh_proposal["alpha"]),   parh_proposal["gamma"] * cos(parh_proposal["alpha"]))) / parh_proposal["rho1"]
    Sig = solve(t(Tr) %*% Tr)

    # diffusion and damping
    DiffDamp0 = - { colSums(Z$wave * crossprod(Sig, Z$wave) ) + parh_proposal["zeta"] } 
    DiffDamp = dt * DiffDamp0
    Adv = c( dt *  crossprod( c(parh_proposal["muX"], parh_proposal["muY"]), Z$wave ) ) # advection

    # spectrum: matern spatial and temporal ar1 
    rho0_inv = 1/parh_proposal["rho0"]
        indCosOnly = 1:Z$ns
    

    #precompute a few spectral constants
    pid2 = pi^(d/2)
    nunu2 = 2^(nu - 1) * nu 
    nu2 = 2 * nu
    nud2 = nu + d/2

    spec = { nunu2 * rho0_inv^nu2 } / { pid2 * {rho0_inv^2 + w2}^nud2 }
    spec[1:Z$ns] = spec[1:Z$ns]/2
    spec = spec * nn / sum(spec) * {1 - exp(2 * DiffDamp)}/{-2}/DiffDamp0 # normalize and include diffusion/advection

microbenchmark::microbenchmark( 

  {  ffbs_proposal = .C( "ffbs_spectral_stm_oneshot", 
        yh=as.double(spate::TSmat.to.vect(y)), 
        ll=as.double(1),
        tau2=as.double(parh_proposal["tau2"]),
        indCosOnly=as.integer(indCosOnly),
        indCos=as.integer(Z$indFFT$indCos),
        indW=as.integer(Z$indFFT$indW),
        indWCon=as.integer(Z$indFFT$indWCon),
        spec=as.double(spec),
        DiffDamp=as.double(DiffDamp),
        Adv=as.double(Adv),
        T=as.integer(T), 
        n=as.integer(n), 
        NF=as.integer(NF), 
        NFc=as.integer(NFc), 
        ns=as.integer(ns),
        PACKAGE="stm" )
}, 
{
    u = .C("TSreal_fft_stm", n=as.integer(n), T= as.integer(T), yh=as.double(spate::TSmat.to.vect(y)), inverse=as.integer(1), 
                  indCos=as.integer(Z$indFFT$indCos), indW=as.integer(Z$indFFT$indW), indWCon=as.integer(Z$indFFT$indWCon), 
                  NFc=as.integer(NFc), PACKAGE="stm")$yh

        G11C = exp(DiffDamp[indCosOnly])    
        G11 =  exp(DiffDamp[Z$indFFT$indCos]) * cos(Adv[Z$indFFT$indCos])   
        G12 = -exp(DiffDamp[Z$indFFT$indCos]) * sin(Adv[Z$indFFT$indCos])

        mll <- .C("ffbs_spectral_stm", wFT=as.double(u), bw=as.double(FALSE), ll=as.double(TRUE), 
                  specCosOnly=as.double(spec[indCosOnly]), G11C=as.double(G11C), specCosSine= as.double(spec[Z$indFFT$indCos]), 
                  G11=as.double(G11), G12=as.double(G12), specAll=as.double(spec), tau2=as.double(parh_proposal["tau2"]), 
                  T=as.integer(T), NFc=as.integer(NFc), ns=as.integer(Z$ns), PACKAGE="stm" )$ll
        mll

         v = .C("ffbs_spectral_stm", wFT=as.double(u), bw=as.double(TRUE), ll=as.double(FALSE), 
              specCosOnly=as.double(spec[indCosOnly]), G11C=as.double(G11C), specCosSine= as.double(spec[Z$indFFT$indCos]), 
              G11=as.double(G11), G12=as.double(G12), specAll=as.double(spec), tau2=as.double(parh_proposal["tau2"]), 
              T=as.integer(T), NFc=as.integer(NFc), ns=as.integer(Z$ns), PACKAGE="stm" )$wFT

          yhat = spate::vect.to.TSmat( .C("TSreal_fft_stm", n=as.integer(n), T=as.integer(T), yh=as.double(v), inverse=as.integer(0), 
              indCos=as.integer(Z$indFFT$indCos), indW=as.integer(Z$indFFT$indW), indWCon=as.integer(Z$indFFT$indWCon),
              NFc=as.integer(NFc), PACKAGE="stm" )$yh, T=T)
}
, times=1000)



   #   str(ffbs_proposal )
    summary(c(ffbs_proposal$yh))
 summary(c(yhat))


    plot(c(y) ~ c(yhat) )
    
    plot(c(y) ~ c(ffbs_proposal$yh) )


 */
