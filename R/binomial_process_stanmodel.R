
binomial_process_stanmodel = function(){

    rstan::stan_model( model_code= "

       functions{
     
        matrix matern_covariance( int N, matrix dist, real phi, real sigma_sq, real tau_sq, int COVFN) {
          matrix[N,N] S;
          real dist_phi; 
          real sqrt3;
          real sqrt5;
          real eps;
          sqrt3=sqrt(3.0);
          sqrt5=sqrt(5.0);
          eps = 1e-12;
          
          if (COVFN==1) { // exponential == Matern nu=1/2 , (p=0; nu=p+1/2)
            for(i in 1:(N-1)){
              for(j in (i+1):N){
                dist_phi = fabs(dist[i,j])/phi;
                S[i,j] = sigma_sq * exp(- dist_phi ); 
            }}

          
          } else if (COVFN==2) { // Matern nu= 3/2 covariance
            for(i in 1:(N-1)){
              for(j in (i+1):N){
               dist_phi = fabs(dist[i,j])/phi;
               S[i,j] = sigma_sq * (1 + sqrt3 * dist_phi) * exp(-sqrt3 * dist_phi);
            }}
          

          } else if (COVFN==3) { // Matern nu=5/2 covariance
            for(i in 1:(N-1)){
              for(j in (i+1):N){
                dist_phi = fabs(dist[i,j])/phi;
                S[i,j] = sigma_sq * (1 + sqrt5 *dist_phi + 5* pow(dist_phi,2)/3) * exp(-sqrt5 *dist_phi);
            }}
        
          } else if (COVFN==4) { // Matern as nu->Inf become Gaussian (aka squared exponential cov)
            for(i in 1:(N-1)){
              for(j in (i+1):N){
                dist_phi = fabs(dist[i,j])/phi;
                S[i,j] = sigma_sq * exp( -pow(dist_phi,2)/2 ) ;
            }}
          } 

             
          for(i in 1:(N-1)){
          for(j in (i+1):N){
            S[j,i] = S[i,j];  // fill upper triangle
          }}

          // create diagonal: nugget(nonspatial) + spatial variance +  eps ensures positive definiteness
          for(i in 1:N) S[i,i] = sigma_sq + tau_sq + eps;            
          return(S)   ;
        }

      data {
       int<lower=1> ndata; //number of data points
       int<lower=0> npred; // number of new points
       int ncases[ndata]; // number of outcomes
       int ntot[ndata]; // number of observations
       int nexpected[ndata]; // number of expected
       matrix[ndata+npred, ndata+npred] dist; //distances between points
       int<lower=1,upper=4> COVFN;  // Choice of covariance function: 1=exp; 2=Matern32; 3=Matern5/2
      }
      transformed data{
       matrix[ndata+npredict, ndata+npredict] dist_sq;
       int<lower=1> N;
       real sqrt3;
       real sqrt5;
       // N = ndata + npred;
       N = ndata;
       sqrt3=sqrt(3.0);
       sqrt5=sqrt(5.0);
      }
      parameters{
       vector[ndata] plogitobs;
       // vector[npred] plogitpreds;
        real<lower=0> beta;
        real<lower=0> sigma_sq; 
        real<lower=0> tau_sq;
        real<lower=0> phi;
      }
      transformed parameters{
       // vector[ndata+npred] mu;
       vector[ndata] mu;
       real dist_phi; 
       for(i in 1:N) mu[i] = beta;
      }
      model{
       vector[N] plogit;
       matrix[N,N] K;  // Sigma and Cholesky of Sigma
       K = matern_covariance( N, dist, phi, sigma_sq, tau_sq, COVFN) ; 
       plogit ~ multi_normal_cholesky(mu,  cholesky_decompose(K) );
       tau_sq ~ cauchy(0, 2.5) ;  // half-cauchy
       sigma_sq ~ cauchy(0, 2.5) ;  // half-cauchy
       phi ~ cauchy(0, 2.5) ; // half-cauchy
       beta ~ cauchy(0, 2.5);
       ncases ~ binomial_logit(ntot, plogitobs);
      }
      generated quantities{
    //  vector[npred] y_pred;
    //  for(i in 1:npred) y_pred[i] = inv_logit(beta+plogitpreds[i]);
        vector[ndata] y_pred;
        for(i in 1:ndata) y_pred[i] = inv_logit(beta+plogitobs[i]);
      }
    "

  )
}
