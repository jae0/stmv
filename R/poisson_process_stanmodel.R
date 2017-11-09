
poisson_process_stanmodel = function(){

    rstan::stan_model( model_code= "
     functions {
            matrix matern_covariance( int N, matrix dist, real phi, real sigma_sq, int COVFN) {
              matrix[N,N] S;
              real dist_phi; 
              real sqrt3;
              real sqrt5;
              sqrt3=sqrt(3.0);
              sqrt5=sqrt(5.0);

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
              for(i in 1:N) S[i,i] = sigma_sq ;            
              return(S)   ;
            }

      }

      data {
        int<lower=1> N;
        int<lower=1> Ncov; // number of predictors/covariates including the constant
        int<lower=0> Np; //number of data locations to predict
        int Y[N] ; // observations
        vector[N] logOffset; // expected
        real lambda_mean ; // mean intensity .. poisson: mean==var 
        // real logOffset; // expected
        // matrix[N,Ncov] X; // covariates (predictors)
        matrix[N,N] dist; //distances between points
        int<lower=1,upper=4> COVFN;  // Choice of Matern covariance function: 
        // 1:nu=1/2 (exponential); 2:nu=3/2; 3:nu=5/2; 4:nu->Inf (gaussian)
      }

      transformed data{
        vector[N] zeros;
        zeros = rep_vector(0, N);
      }

      parameters{
        real<lower=0> sigma_sq;  // spatial variance
        real<lower=0> tau_sq; // nugget, nonspatial varaince
        real<lower=0> phi;    // range scale
        // vector[Ncov] beta ; // coeff of covariates
        //--vector[N] eta; // faster way of generating multi_normal_cholesky .. spatialError
        vector[N] spatialError;  // realization of the GP(0,Cov)
      }

      transformed parameters{
        vector[N] log_lambda;
        log_lambda = logOffset + spatialError;
        // log_lambda = logOffset + X*beta + spatialError;
        {
          // keep the following local so that they are not saved
          //--vector[N] spatialError;  // realization of the GP(0,Cov)
          //-- spatialError = L * eta ; // this is a faster way of doing the multi_normal_cholesky, commented out below
        }
      }

      model{
        matrix[N,N] K;
        K = matern_covariance( N, dist, phi, sigma_sq, COVFN);
        // create diagonal: nugget(nonspatial) + spatial variance, also: +  eps ensures positive definiteness
        for(i in 1:N) K[i,i] = K[i,i] + tau_sq + 1e-12 ;            
        spatialError ~ multi_normal_cholesky( zeros, cholesky_decompose(K )) ;  // w(s) ~MVN(0, K); 
        Y ~ poisson_log(log_lambda);
        tau_sq ~ normal( 0, 1) ;  // vague with fat tails 
        sigma_sq ~ normal( 0, 1)  ;
        phi ~ normal( 0, 1) ;  // distance scaled to 1 so phi should be <1 
      }

      generated quantities{
      //  vector[N] y_pred;
      //  for(i in 1:N)  y_pred[i] = (beta+plogitobs[i]);
      //  for(i in 1:Np) y_pred[N+i] = (beta+plogitpreds[i]);
        vector[N] lambda;
        lambda = exp(log_lambda);
      }
   "
  )
}
