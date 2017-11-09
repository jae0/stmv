
gaussian_process_stanmodel = function(){

    rstan::stan_model( model_code= "

      functions{
        
        matrix matern_covariance( int N, matrix dist, real phi, real sigma_sq, real tau_sq, real eps, int COVFN) {
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
          for(i in 1:N) S[i,i] = sigma_sq + tau_sq + eps;            
          return(S)   ;
        }

      }

    // ---------------
      data {
        int<lower=1> N;
        // int<lower=0> Np; //number of data locations to predict
        vector[N] Y; // observations
        // vector[N] X; // covariates (predictors)
        matrix[N, N] dist; //distances between points
        real eps;
        // real sigma_sq0;
        // real tau_sq0;
        int<lower=1,upper=4> COVFN;  // Choice of Matern covariance function: 
        // 1:nu=1/2 (exponential); 2:nu=3/2; 3:nu=5/2; 4:nu->Inf (gaussian)
      }

     // ---------------
      transformed data{
        vector[N] zeros;
        zeros = rep_vector(0, N);
      }

     // ---------------
      parameters{
        real<lower=0, upper=10> sigma_sq; 
        real<lower=0, upper=10> tau_sq; 
        real<lower=0, upper=10> phi;
        vector [N] spatialError;
      }
  
     // ---------------
      transformed parameters{
      }
  
     // ---------------
      model{
        matrix[N,N] S; // Covariance
        S = matern_covariance( N, dist, phi, sigma_sq, tau_sq, eps, COVFN );
        Y ~ multi_normal_cholesky( spatialError, cholesky_decompose(S) ) ;
        tau_sq ~ cauchy( 0, 0.5) ;  
        sigma_sq ~ cauchy( 0, 0.5) ;  
        phi ~ normal( 1, 1 ) ; 
      }

     // ---------------
      generated quantities{
      //  vector[N] y_pred;
      //  for(i in 1:N)  y_pred[i] = (beta+plogitobs[i]);
      //  for(i in 1:Np) y_pred[N+i] = (beta+plogitpreds[i]);
      }
    "

  )
}

