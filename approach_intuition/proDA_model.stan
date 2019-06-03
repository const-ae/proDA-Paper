
data {
  int nsamples;
  int nrows;
  int nparams;
  int totalmissing;
  matrix[nsamples, nparams] design_matrix;
  matrix[nrows, nsamples] Y;
  int location_prior_df;
  
  real zeta[nsamples];
  real rho[nsamples];
  real mu0;
  real sigma20;
  real eta;
  real nu;
}
parameters {
  matrix[nrows, nparams] Beta;
  real<lower=0> sigma2[nrows];
}
transformed parameters{
  real zetastar[totalmissing];
  {
    int counter = 1;
    for(i in 1:nrows){
      for(j in 1:nsamples){
        if(is_inf(Y[i, j])){
          zetastar[counter] = zeta[j] * sqrt(1 + sigma2[i]/zeta[j]^2);
          counter = counter + 1;
        }
      }
    }
  }
}
model {
  sigma2 ~ scaled_inv_chi_square(nu, sqrt(eta));
  {
    int counter = 1;
    for(i in 1:nrows){
      for(j in 1:nsamples){
        design_matrix[j, ] * Beta[i, ]' ~ student_t(location_prior_df, mu0, sqrt(sigma20));
        if(is_inf(Y[i, j])){
          target += normal_lccdf(design_matrix[j, ] * Beta[i, ]'  | rho[j], fabs(zetastar[counter]));
          counter += 1;
        }else{
          Y[i, j] ~ normal(design_matrix[j, ] * Beta[i, ]', sqrt(sigma2[i]));
}
      }
    }
  }
}


