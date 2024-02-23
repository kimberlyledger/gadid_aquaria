// This is a model file for calculating the qPCR results for Skagit eDNA sampling in 2017.

data { /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Number of observations in various categories
    int N_site ;   // Number of Sites
    int N_month ;  // Number of months observed
    int N_site_month; // Number of site-month combinations observed.
    int N_bottle ; // Number of individual bottles observed.
    int N_pcr ;    // Number of PCR plates
    
    int N_bin_stand ;   // Number of observations for binomial part of the standards model
    int N_count_stand ; // Number of observations for count part of the standards model
    int N_bin_samp   ;  // Number of observations for binomial part of the sample model
    int N_count_samp ;  // Number of observations for count part of the sample model

    // Observations
    int bin_stand[N_bin_stand]     ;
    vector[N_count_stand] count_stand ;
    int bin_samp[N_bin_samp]      ; 
    vector[N_count_samp] count_samp   ;
    
    // Covariates
    vector[N_bin_stand] D_bin_stand     ;
    vector[N_count_stand] D_count_stand ;

    // Standard Indices
    int pcr_stand_bin_idx[N_bin_stand] ;
    int pcr_stand_count_idx[N_count_stand] ;
    int pcr_samp_bin_idx[N_bin_samp] ;
    int pcr_samp_count_idx[N_count_samp] ;

    // Site-month and bottle indices
    int site_month_idx[N_bottle]    ;
    int bottle_idx[N_bottle];
    int gamma_idx[N_site_month] ;
    
    // Sample related indices
    int site_bin_idx[N_bin_samp]      ;
    int site_count_idx[N_count_samp]  ;
    int site_month_bin_idx[N_bin_samp]   ;
    int site_month_count_idx[N_count_samp] ;
    int month_bin_idx[N_bin_samp]     ;
    int month_count_idx[N_count_samp] ;
    int bottle_bin_idx[N_bin_samp]    ;
    int bottle_count_idx[N_count_samp];
    
    // counter for number of sites estimated for each month  
    //vector[N_month] counter ;
    
    //Offset
    real OFFSET ;
}
transformed data{
}
parameters { /////////////////////////////////////////////////////////////////////////////////////////////
    // Standards regression and logit coeffs

      real beta_0[N_pcr] ;
      real beta_1[N_pcr] ;

      real phi_0[N_pcr] ;
      real phi_1[N_pcr];



    // Variance parameters for observation and random effects
      real sigma_stand_int ;

      real<lower=0> tau_bottle ;
      real<lower=0> sigma_pcr ; 
    
    // Effects of sites and bottles
      real gamma[N_site_month] ;
      real delta[N_bottle] ;
}
transformed parameters { ////////////////////////////////////////////////////////////////////////////////
    // Latent variables for each log Density 

      vector[N_bottle] D ;
      vector[N_bin_stand] theta_stand ;
      vector[N_bin_samp] theta_samp ;
      vector[N_count_stand] kappa_stand ;
      vector[N_count_samp] kappa_samp ;

      real sigma_all_stand ;
      real sigma_all_samp ;
      
    // Latent variables for each site
    for(i in 1:N_bottle){
      D[i] = gamma[site_month_idx[i]] + delta[bottle_idx[i]] ;
    }
    
    // Presence-Absence component of model.
    for(i in 1:N_bin_stand){
       theta_stand[i] = phi_0[pcr_stand_bin_idx[i]] + phi_1[pcr_stand_bin_idx[i]] *  (D_bin_stand[i] - OFFSET) ;

    }
    for(i in 1:N_bin_samp){
       theta_samp[i]  = phi_0[pcr_samp_bin_idx[i]] + phi_1[pcr_samp_bin_idx[i]] * (D[bottle_bin_idx[i]] - OFFSET) ;
    }
    
    // Positive Comonent of the model
    for(i in 1:N_count_stand){
       kappa_stand[i] = beta_0[pcr_stand_count_idx[i]] + beta_1[pcr_stand_count_idx[i]] * (D_count_stand[i] - OFFSET) ;
    }

    for(i in 1:N_count_samp){
      kappa_samp[i]  = beta_0[pcr_samp_count_idx[i]] + beta_1[pcr_samp_count_idx[i]] * (D[bottle_count_idx[i]] - OFFSET) ;
    }
    
    // 
      sigma_all_stand = sigma_stand_int;
      sigma_all_samp = pow(sigma_stand_int^2 +
                            sigma_pcr^2,-2) ;

}
model {////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Likelihood components
    bin_stand  ~ bernoulli( inv_logit(theta_stand) ) ;
    bin_samp   ~ bernoulli( inv_logit(theta_samp) ) ;
    count_stand ~ normal(kappa_stand, sigma_all_stand) ;
    count_samp ~ normal(kappa_samp, sigma_all_samp) ;

    // Random effects
    gamma ~ normal(-4,8);
    delta ~ normal(0,tau_bottle) ;
    
    // Priors
    sigma_stand_int ~ gamma(1,1) ;

    tau_bottle ~ gamma(1,1) ;
    sigma_pcr ~ gamma(1,1) ;

    beta_0 ~ normal(35,10) ;
    beta_1 ~ normal(-5,5) ;
    
    phi_0 ~ normal(0 , 10) ;
    phi_1 ~ normal(5, 5) ;
    
}
generated quantities{
}
