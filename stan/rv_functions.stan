real surv_lpmf(array[] int y, matrix Xc_evdeath, vector mu_evdeath, vector b_evdeath, vector mu_vol, real b_vol, int n_subject_time){
  real lpmf = 0;
  vector[n_subject_time] lp_n;
  // array[n_subject_time] int Y_evdeath;
  for (subject_time in 1:n_subject_time){
    int n_region = 72;
    int i_start = (subject_time - 1)*n_region + 1;
    int i_stop = subject_time*n_region;
    // Y_evdeath[subject_time] = y[i_start];
    // if (y[i_start] != mean( y[i_start:i_stop])) reject("BUG found!");
    lp_n[subject_time] = sum(mu_evdeath[i_start:i_stop] +
      Xc_evdeath[i_start:i_stop,:] * b_evdeath + 
      mu_vol[i_start:i_stop] * b_vol)/n_region;
    // lp_n[subject_time] += b_vol * sum(mu_vol[i_start:i_stop]) / 72;
  }
  // lpmf += bernoulli_logit_lpmf(Y_evdeath | lp_n);
  lpmf += bernoulli_logit_lupmf(y | lp_n);
  return lpmf;
} 


real surv_glm2_lpmf(array[] int y, matrix Xc_evdeath, vector mu_evdeath, vector b_evdeath, vector mu_vol, real b_vol, int n_subject_time){
  real lpmf;
  int n_region = 72;
  matrix[n_subject_time, num_elements(b_evdeath) + 1] X;
  vector[num_elements(b_evdeath) + 1] b = append_row(b_evdeath, b_vol);
  // perform col sum, divided by scale_one to make an average
  row_vector[n_region] v_ones = rep_row_vector(1.0/n_region, n_region);
  //intecept of surv
  vector[n_subject_time] mu;

  for (subject_time in 1:n_subject_time){
    
    int i_start = (subject_time - 1)*n_region + 1;
    int i_stop = subject_time*n_region;
    
    mu[subject_time] = sum(mu_evdeath[i_start:i_stop]);
    // average the matrix of X
    X[subject_time, :]  = v_ones * append_col(Xc_evdeath[i_start:i_stop,:], mu_vol[i_start:i_stop]);
  }
  mu /= n_region; 
  lpmf = bernoulli_logit_glm_lupmf(y | X, mu, b);
  return lpmf;
} 

real surv_glm_lpmf(array[] int y, matrix Xc_evdeath, vector mu_evdeath, vector b_evdeath){
  real lpmf;
  int n_subject_time = num_elements(y);
  int n_region = 72;
  //intecept of surv
  vector[n_subject_time] mu;
  matrix[n_subject_time, num_elements(b_evdeath)] X;

  for (subject_time in 1:n_subject_time){
    
    int i_start = (subject_time - 1)*n_region + 1;
    int i_stop = subject_time*n_region;
    
    mu[subject_time] = sum(mu_evdeath[i_start:i_stop]);
    X[subject_time, :]  = Xc_evdeath[i_start,:];
  }
  mu /= n_region; 
  lpmf = bernoulli_logit_glm_lupmf(y | X, mu, b_evdeath);
  return lpmf;
} 

real partial_surv_glm_lpmf(array[] int y, int start, int end, matrix Xc_evdeath, vector mu_evdeath, vector b_evdeath){
  int n_region = 72;
  int Start = (start - 1)*n_region + 1;
  int End = end*n_region;
  matrix[End-Start+1, dims(Xc_evdeath)[2]] Xc_evdeath2 = Xc_evdeath[Start:End,:];
  vector[End-Start+1] mu_evdeath2 = mu_evdeath[Start:End];
  return surv_glm_lpmf(y | Xc_evdeath2, mu_evdeath2, b_evdeath);
}

real partial_t_lpdf(array[] real y, int start, int end, vector nu, vector mu, vector sigma){
  real tot = student_t_lupdf(y | nu[start:end], mu[start:end], sigma[start:end]);
  
  return tot;
}

real partial_normal_lpdf(array[] real y, int start, int end, vector mu, vector sigma){
  real tot = normal_lupdf(y | mu[start:end], sigma[start:end]);
  
  return tot;
}

real partial_sgt(array[] real y, int start, int end, vector mu, vector sigma, real lambdap1half, real p, real q){
  real tot = 0;
  int k = 1;
  for (n in start:end){
    tot += sgt_lpdf(y[k] | mu[n], sigma[n], lambdap1half, p, q);
    k +=1;
  }
  return tot;
}

real partial_sym_gt(array[] real y, int start, int end, vector mu, vector sigma, real p, real q){
  real tot = 0;
  int k = 1;
  for (n in start:end){
    tot += sym_gt_lpdf(y[k] | mu[n], sigma[n], p, q);
    k +=1;
  }
  return tot;
}

real partial_skewt(array[] real y, int start, int end, vector mu, vector sigma, real lambdap1half, real qmhalf){
  real tot = 0;
  int k = 1;
  for (n in start:end){
    tot += constrained_skew_t_lpdf(y[k] | mu[n], sigma[n], lambdap1half, qmhalf);
    k +=1;
  }
  return tot;
}
