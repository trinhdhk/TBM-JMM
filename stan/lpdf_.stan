real partial_sgt(array[] real y, int start, int end, vector mu, real sigma, real lambdap1half, real p, real q){
  real tot = 0;
  int k = 1;
  for (n in start:end){
    tot += sgt_lpdf(y[k] | mu[n], sigma, lambdap1half, p, q);
    k +=1;
  }
  return tot;
}

real partial_sym_gt(array[] real y, int start, int end, vector mu, real sigma, real p, real q){
  real tot = 0;
  int k = 1;
  for (n in start:end){
    tot += sym_gt_lpdf(y[k] | mu[n], sigma, p, q);
    k +=1;
  }
  return tot;
}

real partial_skewt(array[] real y, int start, int end, vector mu, real sigma, real lambdap1half, real qmhalf){
  real tot = 0;
  int k = 1;
  for (n in start:end){
    tot += constrained_skew_t_lpdf(y[k] | mu[n], sigma, lambdap1half, qmhalf);
    k +=1;
  }
  return tot;
}

