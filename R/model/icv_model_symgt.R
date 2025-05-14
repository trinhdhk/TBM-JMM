# Formula creating ----
source('R/model/icv_formulae.R')
Sys.setenv(STAN_NUM_THREADS=5)
rstan::rstan_options(threads_per_chain=5)

# Modeling ----

## Model specs ----
icv_fml <- 
  bf(icv_fml,
     family=sgtbrms::sym_gt(link='identity')) +
  surv_fml + 
  set_rescor(FALSE)


## Make dummy code ----
icv_stancode <- 
  sgtbrms::make_stancode_sgt(
    icv_fml,
    prior = c(
      prior(logistic(0,1), class='Intercept', resp='evdeath'),
      prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
      prior(exponential(1), class='sds', resp='evdeath'),
      
      prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      prior(exponential(1), class='sigma', resp='vol'),
      prior(gamma(3, 0.1), class='p', resp='vol'),
      prior(gamma(3, 0.1), class='q', resp='vol')
    ),
    data=icv)

### Run real model ----
icv_model <-
  brm(
    icv_fml,
    prior = c(
      prior(logistic(0,1), class='Intercept', resp='evdeath'),
      prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
      prior(exponential(1), class='sds', resp='evdeath'),
      
      prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      prior(exponential(1), class='sigma', resp='vol'),
      prior(gamma(3, 0.1), class='p', resp='vol'),
      prior(gamma(3, 0.1), class='q', resp='vol')
    ),
    data=icv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar(), icv.stanvar(icv_stancode, 'sym_gt')),
    
    save_pars = save_pars(all=TRUE),
    cores=2, chains=2,
    sample_file = '.cache/icv_model_symgt.csv', save_warmup=TRUE,
    iter=6000, warmup=3000, thin=1, seed=seed, refresh=200, init_r=0.75,
    control = list(adapt_delta=.88, max_treedepth=14)
  )

### Save model ----
saveRDS(list(data=icv, model=icv_model), file='results/icv_model_symgt.RDS')
