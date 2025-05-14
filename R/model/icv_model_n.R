# Formula creating ----
source('R/model/icv_formulae.R')

# Modeling ----

## Model specs ----
icv_fml_n <- 
  bf(icv_fml,
     family=gaussian(link='identity')) +
  surv_fml + 
  set_rescor(FALSE)


## Make dummy code ----
icv_stancode_n <- 
  make_stancode(
    icv_fml_n,
    prior = c(
      prior(logistic(0,1), class='Intercept', resp='evdeath'),
      # prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
      prior(logistic(0,1), class='b', resp='evdeath'),
      prior(exponential(1), class='sds', resp='evdeath'),
      prior(exponential(1), class='sds', resp='vol'),
      
      prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      prior(exponential(1), class='sigma', resp='vol')
    ),
    data=icv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar()))

### Run real model ----
icv_model_n <-
  brm(
    icv_fml_n,
    prior = c(
      prior(logistic(0,1), class='Intercept', resp='evdeath'),
      prior(logistic(0,1), class='b', resp='evdeath'),
      # prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
      prior(exponential(1), class='sds', resp='evdeath'),
      prior(exponential(1), class='sds', resp='vol'),
      
      prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      prior(exponential(1), class='sigma', resp='vol')
    ),
    data=icv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar(), 
                 icv.stanvar(icv_stancode_n, 'gaussian')),
    
    save_pars = save_pars(all=TRUE),
    cores=2, chains=2,
    iter=6000, warmup=3000, thin=1, seed=seed, refresh=200, init_r=0.75,
    # algorithm='meanfield',
    # cores=1, chains=1,
    # tol_rel_obj = .0001, seed=seed, iter=40000, 
    control = list(adapt_delta=.85, max_treedepth=12)
  )

### Save model ----
saveRDS(list(data=icv, model=icv_model_n),
        file='results/icv_model_n_smooth2.RDS')
