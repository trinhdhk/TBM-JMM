# Formula creating ----
source('R/model/rv_formulae.R')
sigma_fml <- sigma ~ s(zmuvol, k=3) #s(zsdvol, k=3)
# Sys.setenv(STAN_NUM_THREADS=5)
# rstan::rstan_options(threads_per_chain=20)
print(Sys.time())
source('R/helper/get_expect_value.R')
# Get init value
# long.prefit <- readRDS('results/rv_model_n_long.RDS')
# surv.prefit <- readRDS('results/rv_model_surv.RDS')

# long.init <- get_expect_value(long.prefit)
# surv.init <- get_expect_value(surv.prefit)

# combined.init <- add_resp_to_par(list(surv.init, long.init), 
                                #  resp=c('evdeath','vol'),
                                #  re_count_init=1) |>
  # unlist(recursive=F) #|>
  # modifyList(
  #   list(sd_1=rexp(1),  z_1=rnorm(0.5))
  # )

# Modeling ----

## Model specs ----
rv_fml_n <- 
  surv_fml + 
  bf(rv_fml,
     sigma_fml,
     family=gaussian(link='identity')) +
  set_rescor(FALSE)


## Make dummy code ----
rv_stancode_n_thread <- 
  make_stancode(
    rv_fml_n,
    prior = c(
      prior(logistic(0,1), class='Intercept', resp='evdeath'),
      prior(logistic(0,1), class='b', resp='evdeath'),
      # prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
      # prior(student_t(6, 0, 1), class='b', resp='evdeath'),
      prior(exponential(1), class='sds', resp='evdeath'),
      prior(exponential(1), class='sd', resp='evdeath'),
      
      prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      # prior(student_t(6, 0, 1), class='b', resp='vol'),
      prior(exponential(1), class='sds', resp='vol'),
      prior(exponential(1), class='sd', resp='vol'),
      
      prior(std_normal(), class='Intercept', dpar='sigma', resp='vol'),
      prior(std_normal(), class='b', dpar='sigma', resp='vol'),
      prior(exponential(1), class='sds', dpar='sigma', resp='vol')
    ),
    data=rv, 
    threads=2,
    stanvars = c(sgtbrms::expose_sgt_stanvar()))

### Run real model ----
rv_model_n <-
  brm(
    rv_fml_n,
    prior = c(
      prior(logistic(0,1), class='Intercept', resp='evdeath'),
      prior(logistic(0,1), class='b', resp='evdeath'),
      # prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
      # prior(student_t(6, 0, 1), class='b', resp='evdeath'),
      prior(exponential(1), class='sds', resp='evdeath'),
      prior(exponential(1), class='sd', resp='evdeath'),
      
      prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      # prior(student_t(6, 0, 1), class='b', resp='vol'),
      prior(exponential(1), class='sds', resp='vol'),
      prior(exponential(1), class='sd', resp='vol'),
      
      prior(std_normal(), class='Intercept', dpar='sigma', resp='vol'),
      prior(std_normal(), class='b', dpar='sigma', resp='vol'),
      prior(exponential(1), class='sds', dpar='sigma', resp='vol')
    ),
    data=rv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar(), 
                 rv_stanvar_thread(rv_stancode_n_thread, 'gaussian', n_subject_time)),
    
    threads=threading(60, 25),
    save_pars = save_pars(all=TRUE),
    cores=1, chains=1,
    output_dir = '.cache',
    output_basename = 'rv_model_n_thread4',
    save_warmup=TRUE,
    
    seed=seed,
    iter=6000, warmup=1000, thin=1,
    refresh=200, 
    # init_r=0.5,
    backend = 'cmdstanr',
    init=0.5,
    # init = list(combined.init),
    control = list(adapt_delta=.78, max_treedepth=10)
  )

### Save model ----
saveRDS(list(data=rv, model=rv_model_n), file='results/rv_model_n3.RDS')
