# Formula creating ----
source('R/model/rv_formulae.R')
sigma_fml <- sigma ~ s(zmuvol, k=3)
q_fml <- q ~ s(zmuvol, k=3)
# Sys.setenv(STAN_NUM_THREADS=15)
# rstan::rstan_options(threads_per_chain=20)

source('R/helper/get_expect_value.R')
# Get init values ----
long.prefit <- readRDS('results/rv_model_symgt_long.RDS')
surv.prefit <- readRDS('results/rv_model_surv.RDS')

long.init <- get_expect_value(long.prefit)
surv.init <- get_expect_value(surv.prefit)

combined.init <- add_resp_to_par(list(surv.init, long.init), 
                                 resp=c('evdeath','vol'),
                                 re_count_init=1) |>
  unlist(recursive=F) #|>

combined.init$bsp_evdeath <- c(rnorm(1), combined.init$bsp_evdeath)
  # modifyList(
  #   list(sd_1=rexp(1),  z_1=rnorm(0.5))
  # )


# Modeling ----

## Model specs ----
rv_fml <- 
  surv_fml + 
  bf(rv_fml,
     sigma_fml,
     q_fml,
     family=sgtbrms::sym_gt(link='identity')) +
  set_rescor(FALSE)


## Make dummy code ----
rv_stancode <- 
  sgtbrms::make_stancode_sgt(
    rv_fml,
    prior = c(
      prior(logistic(0,1), class='Intercept', resp='evdeath'),
      prior(logistic(0,1), class='b', resp='evdeath'),
      # prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
      # prior(student_t(6, 0, 1), class='b', resp='evdeath'),
      prior(exponential(1), class='sds', resp='evdeath'),
      prior(exponential(1), class='sd', resp='evdeath'),
      
      prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      prior(exponential(1), class='sds', resp='vol'),
      # prior(student_t(6, 0, 1), class='b', resp='vol'),
      prior(exponential(1), class='sd', resp='vol'),
      
      prior(std_normal(), class='Intercept', dpar='sigma', resp='vol'),
      prior(std_normal(), class='b', dpar='sigma', resp='vol'),
      prior(exponential(1), class='sds', dpar='sigma', resp='vol'),
      prior(gamma(3, 0.1), class='p', resp='vol'),
      prior(std_normal(), class='Intercept', dpar='q', resp='vol'),
      prior(std_normal(), class='b', dpar='q', resp='vol'),
      prior(exponential(1), class='sds', dpar='q', resp='vol')
      # prior(gamma(3, 0.1), class='q', resp='vol')
    ),
    threads=2,
    data=rv)
    # stanvars = c(sgtbrms::expose_sgt_stanvar()))

### Run real model ----
rv_model <-
  brm(
    rv_fml,
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
      prior(exponential(1), class='sds', dpar='sigma', resp='vol'),
      prior(gamma(3, 0.1), class='p', resp='vol'),
      prior(std_normal(), class='Intercept', dpar='q', resp='vol'),
      prior(std_normal(), class='b', dpar='q', resp='vol'),
      prior(exponential(1), class='sds', dpar='q', resp='vol')
      # prior(gamma(3, 0.1), class='q', resp='vol')
      # prior(exponential(1), class='sds', dpar='sigma', resp='vol')
    ),
    data=rv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar(), 
                 rv_stanvar_thread(rv_stancode, 'sym_gt', n_subject_time)),
    threads=threading(60, 20), # 60 17
    save_pars = save_pars(all=TRUE),
    cores=1, chains=1,
    output_dir = '.cache',
    output_basename = 'rv_model_symgt_thread4',
    # sample_file = '.cache/rv_model_symgt_thread.csv',
    save_warmup=TRUE,
    seed=seed,
    iter=5000, warmup=1500, thin=1,
    refresh=200, 
    # init_r=0.5,
    # init='random',
    # init=0.5,
    init = list(combined.init),
    backend = 'cmdstanr',
    # algorithm='meanfield',
    # tol_rel_obj = .00001,
    # iter=20000, 
    control = list(adapt_delta=.73, max_treedepth=12)
  )

### Save model ----
saveRDS(list(data=rv, model=rv_model), file='results/rv_model_symgt4.RDS')
