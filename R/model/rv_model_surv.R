# Formula creating ----
source("R/model/rv_formulae.R")

# make surv dataset
survdata <- rv[, head(.SD, 1), by=.(obj, tstop)]


brmobj <-
  brm(
    surv_fml0,
    prior = 
      c(
        prior(logistic(0,1), class='Intercept'),
        prior(logistic(0,1), class='b'),
        prior(exponential(1), class='sds')
      ),
    data=survdata,
    save_pars = save_pars(all = TRUE),
    cores = 1, chains = 0,
    empty=TRUE,
    seed = seed,
    thread=60,
    iter = 3000, warmup = 1000, thin = 1,
    refresh = 200, init_r = 0.75,
    control = list(adapt_delta = .85, max_treedepth = 12)
  )

stancode <- make_stancode(brmobj)
standata <- make_standata(brmobj)

stanmodel <- rstan::stan_model(model_code = stancode, 
                               auto_write = rstan::rstan_options("auto_write"))

rstan::rstan_options(threads_per_chain=60)

stanfit <- rstan::sampling(
  stanmodel,  
  data = standata,
  seed = seed,
  cores = 1, chains = 1,
  iter = 4000, warmup = 1000, thin = 1,
  refresh = 200, init_r = 0.75,
  control = list(adapt_delta = .85, max_treedepth = 12)
)

saveRDS(stanfit, file='results/rv_model_surv.RDS' )
