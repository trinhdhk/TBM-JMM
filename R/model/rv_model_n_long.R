# Formula creating ----
source("R/model/rv_formulae.R")
sigma_fml <- sigma ~ s(zmuvol, k=3)  # s(zsdvol, k=3)

# Modeling ----

## Model specs ----
rv_longfml <-
  bf(
    update(rv_fml, vol ~ .),
    sigma_fml,
    family = gaussian()
  )

brmobj <-
  brm(
    rv_longfml,
    prior = c(
      prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class = "b"),
      # prior(exponential(1), class='sds'),
      prior(exponential(1), class = "sd"),
      prior(std_normal(), class = "Intercept", dpar = "sigma"),
      prior(std_normal(), class = "b", dpar = "sigma"),
      prior(exponential(1), class='sds', dpar='sigma')
    ),
    data = rv[!is.na(vol)],
    stanvars = c(
      sgtbrms::expose_sgt_stanvar()
    ),
    save_pars = save_pars(all = TRUE),
    cores = 1, chains = 0,
    empty=TRUE,
    save_warmup = TRUE,
    seed = seed,
    threads = threading(60),
    iter = 3000, warmup = 1000, thin = 1,
    refresh = 200, init_r = 1,
    control = list(adapt_delta = .85, max_treedepth = 12)
  )

stancode <- make_stancode(brmobj) |> 
  strsplit('\n') |>
  getElement(1)
mu_n <- grep('mu[n] +=', stancode, fixed=T) 
mu_nstr <- stancode[[mu_n]] |>
  strsplit('+=', fixed=T) |>
  getElement(1)
mu_nstr[[2]] <-
  gsub('[n]', '[nn]', x=mu_nstr[[2]], fixed=T)

stancode[[mu_n]] <-
  paste(mu_nstr, collapse='+=')

stancode <- paste(stancode, collapse='\n')

standata <- make_standata(brmobj)

stanmodel <- rstan::stan_model(model_code = stancode, 
                               auto_write = rstan::rstan_options("auto_write"))

rstan::rstan_options(threads_per_chain=60)
stanfit <- rstan::sampling(
  stanmodel,  
  data = standata,
  seed = seed,
  sample_file = ".cache/rv_model_n_long.csv",
  cores = 1, chains = 1,
  iter = 4000, warmup = 1000, thin = 1,
  refresh = 200, init_r = 1,
  control = list(adapt_delta = .85, max_treedepth = 12)
)

saveRDS(stanfit, file='results/rv_model_n_long.RDS' )
