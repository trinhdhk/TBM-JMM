library(tidytable)
library(brms)

seed = 1208

# Helper function ----

icv.stanvar = function(model_code, family){
  pattern <- "// likelihood including constants(.*?)// priors including constants"
  likelihood <- regmatches(model_code, regexec(pattern, model_code))[[1]][2]
  # likelihood <- gsub('if (!prior_only) {', '', likelihood, fixed=TRUE)
  likelihood <- substr(likelihood, 1, nchar(likelihood)-4)
  
  lpdf_ <- 
    brms::stanvar(
      name='lpdf_',
      scode=readLines('stan/lpdf_.stan'),
      block='functions'
    )
  
  old_lpdf <- 
    switch(
      family,
      gaussian = 'target += normal_id_glm_lpdf(Yl_vol | X_vol, mu_vol, b_vol, sigma_vol)',
      skew_normal = 'target += skew_normal_lpdf(Yl_vol | mu_vol, omega_vol, alpha_vol)',
      t = 'target += student_t_lpdf(Yl_vol | nu_vol, mu_vol, sigma_vol)',
      st = 'for (n in 1:N_vol) {\n      target += constrained_skew_t_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol, lambdap1half_vol, qmhalf_vol);\n    }',
      sgt = 'for (n in 1:N_vol) {\n      target += sgt_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol, lambdap1half_vol, p_vol, q_vol);\n    }',
      sym_gt = 'for (n in 1:N_vol) {\n      target += sym_gt_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol, p_vol, q_vol);\n    }'
    )
  new_lpdf <-
    switch(
      family,
      gaussian = '
      mu_vol += X_vol * b_vol;
      target += normal_lpdf(Yl_vol[Jobs_vol] | mu_vol[Jobs_vol], sigma_vol);',
      skew_normal = '
      target += skew_normal_lpdf(Yl_vol[Jobs_vol] | mu_vol[Jobs_vol], omega_vol, alpha_vol);
      mu_vol += omega_vol * delta_vol * sqrt(2/pi());',
      t = 'target += student_t_lpdf(Yl_vol[Jobs_vol] | nu_vol, mu_vol[Jobs_vol], sigma_vol);',
      st = 'for (n in Jobs_vol) target += constrained_skew_t_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol, lambdap1half_vol, qmhalf_vol);',
      sgt = 'for (n in Jobs_vol) target += sgt_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol, lambdap1half_vol, p_vol, q_vol);',
      sym_gt = 'target += reduce_sum(partial_sym_gt, to_array_1d(Yl_vol[Jobs_vol]), 1, mu_vol[Jobs_vol], sigma_vol, p_vol, q_vol);'
    )
  
  surv.code <- '
  // mu_evdeath[Jmi_vol] += (bsp_evdeath[1]) * (mu_vol[Jmi_vol]-Yl_vol[Jmi_vol]) ;
  mu_evdeath += bsp_evdeath[1] * (mu_vol-Yl_vol) ;
  target += bernoulli_logit_glm_lupmf(Y_evdeath | Xc_evdeath, mu_evdeath, b_evdeath);
  target += std_normal_lupdf(Ymi_vol);'
  likelihood <- gsub('target += bernoulli_logit_glm_lpmf(Y_evdeath | Xc_evdeath, mu_evdeath, b_evdeath);', '', likelihood, fixed=TRUE)
  likelihood <- gsub(old_lpdf, paste(new_lpdf, surv.code, sep='\n'), likelihood, fixed=TRUE)
  likelihood <- gsub('Yl_vol[Jmi_vol] = Ymi_vol', 'Yl_vol[Jmi_vol] = mu_vol[Jmi_vol]', likelihood, fixed=TRUE)
  
  stopifnot(length(new_lpdf) > 0)
  tdata_code_def <- 'int Jobs_vol[N_vol - Nmi_vol];
  matrix[N_vol - Nmi_vol,K_vol] Xobs_vol;' |>
    brms::stanvar(name='obs_tdata_def', scode=_, block='tdata', position='start')
  tdata_code_end <- '
  {
    int k = 1;
    for (n in 1:N_vol) {
      int mi = 0;
      for (m in 1:Nmi_vol) {
       
        if (Jmi_vol[m] == n) {
          mi = 1;
          break;
        }
      }
      if (mi == 0){
        Jobs_vol[k] = n;
        k += 1;
      }
    }
    Xobs_vol = X_vol[Jobs_vol,:];
  }' |>  brms::stanvar(name='obs_tdata', scode=_, block='tdata', position='end')
  ll_1 <- '/*' |>  brms::stanvar(name='cancel_mid', scode=_, block='model', position='start')
  ll_2 <- '*/' |>  brms::stanvar(name='cancel_out', scode=_, block='likelihood', position='end')
  ll <- likelihood |> brms::stanvar(name='new', scode=_, block='likelihood', position='end')
  
  c(tdata_code_def, tdata_code_end, ll_1, ll_2, ll, lpdf_)
}

gq <- stanvar(
  name='n_gq',
  block='genquant',
  position='end',
  scode=
    '// vector combining observed and missing responses
    vector[N_vol] Yl_vol = Y_vol;
    // initialize linear predictor term
    vector[N_vol] mu_vol = rep_vector(0.0, N_vol);
    // initialize linear predictor term
    vector[N_evdeath] mu_evdeath = rep_vector(0.0, N_evdeath);
    vector[N_evdeath] H = rep_vector(0.0, N_evdeath);
    Yl_vol[Jmi_vol] = mu_vol[Jmi_vol];
    mu_vol += X_vol * b_vol;
    mu_evdeath += Intercept_evdeath + Xs_evdeath * bs_evdeath + Zs_evdeath_1_1 * s_evdeath_1_1 + Zs_evdeath_2_1 * s_evdeath_2_1;
    for (n in 1:N_vol) {
      // add more terms to the linear predictor
       mu_vol[n] += (bsp_vol[1]) * mo(simo_vol_1, Xmo_vol_1[n]) + (bsp_vol[2]) * mo(simo_vol_2, Xmo_vol_2[n]) * Csp_vol_1[n] + (bsp_vol[3]) * mo(simo_vol_3, Xmo_vol_3[n]) * Csp_vol_2[n] + r_1_vol_1[J_1_vol[n]] * Z_1_vol_1[n];
    }
    for (n in 1:N_evdeath) {
      // add more terms to the linear predictor
      mu_evdeath[n] += (bsp_evdeath[1]) * Yl_vol[n] + (bsp_evdeath[2]) * mo(simo_evdeath_1, Xmo_evdeath_1[n]) + (bsp_evdeath[3]) * mo(simo_evdeath_2, Xmo_evdeath_2[n]) * Csp_evdeath_1[n];
    }
    mu_evdeath += bsp_evdeath[1] * (mu_vol-Yl_vol) ;
    H = inv_logit(Xc_evdeath * b_evdeath + mu_evdeath);'
)

# Data -----
# approx <- 15
approx <- 24
source('R/02_prepare_data.R')
icv <- clin_vols |>
  filter(grepl('(ventricular)|(entricl)|(Vent)', region)) |>
  mutate(vol = sum(vol), .by=c(timeInt, obj)) |>
  select(-region, -number) |>
  unique() 

mu_cd4 <- mean(log2(icv[hiv==TRUE, cd4]), na.rm=TRUE)  
sd_cd4 <- sd(log2(icv[hiv==TRUE, cd4]), na.rm=TRUE)
icv <- mutate(icv,
              mrc = factor(mrc, levels=1:3, ordered = TRUE),
              muvol = mean(log(vol), na.rm=T),
              sdvol = sd(log(vol), na.rm=T),
              vol = scale(log(vol)) |> as.vector(),
              weight = scale(weight) |> as.vector(),
              csflym = scale(log10(csflym+1)) |> as.vector(),
              cd4 = ifelse(hiv, (log2(cd4) - mu_cd4)/sd_cd4, 0)) |>
  mutate(cd4 = replace_na(cd4, mu_cd4)) |>
  select(-ev_death) |>
  rename(ev_death=y) |>
  mutate(strata = paste(arm, as.numeric(hiv), mrc, genotype, sep='-')) |>
  mutate(
    severe=as.numeric(mrc)>1,
    genohiv = paste(hiv, genotype, sep='-'),
    sd_wk = sd(mri_wk, na.rm=TRUE),
    mu_wk = mean(mri_wk, na.rm=TRUE),
    mu_cd4 = mu_cd4,
    sd_cd4 = sd_cd4,
    mri_wk = mri_wk/sd_wk,
    arm = as.numeric(arm)-1,
    all_strata = paste0(mrc, arm, hiv, genotype, sep='-'),
    armhiv = paste(hiv, arm, sep='-')) |>
  mutate(across(c(tt_death, tstart, tstop), ~ (.x)/sd_wk))


# Formulae ----
surv_fml <-
  bf(ev_death ~ mi(vol) +
       s(tstop, k=3, by=genohiv, fx=TRUE) + 
       csflym + cd4 + arm + 
       arm:hiv + csflym:hiv, #+ 
       # mo(genotype) + mo(genotype):hiv, 
     family=bernoulli(link='logit')) 

icv_fml <- 
  bf(
    vol | mi() ~ 0+Intercept +
      s(mri_wk, k=3, by=armhiv, fx=TRUE) +
      mo(genotype) + mo(mrc) + hiv + 
      mo(mrc):mri_wk + mo(genotype):mri_wk + 
      arm:mo(mrc):mri_wk + arm:mo(genotype):mri_wk + 
      arm:mo(mrc):hiv:mri_wk + arm:mo(genotype):hiv:mri_wk + 
      # mrc:hiv:mri_wk + arm:mrc:hiv:mri_wk +
      # I(mrc>1):hiv:mri_wk +  arm:I(mrc>1):hiv:mri_wk +
      (1|obj)
  )

icv_fml <- 
  bf(
    vol | mi() ~ 0+Intercept +
      s(mri_wk, k=3, by=genohiv) +
      mo(mrc) + mo(mrc):mri_wk +
      arm:mri_wk + 
      arm:mo(mrc):mri_wk + arm:mo(genotype):mri_wk + arm:hiv:mri_wk + 
      arm:mo(mrc):hiv:mri_wk + arm:mo(genotype):hiv:mri_wk + 
      # arm:mo(mrc):mo(genotype):hiv:mri_wk + 
      # mrc:hiv:mri_wk + arm:mrc:hiv:mri_wk +
      # I(mrc>1):hiv:mri_wk +  arm:I(mrc>1):hiv:mri_wk +
      (1|obj)
  )

icv_fml <- 
  bf(
    vol | mi() ~ 0+Intercept +
      s(mri_wk, k=3, by=genohiv) +
      mo(mrc) + mo(mrc):mri_wk + cd4 +
      arm:mri_wk +
      arm:mo(mrc):mri_wk + arm:mo(genotype):mri_wk + arm:hiv:mri_wk + 
      arm:mo(mrc):hiv:mri_wk + arm:mo(genotype):hiv:mri_wk + 
      arm:cd4:mri_wk + arm:cd4:mo(mrc):mri_wk + 
      (1|obj)
  )

# icv_fml <- 
#   bf(
#     vol | mi() ~ 0+Intercept +
#       s(mri_wk, k=4, by=all_strata) + 
#       # mo(mrc):mri_wk +
#       (1|obj)
#   )



