library(brms)
library(ggplot2)
library(tidytable)

obj <- icv_model_t2

condition <- expand.grid(genotype=c('CC', 'CT', 'TT'),
                         # cd4 = obj$data[hiv==T, t(ggdist::median_qi(cd4)[, 1:3])][,1],
                         hiv = c(F,T), mrc=1:3) |>
  as.data.frame() |>
  mutate(genohiv = paste(hiv, genotype, sep='-'))
icv_pred = conditional_effects(
  obj$model, 
  'arm',
  conditions = condition)$vol.vol_arm %>%
  filter(arm %in% c(0,1))
icv_pred <- dplyr::select(icv_pred,-dplyr::ends_with('__'))
icv_pred$.id <- 1:nrow(icv_pred)

sd_wk <- obj$data$sd_wk[[1]]
sd_vol <-  obj$data$sdvol[[1]]
mu_vol <- obj$data$muvol[[1]]
icv_pred <- 
  rbind(
    icv_pred |> dplyr::mutate(mri_wk = 0, t = 0),
    icv_pred |> dplyr::mutate(mri_wk = (60/7)/sd_wk, t = 1)
  ) |> tidybayes::add_epred_draws(obj$model, ndraws=500, resp='vol') 

icv_pred <- 
  icv_pred |>
  #dplyr::filter(.category=='vol') |>
  dplyr::mutate(.epred = .epred * sd_vol + mu_vol) |>
  dplyr::group_by(t) |>
  tidyr::pivot_wider(id_cols=c(.id, .draw, mrc, genotype, genohiv, arm, hiv), names_from='t', values_from='.epred', names_prefix='t.') |>
  dplyr::mutate(volChange=exp(t.1 - t.0),
                arm = ifelse(arm==0, 'Placebo', 'Dexamethasone')) |>
  dplyr::select(-.id, -t.0, -t.1) |>
  # dplyr::group_by(.draw, .mrc, genotype, hiv) |>
  tidyr::pivot_wider(id_cols = c(.draw, mrc, genotype, hiv),
                     names_from = 'arm', 
                     values_from = 'volChange') |>
  dplyr::mutate(tef = Dexamethasone/Placebo)

# ggplot(icv_pred |> filter(arm%in%c(0,1)) |> mutate(arm=factor(arm), 
#                                                    hiv=factor(hiv, labels=c("HIV-", "HIV+"))),
#        aes(y=volChange, x=genotype, color=arm)) + 
#   ggdist::stat_pointinterval(position=position_dodge(.2)) + facet_grid(hiv~mrc)


ggplot(icv_pred |> mutate(
                     hiv=factor(hiv, labels=c("HIV-", "HIV+"))),
       aes(y=tef, x=genotype, color=mrc)) +  
  ggdist::stat_pointinterval(position=position_dodge(.2)) + facet_grid(~hiv) + 
  scale_color_brewer(type='qual', palette=2)
