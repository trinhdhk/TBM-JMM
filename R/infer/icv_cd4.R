condition2 <- expand.grid(
  genotype=c('CC', 'CT', 'TT'),
  hiv = T,
  # cd4 = obj$data[hiv==T, t(ggdist::median_qi(cd4, .width = .5)[, 1:3])][,1],
  cd4 = obj$data[hiv==T, seq(quantile(cd4,.2), quantile(cd4, .8), length.out=30)],
  # arm = c(0,1),
  mrc=1:3) |>
  as.data.frame() |>
  mutate(genohiv = paste(hiv, genotype, sep='-'))
icv_cd4 <-
  conditional_effects(
    obj$model, 
    'arm',
    conditions = condition2)$vol.vol_arm %>%
  filter(arm %in% c(0,1))

icv_cd4 <- dplyr::select(icv_cd4,-dplyr::ends_with('__'))
icv_cd4$.id <- 1:nrow(icv_cd4)
icv_cd4 <- 
  rbind(
    icv_cd4 |> dplyr::mutate(mri_wk = 0, t = 0),
    icv_cd4 |> dplyr::mutate(mri_wk = (60/7)/sd_wk, t = 1)
  ) |> tidybayes::add_epred_draws(obj$model, ndraws=500, resp='vol') 

icv_cd4 <- 
  icv_cd4 |>
  #dplyr::filter(.category=='vol') |>
  dplyr::mutate(.epred = .epred * sd_vol + mu_vol) |>
  dplyr::group_by(t) |>
  tidyr::pivot_wider(id_cols=c(.id, .draw, mrc, genotype, arm, cd4),
                     names_from='t', values_from='.epred', names_prefix='t.') |>
  dplyr::mutate(volChange=exp(t.1 - t.0),
                arm = ifelse(arm==0, 'Placebo', 'Dexamethasone')) |>
  dplyr::select(-.id, -t.0, -t.1) |>
  # dplyr::group_by(.draw, .mrc, genotype, hiv) |>
  tidyr::pivot_wider(id_cols = c(.draw, mrc, genotype, cd4),
                     names_from = 'arm', 
                     values_from = 'volChange') |>
  dplyr::mutate(tef = Dexamethasone/Placebo) |>
  mutate(cd4 = round(2^(cd4 * sd_cd4 + mu_cd4),0))

ggplot(icv_cd4,
  aes(y=tef, x=cd4)) +  
  ggdist::stat_lineribbon(alpha=.2) + facet_grid(~mrc) + 
  # ggdist::stat_pointinterval(position=position_dodge(.2)) + facet_grid(~genotype) + 
  # scale_x_continuous(trans='log10') +  
  scale_color_brewer(type='qual', palette=2)
