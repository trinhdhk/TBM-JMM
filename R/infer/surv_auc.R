clin_short_vols <- readRDS('data/imported/clin_vols.RDS' ) |>
  filter(grepl('(ventricular)|(entricl)|(Vent)', region)) |>
  mutate(vol = sum(vol), .by=c(timeInt, subject)) |>
  select(-region, -number) |>
  unique() |>
  filter(!is.na(vol))

mu_cd4 <- mean(log2(clin_short_vols[hiv==TRUE, cd4]), na.rm=TRUE)  
sd_cd4 <- sd(log2(clin_short_vols[hiv==TRUE, cd4]), na.rm=TRUE)

clin_short_vols <-
  clin_short_vols |>
  mutate(
    obj = order(subject),
    mrc = factor(mrc, levels=1:3, ordered = TRUE),
    muvol = mean(log(vol), na.rm=T),
    sdvol = sd(log(vol), na.rm=T),
    vol = scale(log(vol)) |> as.vector(),
    weight = scale(weight) |> as.vector(),
    csflym = scale(log10(csflym+1)) |> as.vector(),
    cd4 = ifelse(hiv, (log2(cd4) - mu_cd4)/sd_cd4, 0)) |>
  mutate(cd4 = replace_na(cd4, mu_cd4)) |>
  mutate(strata = paste(arm, as.numeric(hiv), mrc, genotype, sep='-')) |>
  mutate(
    severe=as.numeric(mrc)>1,
    genohiv = paste(hiv, genotype, sep='-'),
    sd_wk = sd(mri_wk, na.rm=TRUE),
    mu_wk = mean(mri_wk, na.rm=TRUE),
    mri_wk = mri_wk/sd_wk,
    arm = as.numeric(arm)-1) |>
  mutate(across(c(tt_death, tstart, tstop), ~ (.x)/sd_wk))

pred <-
  predict(icv_model_n$model,
          newdata = clin_short_vols,
          re_formula = NA,
          resp = 'evdeath')
cidx = discSurv::cIndex(marker = pred[,'Estimate'], 
       testTime = clin_short_vols$timeDisc, 
       testEvent = clin_short_vols$ev_death, 
       trainTime = clin_short_vols$timeDisc, 
       trainEvent = clin_short_vols$ev_death
    )

predlong <-
  predict(icv_model_n$model,
          newdata = icv_model_n$data |> filter(!is.na(vol)),
          re_formula = NA,
          resp = 'evdeath')

discSurv::calibrationPlot(
  predlong[,'Estimate'],
  icv_model_n$data |> filter(!is.na(vol)) |> mutate(y=ev_death),
)
