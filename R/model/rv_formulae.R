library(tidytable)
library(brms)
source('R/helper/stancode.R')
seed = 1208

# Helper function ----

rv_stanvar = function(model_code, family, n_subject_time){
  # Remove Ymi_vol from parameter
  model_code_split <-
    strsplit(model_code, '\n')[[1]]
  
  param_start <- grep('parameters \\{', model_code_split)[[1]] + 1
  param_end <- grep('transformed parameters \\{', model_code_split)[[1]] - 2
  param_code <- model_code_split[param_start:param_end]
  param_code <- param_code[!grepl("Ymi_", param_code, fixed=TRUE)]
  param_code <- paste(param_code, collapse='\n')
  
  # ---
  # Remove Ymi_vol from likelihood
  ll_start <- grep('// likelihood including constants', model_code_split)[[1]] + 2
  ll_end <- grep('// priors including constants', model_code_split)[[1]] - 2
  ll_code <- model_code_split[ll_start:ll_end]
  Yl_def <- grep('Ymi_', ll_code, fixed=TRUE, value=TRUE)
  newYl_def <- gsub('Ymi_vol', 'rep_vector(0, Nmi_vol)', Yl_def)
  
  pattern <- "// likelihood including constants(.*?)// priors including constants"
  likelihood <- regmatches(model_code, regexec(pattern, model_code))[[1]][2]
  # likelihood <- gsub('if (!prior_only) {', '', likelihood, fixed=TRUE)
  likelihood <- substr(likelihood, 1, nchar(likelihood)-4)
  likelihood <- gsub(Yl_def, newYl_def, likelihood, fixed=TRUE)
  
  # adding subject-time list to stan
  subject_time <-
    brms::stanvar(
      x=n_subject_time,
      name='n_subject_time',
      block='data'
    )
  
  surv_y_def <-
    brms::stanvar(
      name = 'surv_y_def',
      block = 'tdata',
      scode =' array[n_subject_time] int Y_evdeath_subject;',
      position = 'start'
    )
  
  surv_y_assign <- 
    brms::stanvar(
      name='surv_y_assign',
      block='tdata',
      position = 'end',
      scode='
      for (subject_time in 1:n_subject_time){
        int i_start = (subject_time - 1)*72 + 1;
        int i_stop = subject_time*72;
        Y_evdeath_subject[subject_time] = Y_evdeath[i_start];
        if (Y_evdeath[i_start] != mean(Y_evdeath[i_start:i_stop])) reject("BUG found!");
      }
      '
    )
  
  surv_lpmf <- 
    brms::stanvar(
      name='surv_lpmf',
      scode=readLines('stan/rv_functions.stan'),
      block='functions'
    )
  
  # Add cholesky for correlation
  # chol_param <- 
  #   brms::stanvar(
  #     name='chol_param',
  #     scode='
  #     //corr_matrix[N_1] Sigma;
  #     cholesky_factor_corr[N_1] L_Omega; //correlation matrix of region effect on death',
  #     block='parameters'
  #   )
  
  # chol_prior <-
  #   brms::stanvar(
  #     name='chol_prior',
  #     scode='lprior += lkj_corr_cholesky_lpdf(L_Omega | 4);',
  #     block='tparameters',
  #     position='end'
  #   )
  
  # chol_z1_lpdf <- # Assuming z1 is the random effect of region on survival. 
  #   brms::stanvar(
  #     name='chol_z1_lpdf',
  #     block='model',
  #     position='end',
  #     scode='
  #     //target += lkj_corr_lupdf(Sigma | 4);
  #     //target += multi_normal_lupdf(z_1[1] | rep_vector(0, N_1), Sigma) - std_normal_lpdf(z_1[1]);
  #     target += lkj_corr_cholesky_lupdf(L_Omega | 4);
  #     target += multi_normal_cholesky_lupdf(z_1[1]| rep_vector(0, N_1), L_Omega) - std_normal_lpdf(z_1[1]);'
  #   )
  
  old_lpdf <- 
    switch(
      family,
      gaussian = 'target += normal_lpdf(Yl_vol | mu_vol, sigma_vol)',
      skew_normal = 'target += skew_normal_lpdf(Yl_vol | mu_vol, omega_vol, alpha_vol)',
      t = 'target += student_t_lpdf(Yl_vol | nu_vol, mu_vol, sigma_vol)',
      st = 'for (n in 1:N_vol) {\n      target += constrained_skew_t_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], lambdap1half_vol, qmhalf_vol);\n    }',
      sgt = 'for (n in 1:N_vol) {\n      target += sgt_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], lambdap1half_vol, p_vol, q_vol);\n    }',
      sym_gt = 'for (n in 1:N_vol) {\n      target += sym_gt_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], p_vol, q_vol);\n    }'
    )
  new_lpdf <-
    switch(
      family,
      # gaussian = 'target += normal_lupdf(Yl_vol[Jobs_vol] | mu_vol[Jobs_vol], sigma_vol[Jobs_vol]);',
      gaussian = 'target += reduce_sum(partial_normal_lpdf, to_array_1d(Yl_vol[Jobs_vol]), 1, mu_vol[Jobs_vol], sigma_vol[Jobs_vol]);',
      skew_normal = '
      target += skew_normal_lupdf(Yl_vol[Jobs_vol] | mu_vol[Jobs_vol], omega_vol[Jobs_vol], alpha_vol);
      mu_vol += omega_vol * delta_vol * sqrt(2/pi());',
      # t = 'target += student_t_lupdf(Yl_vol[Jobs_vol] | nu_vol, mu_vol[Jobs_vol], sigma_vol[Jobs_vol]);',
      t = 'target += reduce_sum(partial_t_lpdf, to_array_1d(Yl_vol[Jobs_vol]), 1, nu_vol[Jobs_vol], mu_vol[Jobs_vol], sigma_vol[Jobs_vol]);',
      st = 'target += reduce_sum(partial_skewt, to_array_1d(Yl_vol[Jobs_vol]), 1, mu_vol[Jobs_vol], sigma_vol[Jobs_vol], lambdap1half_vol, qmhalf_vol);',
      sgt = 'target += reduce_sum(partial_sgt, to_array_1d(Yl_vol[Jobs_vol]), 1, mu_vol[Jobs_vol], sigma_vol[Jobs_vol], lambdap1half_vol, p_vol, q_vol);',
      sym_gt = 'target += reduce_sum(partial_sym_gt, to_array_1d(Yl_vol[Jobs_vol]), 1,  mu_vol[Jobs_vol], sigma_vol[Jobs_vol], p_vol, q_vol);'
      # st = 'for (n in Jobs_vol) target += constrained_skew_t_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], lambdap1half_vol, qmhalf_vol);',
      # sgt = 'for (n in Jobs_vol) target += sgt_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], lambdap1half_vol, p_vol, q_vol);',
      # sym_gt = 'for (n in Jobs_vol) target += sym_gt_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], p_vol, q_vol);'
    )
  
  # Separate fixed mu_vol and random region effect to mu_vol
  old_vol_region_re <- 
    likelihood |>
    strsplit('\n') |> 
    getElement(1) |> 
    grep('mu_vol\\[n\\] \\+\\=', x=_, value=T) 
  
  new_vol_region_re <- 
    old_vol_region_re |> 
    gsub('mu_vol\\[n\\] \\+= ', 'mu_vol_fixed[n] += mu_vol[n] +' ,x=_) |> 
    # gsub('r_1_vol_1[J_1_vol[n]] * Z_1_vol_1[n]', '',x=_, fixed=T)
    gsub('r_2_vol_1[J_2_vol[n]] * Z_2_vol_1[n] +', 
         'r_2_vol_1[J_2_vol[n]] * Z_2_vol_1[n];\n mu_vol[n] = mu_vol_fixed[n] +',
         x=_, fixed=T) 
  
  likelihood <- gsub(old_vol_region_re, new_vol_region_re, likelihood, fixed=TRUE)
  
  surv.code <- glue::glue('
  //mu_evdeath += (bsp_evdeath[1] + r_1_evdeath_1[J_1_evdeath]) .* (-Yl_vol) +  (bsp_evdeath[1] * exp(r_1_evdeath_1[J_1_evdeath]-r_1_evdeath_1[J_1_evdeath[1]])) .* mu_vol;
  //mu_evdeath += bsp_evdeath[1] .* (mu_vol[Jmi_vol]-Yl_vol[Jmi_vol]);
  
  //mu_evdeath +=  bsp_evdeath[1] .* (-Yl_vol);
  //mu_evdeath += r_1_evdeath_1[J_1_evdeath] .* (mu_vol-Yl_vol);
  
  mu_evdeath += bsp_evdeath[1] .* (mu_vol_fixed - Yl_vol);
  mu_evdeath += 72 * r_1_evdeath_1[J_1_evdeath] .* (mu_vol - mu_vol_fixed - Yl_vol);
  //0302 mu_evdeath += bsp_evdeath[1] .* (mu_vol - Yl_vol);
  
  //target += bernoulli_logit_glm_lpmf(Y_evdeath | Xc_evdeath, mu_evdeath, b_evdeath);
  //target += surv_lpmf(Y_evdeath_subject | Xc_evdeath, mu_evdeath, b_evdeath, mu_vol, bsp_evdeath[1], n_subject_time);
  
  //target += surv_glm2_lpmf(Y_evdeath_subject | Xc_evdeath, mu_evdeath, b_evdeath, mu_vol, bsp_evdeath[1], n_subject_time);
  //target += surv_glm_lpmf(Y_evdeath_subject | Xc_evdeath, mu_evdeath, b_evdeath, n_subject_time);
  target += reduce_sum(partial_surv_glm_lpmf, Y_evdeath_subject, 1, Xc_evdeath, mu_evdeath, b_evdeath);
  
  // 0302 target += bernoulli_logit_glm_lpmf(Y_evdeath | Xc_evdeath, mu_evdeath, b_evdeath)/72;
  //target += std_normal_lupdf(Ymi_vol);')
  likelihood <- gsub('target += bernoulli_logit_glm_lpmf(Y_evdeath | Xc_evdeath, mu_evdeath, b_evdeath);', '', likelihood, fixed=TRUE)
  likelihood <- gsub(old_lpdf, paste(new_lpdf, surv.code, sep='\n'), likelihood, fixed=TRUE)
  likelihood <- gsub('Yl_vol[Jmi_vol] = Ymi_vol', 'Yl_vol[Jmi_vol] = mu_vol[Jmi_vol]', likelihood, fixed=TRUE)
  
  stopifnot(length(new_lpdf) > 0)
  tdata_code_def <- 'array[N_vol - Nmi_vol] int Jobs_vol;
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
  ll_2 <- '*/\n vector[N_vol] mu_vol_fixed = rep_vector(0.0, N_vol);\n' |>  brms::stanvar(name='cancel_out', scode=_, block='likelihood', position='end')
  ll <- likelihood |> brms::stanvar(name='new', scode=_, block='likelihood', position='end')

  # very annoying that cannot inject to start of parameters so had to change internal function of brms
  param_cancel_in <- brms::stanvar(name='param_cancel_in', scode='/*', block='parameters', position='start')
  param_cancel_out <- brms::stanvar(name='param_cancel_out', scode='*/\n', block='parameters', position='end')
  new_param  <- brms::stanvar(name='param_code', scode=param_code, block='parameters', position='end')
  c(tdata_code_def, tdata_code_end, surv_y_assign,
    param_cancel_in, param_cancel_out, new_param,
    ll_1, ll_2, ll, subject_time, surv_lpmf,
    surv_y_def)
    # chol_param, 
    # chol_prior,
    # chol_z1_lpdf)
}

get_param_of_fn <- function(fn_def){
  regmatches(fn_def, regexec('\\(.*?\\)', fn_def))[[1]] |> 
    gsub('\\(|\\)','',x=_) |>
    strsplit(',\\s?') |> getElement(1)
}

get_random_term <-
  function(re_string){
    stringr::str_extract_all(re_string, 'r_\\d+.*?\\d+')[[1]]
  }

get_random_var <- 
  function(re_string, id){
    stringr::str_extract_all(re_string, glue::glue('r_{id}+.*?\\d+\\[J_{id}.*?\\[.*?\\]]'))[[1]]
  }

remove_random_term <-
  function(re_string, id){
    re_split <- strsplit(re_string, '\\s*?\\+\\s*?')[[1]]
    new_str <- sapply(re_split, 
           \(rr){
             gsub(glue::glue('r_{id}+.*?\\d+\\[J_{id}.*?\\[.*?\\]]'), '0', rr) |>
               gsub('0\\s*?\\*.*$', '', x=_)
           }) 
    new_str[!grepl('^\\s*?$', new_str)] |>
      paste(collapse=' + ')
  }

get_max_id <-
  function(rterm){
    strsplit(rterm, '_') |> 
      sapply(getElement,2) |> as.numeric() |> max()
  }

rv_stanvar_thread = function(model_code, family, n_subject_time){
  
  model_code_split <-
    strsplit(model_code, '\n')[[1]]
  
  # Remove Ymi_vol from parameter
  param_start <- grep('parameters \\{', model_code_split)[[1]] + 1
  param_end <- grep('transformed parameters \\{', model_code_split)[[1]] - 2
  param_code <- model_code_split[param_start:param_end]
  param_code <- param_code[!grepl("Ymi_", param_code, fixed=TRUE)]
  param_code <- paste(param_code, collapse='\n')
  
  # very annoying that cannot inject to start of parameters so had to change internal function of brms
  param_cancel_in <- brms::stanvar(name='param_cancel_in', scode='/*', block='parameters', position='start')
  param_cancel_out <- brms::stanvar(name='param_cancel_out', scode='*/\n', block='parameters', position='end')
  new_param  <- brms::stanvar(name='param_code', scode=param_code, block='parameters', position='end')
  
  # ---
  # Remove Ymi_vol from likelihood
  ll_start <- grep('// likelihood including constants', model_code_split)[[1]] + 1
  ll_end <- grep('generated quantities \\{', model_code_split)[[1]] - 2
  ll_code <- model_code_split[ll_start:ll_end]
  Yl_def <- grep('Ymi_', ll_code, fixed=TRUE)
  ll_code[Yl_def] <- gsub('Ymi_vol', 'rep_vector(0, Nmi_vol)', ll_code[[Yl_def]])
  
  reduce_sum <- grep('reduce_sum', ll_code, fixed=TRUE, value=TRUE)
  reduce_par <-
    sapply(
      rev(reduce_sum),
      \(x) get_param_of_fn(x)[-1:-2]
    ) |> unlist() |> unique() 
  reduce_par <- c('joint_model_log_lik_lpmf, seq_subject_time', reduce_par[1], 'I_obs, Y_evdeath_subject', reduce_par[-1])
  reduce_call <-
    paste(
      'target += reduce_sum(',
      paste(reduce_par, collapse=', '),
      ');'
    )
  
  reduce_sum_pos <- grep('reduce_sum', ll_code, fixed=TRUE)[[1]]-1
  ll_code <- ll_code[-grep('reduce_sum', ll_code, fixed=TRUE)]
  ll_code <- c(ll_code[1:reduce_sum_pos], reduce_call, ll_code[(1+reduce_sum_pos):length(ll_code)])
  ll_1 <- '/*' |>  brms::stanvar(name='cancel_mid', scode=_, block='model', position='start')
  ll_2 <- '*/\n ' |>  brms::stanvar(name='cancel_out', scode=_, block='model', position='end')
  ll <- ll_code |>
    brms::stanvar(name='new', scode=_, block='model', position='end')
  
  # ---
  # manipulate partial functions
  fn_start <- grep('real partial_log_lik_', model_code_split)
  fn_end <- grep('return ptarget', model_code_split)
  fn_param <- 
    sapply(
      rev(seq_along(fn_start)),
      function(i){
        fn_start <- fn_start[[i]]
        get_param_of_fn(model_code_split[[fn_start]])
      }
    ) |> unlist() |> unique()
  fn_param <- fn_param[!grepl('seq\\_', fn_param)]
  fn_param <- c('array[] int seq_subject_time', fn_param[1:2], 'array[] int I_obs', 'array[] int Y_evdeath_subject',fn_param[-1:-2])
  fn_new <- paste(
    'real joint_model_log_lik_lpmf(',
    paste(fn_param, collapse=', '),
    ')'
  )
  
  fn_content <-
    sapply(
      rev(seq_along(fn_start)),
      function(i){
        fn_start <- fn_start[[i]]
        fn_end <- fn_end[[i]]
        code <- model_code_split[(fn_start+1):fn_end]
        # if (i == 2){
          code <- lapply(code, gsub, pattern='start', replacement='Start', fixed=T) |>
            lapply(gsub, pattern='end', replacement='End', fixed=T) |>
            unlist()
        # }
        code
      }
    ) |> unlist()
  fn_def_i <- (grepl('^\\s*(real|int|vector|array|matrix)', fn_content) & !grepl('int nn', fn_content))
  fn_return_i <- grepl('return ptarget', fn_content) 
  fn_def <- fn_content[fn_def_i] |> unique() 
  fn_def <- fn_def[!grepl('int N = ',fn_def, fixed=T)]
  fn_return <- fn_content[fn_return_i] |> unique()
  fn_cal <- fn_content[!fn_def_i & !fn_return_i]
  
  # We will process fn_cal, where all the work stay
  # change the lpdf
  
  # add the def for Jobs_vol
  fn_def <- c(
    fn_def,
    'array[sum(I_obs[Start:End])] int Jobs_vol = which(I_obs[Start:End]);'
  )
  
  old_lpdf <- 
    switch(
      family,
      gaussian = 'ptarget += normal_lpdf(Yl_vol[Start:End] | mu_vol, sigma_vol);',
      skew_normal = 'ptarget += skew_normal_lpdf(Yl_vol[Start:End] | mu_vol, omega_vol, alpha_vol);',
      t = 'ptarget += student_t_lpdf(Yl_vol[Start:End] | nu_vol, mu_vol, sigma_vol);',
      st = 'for (n in 1:N_vol) {\n      ptarget += constrained_skew_t_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], lambdap1half_vol, qmhalf_vol);\n    }',
      sgt = 'for (n in 1:N_vol) {\n      ptarget += sgt_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], lambdap1half_vol, p_vol, q_vol);\n    }',
      sym_gt = 'ptarget += sym_gt_lpdf(Yl_vol[nn] | mu_vol[n], sigma_vol[n], p_vol, q_vol[n]);'
    )
  new_lpdf <-
    switch(
      family,
      # gaussian = 'target += normal_lupdf(Yl_vol[Jobs_vol] | mu_vol[Jobs_vol], sigma_vol[Jobs_vol]);',
      gaussian = 'ptarget += normal_lupdf(Yl_vol[Start:End][Jobs_vol] | mu_vol[Jobs_vol], sigma_vol[Jobs_vol]);',
      skew_normal = '
      ptarget += skew_normal_lupdf(Yl_vol[Start:End][Jobs_vol] | mu_vol[Jobs_vol], omega_vol[Jobs_vol], alpha_vol);
      mu_vol += omega_vol * delta_vol * sqrt(2/pi());',
      t = 'ptarget += student_t_lupdf(Yl_vol[Start:End][Jobs_vol] | nu_vol[Jobs_vol], mu_vol[Jobs_vol], sigma_vol[Jobs_vol]);',
      st = 'ptarget += constrained_skew_t_lpdf(Yl_vol[Jobs_vol] | mu_vol[Jobs_vol], sigma_vol[Jobs_vol], lambdap1half_vol, qmhalf_vol);',
      sgt = 'ptarget += sgt_lpdf(Yl_vol[Jobs_vol] | mu_vol[Jobs_vol], sigma_vol[Jobs_vol], lambdap1half_vol, p_vol, q_vol);',
      sym_gt = 'ptarget += sym_gt_lpdf(Yl_vol[Start:End][j] | mu_vol[j], sigma_vol[j], p_vol, q_vol[j]);'
    )
  
  old_ll <- grep(old_lpdf, fn_cal, fixed=TRUE)
  fn_cal[old_ll] <- gsub(old_lpdf, new_lpdf, fn_cal[old_ll], fixed=TRUE)
  if (family %in% c('st', 'sgt', 'sym_gt')){
    fn_cal[old_ll-2] <-  '    for (j in Jobs_vol) {'
    fn_cal[old_ll-1] <- ''
  }
  
  # Separate fixed mu_vol and random region effect to mu_vol
  old_vol_region_re <- 
    fn_cal |> 
    grep('mu_vol\\[n\\] \\+\\=', x=_) 
  ## Get all random terms and determine the max id == id of region
  rhs <- strsplit(fn_cal[old_vol_region_re],'\\s\\+\\=\\s')[[1]][[2]]
  rterms <- get_random_term(rhs)
  max_id <- get_max_id(rterms)
  # replace all random term with 0
  fixed_part <- paste0('mu_vol_fixed[n] = mu_vol[n] + ', remove_random_term(rhs, max_id))
  if (substr(fixed_part, nchar(fixed_part), nchar(fixed_part)) != ';'){
    fixed_part <- paste0(fixed_part, ';')
  } 
  
  fn_cal[old_vol_region_re] <- 
    fn_cal[old_vol_region_re] |> 
    # Fix a bug in random effect of monotonic effects
    gsub('\\[n\\]', '[nn]', x=_) |>
    gsub('mu_vol\\[nn\\] \\+= ', 'mu_vol[n] += ' ,x=_) %>%
    paste0(fixed_part, '\n',.)

  
  # fn_cal[old_vol_region_re] <- 
  #   fn_cal[old_vol_region_re] |> 
  #   # Fix a bug in random effect of monotonic effects
  #   gsub('\\[n\\]', '[nn]', x=_) |>
  #   gsub('mu_vol\\[nn\\] \\+= ', 'mu_vol_fixed[n] += mu_vol[n] +' ,x=_) |> 
  #   # gsub('r_1_vol_1[J_1_vol[n]] * Z_1_vol_1[n]', '',x=_, fixed=T)
  #   gsub('r_2_vol_1[J_2_vol[nn]] * Z_2_vol_1[nn] +', 
  #        'r_2_vol_1[J_2_vol[nn]] * Z_2_vol_1[nn];\n mu_vol[n] = mu_vol_fixed[n] +',
  #        x=_, fixed=T) 
  
  
  # add a definition for mu_vol_fixed
  fn_def <- c(fn_def, 
              'vector[N_vol] mu_vol_fixed = rep_vector(0.0, N_vol);')
  
  # Change the surv code
  # add def for y_evdeath
  fn_def <- c(
    'int n_region = 72;',
    'int Start = (start - 1)*n_region + 1;',
    'int End = end*n_region;',
    'array[end-start+1] int y_evdeath_subject = Y_evdeath_subject[start:end];',
    fn_def
  )
  
  surv_new <-'mu_evdeath += bsp_evdeath[1] .* (mu_vol_fixed); \n
  mu_evdeath += 72 * r_1_evdeath_1[J_1_evdeath[Start:End]] .* (mu_vol - mu_vol_fixed); \n
  ptarget += surv_glm_lpmf(y_evdeath_subject'

  surv_mu <- grep('mu_evdeath[n] +=', fn_cal, fixed=TRUE)
  fn_cal[surv_mu] <- ''
  surv_old <- grep('ptarget += bernoulli_logit_glm_lpmf(Y_evdeath[Start:End]', fn_cal, fixed=T)
  fn_cal[surv_old] <- 
    gsub('ptarget += bernoulli_logit_glm_lpmf(Y_evdeath[Start:End]', surv_new, fn_cal[surv_old], fixed=TRUE)
  
  # return the content
  fn_content <- c(fn_def, 
                  fn_cal,
                  fn_return)
  
  fn_string <- paste(
    c(fn_new,
      '{',
      fn_content,
      '}'),
    collapse='\n'
  )
  
  partial_lpmf <- 
    brms::stanvar(
      name='partial_lpmf',
      scode=fn_string,
      block='functions'
    )
  
  
  # adding subject-time list to stan
  subject_time <-
    brms::stanvar(
      x=n_subject_time,
      name='n_subject_time',
      block='data'
    )
  
  surv_y_def <-
    brms::stanvar(
      name = 'surv_y_def',
      block = 'tdata',
      scode =' array[n_subject_time] int Y_evdeath_subject;',
      position = 'start'
    )
  
  surv_y_assign <- 
    brms::stanvar(
      name='surv_y_assign',
      block='tdata',
      position = 'end',
      scode='
      for (subject_time in 1:n_subject_time){
        int i_start = (subject_time - 1)*72 + 1;
        int i_stop = subject_time*72;
        Y_evdeath_subject[subject_time] = Y_evdeath[i_start];
        if (Y_evdeath[i_start] != mean(Y_evdeath[i_start:i_stop])) reject("BUG found!");
      }
      '
    )
  
  tdata_code_def <- 'array[N_vol - Nmi_vol] int Jobs_vol= which_obs(Jmi_vol, N_vol);
  matrix[N_vol - Nmi_vol,K_vol] Xobs_vol;
  array[N] int I_obs = is_obs(Jmi_vol, N_vol);
  array[n_subject_time] int seq_subject_time = sequence(1, n_subject_time);' |>
    brms::stanvar(name='obs_tdata_def', scode=_, block='tdata', position='start')
  tdata_code_end <- '
    Xobs_vol = X_vol[Jobs_vol,:];' |>  
    brms::stanvar(name='obs_tdata', scode=_, block='tdata', position='end')
  
  
  # Add some functions
  surv_lpmf <- 
    brms::stanvar(
      name='surv_lpmf',
      scode=readLines('stan/rv_functions.stan'),
      block='functions'
    )
  
  obs_function <-
    brms::stanvar(
      name='obs_function',
      scode = '
        array[] int which(array[] int x){
          array[sum(x)] int out;
          int k = 1;
          for (i in 1:num_elements(x)){
            if (x[i] == 1){
              out[k] = i;
              k += 1;
            }
          }
          return out;
        }
      
        array[] int is_obs(array[] int Jmi, int N){
          array[N] int I_obs = rep_array(1, N);
          int Nmi = num_elements(Jmi);
          for (n in 1:N){
            for (m in 1:Nmi){
              if (Jmi[m] == n){
                I_obs[n] = 0;
                break;
              }
            }
          }
          return I_obs;
        }
        
        array[] int which_obs(array[] int Jmi, int N){
          int k = 1;
          int Nmi = num_elements(Jmi);
          array[N-Nmi] int Jobs;
          for (n in 1:N){
            int mi = 0;
            for (m in 1:Nmi){
              if (Jmi[m] == n){
                mi = 1;
                break;
              }
            }
              
            if (mi == 0){
              Jobs[k] = n;
              k += 1;
            }
          }
          return Jobs;
        }
      ',
      block='functions'
    )
 
  
  
  c(tdata_code_def, tdata_code_end, surv_y_assign,
    param_cancel_in, param_cancel_out, new_param,
    ll_1, ll_2, ll, subject_time, surv_lpmf,
    obs_function, partial_lpmf,
    surv_y_def)
  # chol_param, 
  # chol_prior,
  # chol_z1_lpdf)
}

# Data -----
approx <- 20
source('R/02_prepare_data.R')
rv <- 
  clin_vols |>
  # readRDS('data/imported/clin_vols.RDS') |>
  filter(!grepl('(ventricular)|(entricl)|(Vent)', region)) |>
  mutate(vol = sum(vol), .by=c(timeInt, obj, region)) |>
  select(-number) |>
  unique() 

mu_cd4 <- mean(log2(rv[hiv==TRUE, cd4]), na.rm=TRUE)  
sd_cd4 <- sd(log2(rv[hiv==TRUE, cd4]), na.rm=TRUE)
rv <- mutate(rv,
              weight = scale(weight) |> as.vector(),
              csflym = scale(log10(csflym)) |> as.vector(),
              cd4 = ifelse(hiv, (log2(cd4) - mu_cd4)/sd_cd4, 0)) |>
  mutate(cd4 = replace_na(cd4, mu_cd4)) |>
  select(-ev_death) |>
  rename(ev_death=y) |>
  mutate(strata = paste(arm, as.numeric(hiv), mrc, genotype, sep='-')) |>
  mutate(
    severe=as.numeric(mrc)>1,
    # mrc = as.numeric(mrc)>1, #factor(mrc, ordered=FALSE),
    sd_wk = sd(mri_wk, na.rm=TRUE),
    mu_wk = mean(mri_wk, na.rm=TRUE),
    mu_cd4 = mu_cd4,
    sd_cd4 = sd_cd4,
    genohiv = paste(hiv, genotype, sep='-'),
    mri_wk = mri_wk/sd_wk
    # arm = as.numeric(arm)-1
    ) |>
  mutate(across(c(tt_death, tstart, tstop), ~ (.x)/sd_wk)) |>
  mutate(
    muvol = mean(log(vol), na.rm=T),
    sdvol = sd(log(vol), na.rm=T),
    vol = scale(log(vol)) |> as.vector(),
    .by = region
  ) |>
  mutate(
    zmuvol = scale(muvol),
    zsdvol = scale(sdvol)
  ) |>
  arrange(obj, timeInt)

# get the subject time list
subject_time <- with(rv, paste(obj, timeInt)) 
n_subject_time <- length(unique(subject_time))

# Formulae ----
surv_fml <-
  bf(ev_death ~ mi(vol) +
       (0+mi(vol)|region) +
       s(tstop, k=5, fx=TRUE) + 
       hiv + 
       mo(genotype) + 
       csflym + cd4 + 
       arm + 
       arm:hiv + 
       arm:csflym +
       arm:mo(genotype) +
       csflym:hiv +
       mo(genotype):hiv,
      cmc=FALSE,
     family=bernoulli(link='logit')) 

# without vol
surv_fml0 <-
  bf(ev_death ~ 
       s(tstop, k=5, fx=TRUE) + 
       hiv + 
       mo(genotype) + 
       csflym + cd4 + 
       arm + 
       arm:hiv + 
       arm:csflym +
       arm:mo(genotype) +
       csflym:hiv +
       mo(genotype):hiv,
       cmc=FALSE,
     # mo(genotype) + mo(genotype):hiv, 
     family=bernoulli(link='logit')) 

rv_fml <- 
  bf(
    vol | mi() ~ 0+Intercept +
      # s(mri_wk, k=3)  +
      mri_wk +
      mo(genotype) + 
      mo(mrc) + 
      hiv +
      cd4 +
      
      mo(genotype):hiv + 
      mo(mrc):hiv + 
      mo(genotype):mo(mrc):hiv +
      
      hiv:mri_wk + 
      mo(mrc):mri_wk +
      mo(genotype):mri_wk +
      mo(mrc):mo(genotype):mri_wk + 
      arm:mri_wk +
      cd4:mri_wk +
      
      arm:mo(mrc):mri_wk + 
      arm:mo(genotype):mri_wk + 
      arm:hiv:mri_wk + 
      arm:cd4:mri_wk +
      mo(mrc):hiv:mri_wk +
      mo(genotype):hiv:mri_wk +
      
      arm:mo(mrc):hiv:mri_wk +
      arm:mo(genotype):hiv:mri_wk +
      arm:mo(genotype):mo(mrc):mri_wk +
      (1|obj) + 
      # s(mri_wk, k=3, by=region, bs='ts')+
      (1 + mri_wk + 
         mo(genotype) + 
         mo(mrc) + 
         hiv +
         cd4 +
         
         mo(genotype):hiv + 
         mo(mrc):hiv + 
         mo(genotype):mo(mrc):hiv +
         
         hiv:mri_wk + 
         mo(mrc):mri_wk +
         mo(genotype):mri_wk +
         mo(mrc):mo(genotype):mri_wk + 
         arm:mri_wk +
         cd4:mri_wk +
         
         arm:mo(mrc):mri_wk + 
         arm:mo(genotype):mri_wk + 
         arm:hiv:mri_wk + 
         arm:cd4:mri_wk +
         mo(mrc):hiv:mri_wk +
         mo(genotype):hiv:mri_wk +
         
         arm:mo(mrc):hiv:mri_wk +
         arm:mo(genotype):hiv:mri_wk +
         arm:mo(genotype):mo(mrc):mri_wk
         | gr(region, cor=FALSE, id=3)
      )
  )

rv_fml <-
  bf(
    vol | mi() ~ 0 +
      s(mri_wk, k=3)  +
      arm + 
      mo(genotype):hiv +
      mo(mrc):hiv +
      arm:hiv:mri_wk + 
      mo(mrc):hiv:mri_wk +
      arm:mo(genotype):mri_wk +
      arm:mo(mrc):hiv:mri_wk +
      (1|obj) + 
      (0 + mri_wk + arm + 
         mo(genotype):hiv +
         mo(mrc):hiv +
         arm:hiv:mri_wk + 
         mo(mrc):hiv:mri_wk +
         arm:mo(genotype):mri_wk +
         arm:mo(mrc):hiv:mri_wk 
       | gr(region, cor=FALSE, id=3)
       )
  )
  
