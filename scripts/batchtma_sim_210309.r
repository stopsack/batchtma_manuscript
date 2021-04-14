# Plasmode simulation study for batch effects
# Konrad H. Stopsack <stopsack@mskcc.org>, Nov 11, 2020 onwards
#

library(here)       # relative path to project directory
library(khsverse)   # v0.0.1; includes khsmisc v0.3.1; remotes::install_github("stopsack/khsverse")
library(batchtma)   # v0.1.3; remotes::install_github("stopsack/batchtma")
library(simsurv)    # v1.0.0; simulated survival data
library(flexsurv)   # v2.0; flexible parametric survival models
library(rsample)    # tidyverse bootstrapping
library(tictoc)     # benchmarking
library(furrr)      # future_map2() when simulating eventtimes; updated to v0.2.1

comb <- readRDS(file = here("out", "for_simulation_210309.rds"))  # saved by batchtma_sim_210309.Rmd

############################################################################################################
# All function definitions
############################################################################################################
# 3) Simulate event times 

# Requirements:
# a) Exposure must be called "exposure" 
# b) Exposure must not be missing
# c) Confounder changes need to be hard-coded for computational efficiency
 
logcumhaz <- function(t, x, betas, knots) {  # Calculate the log cumulative hazard
  basis <- flexsurv::basis(knots, log(t))
  betas[["gamma0"]]              * basis[[1]] + 
    betas[["gamma1"]]            * basis[[2]] +
    betas[["gamma2"]]            * basis[[3]] +
    betas[["gamma3"]]            * basis[[4]] +
    betas[["gamma4"]]            * basis[[5]] +
    betas[["exposure"]]          * x[["exposure"]] +
    betas[["gleason3+4"]]        * if_else(x[["gleason"]]    == "3+4",     1, 0) +
    betas[["gleason4+3"]]        * if_else(x[["gleason"]]    == "4+3",     1, 0) +
    betas[["gleason8"]]          * if_else(x[["gleason"]]    == "8",       1, 0) +
    betas[["gleason9-10"]]       * if_else(x[["gleason"]]    == "9-10",    1, 0) +
    betas[["ptnm_miss2pT3"]]     * if_else(x[["ptnm_miss2"]] == "pT3",     1, 0) +
    betas[["ptnm_miss2hirisk"]]  * if_else(x[["ptnm_miss2"]] == "hirisk",  1, 0)
}

eventtime_sim <- function(split, truehr, truemodel) {
  coef_edited <- truemodel$coefficients
  coef_edited[["exposure"]] <- log(truehr)
  
  dat <- simsurv(betas        = coef_edited,      # "true" parameter values
                 x            = analysis(split),  # covariate data
                 knots        = truemodel$knots,  # knot locations for splines
                 logcumhazard = logcumhaz,        # definition of log cum hazard
                 idvar        = "id",             # critical with bootstrapping
                 ids          = unique(analysis(split)$id),
                 maxt         = 300,              # no right-censoring = NULL
                 interval     = c(1E-6,1000000))  # interval for root finding
  inner_join(analysis(split), dat, by = "id")
}

poss_eventtime_sim <- possibly(
  .f = eventtime_sim, 
  otherwise = tibble(id = NA, eventlethal = NA, lethalfu = NA, 
                     gleason = NA, ptnm_miss2 = NA, exposure = NA, 
                     eventtime = NA, status = NA))

simulate_eventtimes <- function(realdata, n_hr = 10, replications = 10) {
  realdata <- realdata %>% 
    select(id, eventlethal, lethalfu, gleason, ptnm_miss2, exposure)
  
  true_mod <- realdata %>%  # estimate true survival model
    flexsurv::flexsurvspline(  
      formula = Surv(time = lethalfu, event = eventlethal) ~ exposure + gleason + ptnm_miss2,
      k = 3)
  
  tic()
  
  plan(multisession, workers = 4)  # for furrr::future_map()
  dat <- bootstraps(data = realdata, times = n_hr * replications) %>%
    mutate(truehr = rep(exp(seq(from = log(1/3), to = log(3), length.out = n_hr)), 
                        times = replications),
           data   = future_map2(.x = splits, .y = truehr, 
                                .f = poss_eventtime_sim, truemodel = true_mod,
                                .options = furrr_options(seed = TRUE)),
           isna   = map_lgl(.x = data, .f = ~is.na(.$eventlethal[1]), TRUE, FALSE)) %>%
    filter(!isna) %>% 
    select(-isna)
  toc()
  dat
}

############################################################################################################
# 4) Add batch effects
# batch means/variances as per real data
add_batcheffects <- function(realdata) {
  tma_means_er_alpha <- c("A"  = 1.23,
                          "B"  = 0.71,
                          "C"  = 0.25,
                          "D" = 0.69,
                          "E" = 0.54, 
                          "F" = 0.4, 
                          "G" = -0.06, 
                          "H" =	-0.71, 
                          "I" = -0.34,
                          "J" = -0.74,
                          "K" = -0.98, 
                          "L" = -0.35, 
                          "M" = -0.14, 
                          "N" = -0.5)
  tma_vars_er_alpha  <- c("A"  = 1.48, 
                          "B"  = 1.1,
                          "C"  = 1.02,
                          "D" = 1.16,
                          "E" = 2.31, 
                          "F" = 1.22, 
                          "G" = 1.09, 
                          "H" = 0.25,
                          "I" = 0.7,
                          "J" = 0.45,
                          "K" = 0.57, 
                          "L" = 0.67, 
                          "M" = 1.09, 
                          "N" = 0.89)
  
  realdata %>% 
    # (B) means only, like for ER-alpha
    mutate(exposureBbtc = exposure + tma_means_er_alpha[as.character(tma)]) %>%  
    group_by(tma) %>%
    # (C) means and variance, like for ER-alpha
    mutate(exposureCbtc = (exposure - mean(exposure)) * sqrt(tma_vars_er_alpha[as.character(tma)]) +
             mean(exposure) + tma_means_er_alpha[as.character(tma)],
           # (D) means, variance like for ER-alpha plus effect modification by Gleason
           exposureDbtc = (exposure - mean(exposure)) * 
             sqrt(tma_vars_er_alpha[as.character(tma)]) * 
             (as.numeric(gleason) - mean(as.numeric(gleason)) + 6)/6 +
             mean(exposure) + 
             tma_means_er_alpha[as.character(tma)] * (as.numeric(gleason) - 
                                                        mean(as.numeric(gleason)) + 6)/6) %>%
    ungroup() %>%
    all_adjustments(markers = "exposure") %>%
    all_adjustments(markers = "exposureBbtc") %>%
    all_adjustments(markers = "exposureCbtc") %>%
    all_adjustments(markers = "exposureDbtc")
}

############################################################################################################
# 5) Generate batch effect-corrected values

all_adjustments <- function(data, markers, drop = FALSE) {
  res <- data %>% 
    adjust_batch(markers = all_of(markers), batch = tma, method = simple) %>%
    adjust_batch(markers = all_of(markers), batch = tma, method = standardize,
                confounders = c(gleason, ptnm_miss2)) %>%
    adjust_batch(markers = all_of(markers), batch = tma, method = ipw,
                confounders = c(gleason, ptnm_miss2)) %>%
    adjust_batch(markers = all_of(markers), batch = tma, method = quantreg,
                confounders = c(gleason, ptnm_miss2), quantreg_tau = c(0.25, 0.75)) %>%
    adjust_batch(markers = all_of(markers), batch = tma, method = quantnorm)
  if(drop == TRUE) {
    res <- res %>% 
      select(id, starts_with(markers), tma)
  }
  res
}

############################################################################################################
# 6) Evaluate the final simulated data
cox_adj <- function(data, adjust, multivariable = TRUE, stratified = FALSE) {
  coxph(formula = as.formula(paste0("Surv(time = eventtime, event = status) ~ exposure", 
                                    adjust, 
                                    if_else(multivariable == TRUE, " + gleason + ptnm_miss2", ""),
                                    if_else(stratified    == TRUE, " + strata(tma)",          ""))),
        data = data) %>%
    tidy(exponentiate = TRUE) %>%
    filter(term == paste0("exposure", adjust)) %>% 
    pull(estimate) %>%
    as.numeric()
}

cox_stratified_ivpool <- function(data, adjust, multivariable = TRUE) {
  data %>%
    distinct(tma) %>%
    mutate(res = map(
      .x = tma,
      .f = ~{ 
        coxph(formula = as.formula(paste0("Surv(time = eventtime, event = status) ~ exposure", adjust, 
                                          if_else(multivariable == TRUE, " + gleason + ptnm_miss2", ""))),
              data = data %>% filter(tma == .x)) %>%
          tidy() %>% 
          filter(term == paste0("exposure", adjust)) })) %>%
    unnest_wider(col = res) %>%
    mutate(weight = 1 / (std.error ^ 2)) %>%  # note, na.rm = FALSE by default
    summarize(estimate  = exp(sum(estimate * weight) / sum(weight))) %>%
    pull(estimate)  # will be missing if any batch-specific estimates are missing!
}

poss_cox_stratified_pool <- possibly(.f = cox_stratified_ivpool, otherwise = NA_real_)

gen_adj_hrs <- function(simdata, expdata, suffix, batcheffect, multivariable = TRUE) {
  expdata <- expdata %>% add_batcheffects()
  
  plan(multisession, workers = 4)  # for furrr::future_map()
  
  simdata %>% 
    mutate(
      data_adj = map(.x = data,
                     .f = ~left_join(.x %>% select(-gleason, -ptnm_miss2, -exposure), 
                                     expdata %>% select(-eventlethal, -lethalfu), by = "id")),
      hr_adj1 = map2_dbl(.x = data_adj, .y = paste0(suffix, ""),      
                         .f = cox_adj, multivariable = multivariable),
      hr_adj2 = map2_dbl(.x = data_adj, .y = paste0(suffix, "_adj2"), 
                         .f = cox_adj, multivariable = multivariable),
      hr_adj3 = map2_dbl(.x = data_adj, .y = paste0(suffix, "_adj3"), 
                         .f = cox_adj, multivariable = multivariable),
      hr_adj4 = map2_dbl(.x = data_adj, .y = paste0(suffix, "_adj4"), 
                         .f = cox_adj, multivariable = multivariable),
      hr_adj5 = map2_dbl(.x = data_adj, .y = paste0(suffix, "_adj5"), 
                         .f = cox_adj, multivariable = multivariable),
      hr_adj6 = map2_dbl(.x = data_adj, .y = paste0(suffix, "_adj6"), 
                         .f = cox_adj, multivariable = multivariable),
      hr_adj8 = map2_dbl(.x = data_adj, .y = suffix, 
                         .f = poss_cox_stratified_pool, multivariable = multivariable),
      hr_adj9 = map2_dbl(.x = data_adj, .y = suffix, 
                         .f = cox_adj, multivariable = multivariable,
                         stratified = TRUE)) %>%
    select(truehr, starts_with("hr_adj")) %>%
    mutate(batcheff = batcheffect)
}

all_adj_hrs <- function(realdata, simdata, datalabel) {
  bind_rows(
    gen_adj_hrs(simdata = simdata, expdata = realdata, suffix = "",     batcheffect = 1),
    gen_adj_hrs(simdata = simdata, expdata = realdata, suffix = "Bbtc", batcheffect = 2),
    gen_adj_hrs(simdata = simdata, expdata = realdata, suffix = "Cbtc", batcheffect = 3),
    gen_adj_hrs(simdata = simdata, expdata = realdata, suffix = "Dbtc", batcheffect = 4)) %>%
    mutate(datalabel = datalabel)
}

############################################################################################################
# Run everything with data: 
#  1) Load data
#  2) Induce confounding
#  3) Simulate event times
#  4) Add batch effects
#  5) Remove batch effects
#  6) Evaluate everything
############################################################################################################
set.seed(123)
comb <- comb %>% 
  filter(lethalfu != 0) %>%
  mutate(exposure   = rnorm(n = nrow(.)),
         # must collapse stage pT3 because of non-convergence of flexsurvreg() in dfexp4
         ptnm_miss2 = fct_collapse(ptnm_miss2, pT3 = c("pT3/T3a", "pT3b")))

# (1) Observed biomarker
dfexp1 <- comb

# (2) Biomarker forced to be associated with calendar year, Gleason and stage (modestly)
dfexp2 <- comb %>%
  mutate(exposure = exposure + 
           0.2 * (as.numeric(gleason)    - mean(as.numeric(gleason))) +
           0.1 * (as.numeric(ptnm_miss2) - mean(as.numeric(ptnm_miss2))))

# # (3) Biomarker forced to be associated with calendar year, Gleason and stage (strongly)
dfexp3 <- comb %>%
  mutate(exposure = exposure +
           0.4 * (as.numeric(gleason)    - mean(as.numeric(gleason))) +
           0.2 * (as.numeric(ptnm_miss2) - mean(as.numeric(ptnm_miss2))))

# (4) Strong confounding plus stronger imbalance of confounders across batches
#     a)  Alter biomarker as in (3), plus:
#     b)  Exclude tumors with Gleason >= 8 on TMAs 49, 124, and 125
#     c)  Exclude tumors with stage >= pT3b/M1/missing on TMAs 127-318
dfexp4 <- comb %>%
  filter(!(ptnm_miss2 %in% c("hirisk", "M1", "pT4/N1") & tma %in% c("H", "K", "L", "M", "N"))) %>%
  filter(!(gleason %in% c("5-6") & tma %in% c("A", "E", "F"))) %>%
  mutate(exposure = exposure +
           0.4 * (as.numeric(gleason)    - mean(as.numeric(gleason))) +
           0.2 * (as.numeric(ptnm_miss2) - mean(as.numeric(ptnm_miss2))))

# Run simulations
n_hr <- 7  # 9 steps between HR = 1/3 and HR = 3, including 1.0
replications <- 200
# dfsim1 <- simulate_eventtimes(realdata = dfexp1, n_hr = n_hr, replications = replications)
# dfsim2 <- simulate_eventtimes(realdata = dfexp2, n_hr = n_hr, replications = replications)
# dfsim3 <- simulate_eventtimes(realdata = dfexp3, n_hr = n_hr, replications = replications)
# dfsim4 <- simulate_eventtimes(realdata = dfexp4, n_hr = n_hr, replications = replications)

# save_safely(dfsim1, dfsim2, dfsim3, dfsim4,
#             file = here("out", "sim_data_df1234_500_210309.RData"))
# saveRDS_safely(dfsim2, file = here("out", "sim_data_df2_500_210309.rds"))

load(here("out", "sim_data_df1234_500_210309.RData"))

tic()
res1 <- all_adj_hrs(realdata = dfexp1, simdata = dfsim1, datalabel = 1)
toc()
res2 <- all_adj_hrs(realdata = dfexp2, simdata = dfsim2, datalabel = 2)
toc()
res3 <- all_adj_hrs(realdata = dfexp3, simdata = dfsim3, datalabel = 3)
toc()
res4 <- all_adj_hrs(realdata = dfexp4, simdata = dfsim4, datalabel = 4)
toc()

bind_rows(res1, res2, res3, res4) %>%
  write_csv_safely(file = here("out", "sim_allhrs_500_210309.csv"))

# Environment
sessionInfo()

# Clean up!
rm(list = ls())