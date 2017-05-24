library(dplyr)
library(ggplot2)
library(magrittr)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Read in / process the posterior samples
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

par_location <- "~/Dropbox/Shared Folders/Asymm_soay/results/Posteriors"

file_names <- c(
  bw_lmb = "Lamb Birth weight corr model Posteriors.Rdata",
  gr_lmb = "Lamb growth spring corr model Posteriors.Rdata",
  gr_all = "Growth age square corr model Posteriors.Rdata",
  su_lmb = "Lamb survival spring corr model Posteriors.Rdata",
  su_old = "Survival  age square corr model Posteriors.Rdata",
  rp_all = "Reproduction corr model Posteriors.Rdata",
  tw_old = "Twinning corr model Posteriors.Rdata"
)

processed_post <- list()
for (i in seq_along(file_names)) {
  # 
  load(file.path(par_location, file_names[i]))
  #   
  processed_post[[i]] <- 
    tbl_df(post) %>% 
    filter(
      # remove individual intercept deviations
      !grepl("r_1_", parameter),
      # remove the variance component parameters
      !grepl("sd_" , parameter),
      # remove the posterior log p
      !grepl("lp__", parameter)
    ) %>%
    mutate(
      iteration = as.integer(levels(iteration))[iteration],
      chain     = as.integer(levels(chain))[chain],
      parameter =  as.character(parameter)
    ) %T>% print()
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Pull out samples from the posterior
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

processed_post <- 
  bind_rows(processed_post) %>% 
  rename(p = parameter) %>%
  mutate(
    p = ifelse(p == "intercept",     "b_0_f",  p),
    p = ifelse(p == "Intercept",     "b_0_f",  p),
    p = ifelse(p == "sexM",          "b_0_m",  p),
    p = ifelse(p == "z",             "b_z1_f", p),
    p = ifelse(p == "z_sexM",        "b_z1_m", p),
    p = ifelse(p == "age",           "b_a1_f", p),
    p = ifelse(p == "sexM_age",      "b_a1_m", p),
    p = ifelse(p == "IageE2",        "b_a2_f", p),
    p = ifelse(p == "z_age",         "b_az_f", p),
    p = ifelse(p == "sexM_IageE2",   "b_a2_m", p),
    p = ifelse(p == "b_equivdens_1", "gamma",  p),
    p = ifelse(p == "sex.lambM",     "b_0_m",  p),
    p = ifelse(p == "offaretwins",   "b_0_tw", p),
    p = ifelse(p == "z_birth.lamb",  "b_z1_lm", p)
  ) %>%
  group_by(demog.var, p)

par_fix_post <- filter(processed_post, !grepl("r_2_", p))
par_ran_post <- filter(processed_post,  grepl("r_2_", p))

unique(processed_post$demog.var)
unique(filter(processed_post, demog.var == "survival")$p)

## quick look at a specific posterior distribution
processed_post %>% 
  filter(
    demog.var == "lamb.growth",
    p == "b_a1_f"
  ) %>%
  ggplot(aes(x = value)) + 
  geom_histogram() + 
  geom_vline(xintercept = 0, col = "red")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

par_fix <- 
  par_fix_post %>%
  select(-iteration, -chain) %>%
  summarise_if(is.numeric, mean)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

build_su <- function(int_nod, mod_par) {
  # assign the model parameters + integration nodes
  list2env(int_nod, envir = sys.frame(1))
  list2env(mod_par, envir = sys.frame(1))
  # construct equivalent density function
  calc_equi_dens <- build_equi_dens(x, theta)
  # calculate the terms that don't vary
  fixed_f <- b_0_f + b_z1_f * log(x) + 
  fixed_m <- b_0_m + b_z1_m * log(x) + fixed_f
  # 
  fixed_f_by_age <- vector(mode = "list", length = n_ages)
  fixed_m_by_age <- vector(mode = "list", length = n_ages)
  #
  for (i in seq_len(n_ages)) {
    fixed_f_by_age[[i]] <- 
      fixed_f + b_a1_f * a + b_a2_f * a^2
    fixed_m_by_age[[i]] <- 
      fixed_m + b_a1_m * a + b_a2_m * a^2 + fixed_f_by_age[[i]]
  }
  # 
  function(i, n_t) {
    #
    env_eff <- yr_ef[i,1] + (gamma + yr_ef[i,2]) * calc_equi_dens(n_t)
    # 
    age <- vector(mode = "list", length = n_ages)
    #
    for (i in seq_len(n_ages)) {
      age[[i]] <- function() {
        
      }
    }
  }
  
}


build_bw_lmb <- function(nodes, mod_par) {
  # assign the model parameters + integration nodes
  list2env(nodes,   envir = sys.frame(1))
  list2env(mod_par, envir = sys.frame(1))
  # 
  calc_equi_dens <- build_equi_dens(x, theta)
  # calculate the terms that don't vary
  fixed <- b_0_f + b_a1_f * a + b_z1_f * log(x) + b_az_f * a * log(x) + b_a2_f * a^2
  # 
  fixed_s_f <- log(fixed)
  fixed_s_m <- log(fixed + b_0_m )
  fixed_t_f <- log(fixed + b_0_tw)
  fixed_t_m <- log(fixed + b_0_m + b_0_tw)
  #
  function(i, n_t) {
    #
    env_eff <- yr_ef[i,1] + (gamma + yr_ef[i,2]) * calc_equi_dens(n_t)
    #
    to_state <- vector(mode = "list", length = 4)
    # FIX ME!
    to_state$s_f <- function () { # singleton females
      dnorm(x1, mean = fixed_s_f + env_eff, sd = sigma)
    }
    to_state$s_m <- function () { # singleton males
      dnorm(x1, mean = fixed_s_m + env_eff, sd = sigma)
    }
    to_state$t_f <- function () { # twin females
      dnorm(x1, mean = fixed_t_f + env_eff, sd = sigma)
    }
    to_state$t_m <- function () { # twin males
      dnorm(x1, mean = fixed_t_m + env_eff, sd = sigma)
    }
  }
}



tmp <- list(a = 10, b = 2)
f_test <- function(pars) {
  list2env(tmp, envir = sys.frame(1))
  function(n) n * a
}
f_test(tmp)(4)
a
b
rm(a, b)



build_demog_funs$bw_lmb_f1 <- function(x1, x, a, m, p) {
  mu_part <- p["b_0_f"] + p["b_0_tw"] + p["b_a1_f"]*a + p["b_a1_f"]*a^2 + p["b_az_f"]*a*log(x)
  function() {
    
  }
}

build_demog_funs$bw_lmb_f1 <- function(x1, x, a, m, p) {
  mu_part <- p["b_0_f"] + p["b_0_tw"] + p["b_a1_f"]*a + p["b_a1_f"]*a^2 + p["b_az_f"]*a*log(x)
  function() {
    
  }
}


build_demog_funs$bw_lmb_m <- function(x1, x, a, m, p) {
  mu_part <- p["b_0_f"] + p["b_0_m"] + p["b_0_tw"] + 
  function() {
      
  }
}

  













