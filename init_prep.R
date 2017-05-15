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

unique(processed_post$p)

par_fix_post <- filter(processed_post, !grepl("r_2_", p))
par_ran_post <- filter(processed_post,  grepl("r_2_", p))

unique(processed_post$demog.var)

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



n_funs <- length(file_names)
build_demog_funs <- vector(mode   = "list", length = n_funs)
names(build_demog_funs) <- names(file_names)
rm(n_funs)

build_demog_funs$bw_lmb <- function(p, x1, x) {
  
}
  













