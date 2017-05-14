library(dplyr)
library(ggplot2)
library(magrittr)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## 
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
## 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

processed_post <- 
  bind_rows(processed_post) %>%
  group_by(demog.var, parameter)

par_fix_post <- filter(processed_post, !grepl("r_2_", parameter))
par_ran_post <- filter(processed_post,  grepl("r_2_", parameter))

par_fix <- 
  par_fix_post %>%
  select(-iteration, -chain) %>%
  summarise_if(is.numeric, mean)
  
par_fix$parameter



