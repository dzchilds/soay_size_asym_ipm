library(dplyr)
library(ggplot2)
library(magrittr)

options(stringsAsFactors = FALSE)

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
  summarise_if(is.numeric, mean) %>% 
  split(.$demog.var) %>%
  lapply(function(df) {
    pars <- df$value
    names(pars) <- df$p
    as.list(pars)
  })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# 
# assignment happens in the environment where the function is evaluated
# 
assign_required <- function(int_nod, mod_par) {
  # build the equivalent density function and assign
  equ_den <- list(equ_den = build_equ_den(int_nod$x))
  list2env(equ_den, envir = parent.frame())
  # expand the integration nodes and assign
  int_nod <- do.call(expand.grid, int_nod)
  list2env(int_nod, envir = parent.frame())
  # assign the model parameters
  mod_par <- as.list(mod_par)
  list2env(mod_par, envir = parent.frame())
}

mk_list_array <- function(states) {
  dims <- sapply(states, length)
  data <- vector(mode = "list", length = prod(dims))
  array(data, dims, states)
}

build_equ_den <- function(x) {
  dx <- x[2]-x[1]
  function(n_t, theta) {
    sum(n_t * x^theta) * dx / x^theta
  }
}

build_su <- function(int_nod, mod_par) {
  # store age set and remove from 'int_nod'
  age_set <- int_nod$a
  int_nod <- int_nod[names(int_nod) != "a"]
  # 
  assign_required(int_nod, mod_par)
  #
  fixed_f <- lapply(age_set, function(a) {
    b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2
  })
  fixed_m <- lapply(age_set, function(a) {
    b_0_m + b_z1_m * log(x) + b_a1_m * a + b_a2_m * a^2 + fixed_f[[a+1]]
  })
  # 
  function(i, n_t) {
    #
    env_eff <- yr_ef[i,1] + (gamma + yr_ef[i,2]) * equ_den(n_t, theta)
    # 
    from_sex   <- list()
    from_sex$f <- lapply(age_set, function(a) {
      1/(1+exp(-(env_eff + fixed_f[[a+1]])))
    })
    from_sex$m <- lapply(age_set, function(a) {
      1/(1+exp(-(env_eff + fixed_m[[a+1]])))
    })
    function(s, a) from_sex[[s]][[a+1]]
  }
}

build_gr <- function(int_nod, mod_par) {
  # store age set and remove from 'int_nod'
  age_set <- int_nod$a
  int_nod <- int_nod[names(int_nod) != "a"]
  # 
  assign_required(int_nod, mod_par)
  #
  fixed_f <- lapply(age_set, function(a) {
    b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2
  })
  fixed_m <- lapply(age_set, function(a) {
    b_0_m + b_z1_m * log(x) + b_a1_m * a + b_a2_m * a^2 + fixed_f[[a+1]]
  })
  # 
  function(i, n_t) {
    #
    env_eff <- yr_ef[i,1] + (gamma + yr_ef[i,2]) * equ_den(n_t, theta)
    # 
    from_sex   <- list()
    from_sex$f <- lapply(age_set, function(a) {
      dnorm(log(x1), mean = env_eff + fixed_f[[a+1]], sd = sigma) / x1
    })
    from_sex$m <- lapply(age_set, function(a) {
      dnorm(log(x1), mean = env_eff + fixed_m[[a+1]], sd = sigma) / x1
    })
    function(s, a) from_sex[[s]][[a+1]]
  }
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Testing
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

x_vals <- seq(5, 40, length = 250)
nt <- 500 * dnorm(x_vals, mean = 12, sd = 2)
# survival
mp <- c(par_fix$survival, list(yr_ef = matrix(0, ncol = 2)))
nd <- list(x = x_vals, a = 0:20)
f1 <- build_su(nd, mp)
f2 <- f1(1, nt)
plot(x_vals, f2("f", 10), ylim = c(0,1))
# growth
mp <- c(par_fix$growth, list(yr_ef = matrix(0, ncol = 2)))
nd <- list(x1 = x_vals, x = x_vals, a = 0:20)
f1 <- build_gr(nd, mp)
f2 <- f1(1, nt)
gk <- f2("f", 1)
dim(gk) <- rep(length(x_vals), 2)
image(x_vals, x_vals, gk, col = heat.colors(512))




build_bw_lmb <- function(int_nod, mod_par) {
  # 
  assign_required(int_nod, mod_par)
  #
  to_states <- list(c("f","m"), c("s","t"))
  # calculate the terms that don't vary
  ref <- b_0_f + b_a1_f * a + b_z1_f * log(x) + 
                   b_az_f * a * log(x) + b_a2_f * a^2
  # 
  fixed <- mk_list_array(to_states)
  fixed[["f","s"]] <- ref
  fixed[["m","s"]] <- ref + b_0_m
  fixed[["f","t"]] <- ref + b_0_tw
  fixed[["m","t"]] <- ref + b_0_m + b_0_tw
  #
  function(i, n_t) {
    env_eff <- yr_ef[i,1] + (gamma + yr_ef[i,2]) * equ_den(n_t, theta)
    to <- mk_list_array(to_states)
    to[["f","s"]] <- dnorm(log(x1), fixed[["f","s"]] + env_eff, sigma) / x1
    to[["m","s"]] <- dnorm(log(x1), fixed[["m","s"]] + env_eff, sigma) / x1
    to[["f","t"]] <- dnorm(log(x1), fixed[["f","t"]] + env_eff, sigma) / x1 
    to[["m","t"]] <- dnorm(log(x1), fixed[["m","t"]] + env_eff, sigma) / x1
    function(s, n) to_state[[s, n]]
  }
}


to_state <- vector(mode = "list", length = 4)
dim(to_state) <- c(2, 2)
dimnames(to_state) <- list(c("f","m"), c("s","t"))

states <- list(c("f","m"), c("s","t"))


mk_list_array(states)

lapply(c("f","m"), function(nm) nm)


