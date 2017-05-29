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
  gr_lmb = "Lamb simplified growth spring corr model Posteriors.Rdata",
  gr_all = "Growth age square corr model Posteriors.Rdata",
  su_lmb = "Lamb simplified survival spring model Posteriors.Rdata",
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
## utlity functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# 
# inverse logit
# 
inv_logit <- function(x) {
  1.0/(1.0+exp(-x))
}

# 
# build an empty list in array form
# 
mk_list_array <- function(states) {
  dims <- sapply(states, length)
  data <- vector(mode = "list", length = prod(dims))
  array(data, dims, states)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## utlity functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# 
# assignment happens in the environment where the function is evaluated
# 
assign_required <- function(int_nod, mod_par) {
  # build the equivalent density function and assign
  expect_f <- list(expect_f = build_expect_f(int_nod$x))
  list2env(expect_f, envir = parent.frame())
  # expand the integration nodes and assign
  int_nod <- do.call(expand.grid, int_nod)
  list2env(int_nod, envir = parent.frame())
  # assign the model parameters
  mod_par <- as.list(mod_par)
  list2env(mod_par, envir = parent.frame())
}

# 
# calculate expectation required for equivalent density
# 
build_expect_f <- function(x) {
  dx <- x[2]-x[1]
  function(n_t, theta) sum(n_t * x^theta) * dx
}

# 
# winter survival function (sex by age)
# 
build_su <- function(int_nod, mod_par) {
  # store age set and remove from 'int_nod'
  age_set <- int_nod$a
  int_nod <- int_nod[names(int_nod) != "a"]
  # 
  assign_required(int_nod, mod_par)
  #
  states <- mk_list_array(list(c("f","m"), age_set))
  # 
  fixed       <- states
  fixed["f",] <- lapply(age_set, function(a) {
    b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2
  })
  fixed["m",] <- lapply(age_set, function(a) {
    b_0_m + b_z1_m * log(x) + b_a1_m * a + b_a2_m * a^2 + fixed[["f", a+1]]
  })
  # 
  function(i, n_t) {
    from <- states
    #
    env_eff <- yr_ef[i,1] + 
      (gamma + yr_ef[i,2]) * expect_f(n_t, theta) / x^theta
    # 
    from["f",] <- lapply(age_set, function(a) {
      inv_logit(env_eff + fixed[["f", a+1]])
    })
    from["m",] <- lapply(age_set, function(a) {
      inv_logit(env_eff + fixed[["m", a+1]])
    })
    function(s, a) from[[s, a+1]]
  }
}

# 
# winter growth function (sex by age)
# 
build_gr <- function(int_nod, mod_par) {
  # store age set and remove from 'int_nod'
  age_set <- int_nod$a
  int_nod <- int_nod[names(int_nod) != "a"]
  # 
  assign_required(int_nod, mod_par)
  #
  states <- mk_list_array(list(c("f","m"), age_set))
  #
  fixed       <- states
  fixed["f",] <- lapply(age_set, function(a) {
    b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2
  })
  fixed["m",] <- lapply(age_set, function(a) {
    b_0_m + b_z1_m * log(x) + b_a1_m * a + b_a2_m * a^2 + fixed[["f", a+1]]
  })
  # 
  function(i, n_t) {
    from <- states
    #
    env_eff <- yr_ef[i,1] + 
      (gamma + yr_ef[i,2]) * expect_f(n_t, theta) / x^theta
    # 
    from["f",] <- lapply(age_set, function(a) {
      dnorm(log(x1), mean = env_eff + fixed[["f", a+1]], sd = sigma) / x1
    })
    from["m",] <- lapply(age_set, function(a) {
      dnorm(log(x1), mean = env_eff + fixed[["m", a+1]], sd = sigma) / x1
    })
    function(s, a) from[[s, a+1]]
  }
}

# 
# reproduction function (female only by age)
# 
build_rp <- function(int_nod, mod_par) {
  # 
  assign_required(int_nod, mod_par)
  #
  fixed <- b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2
  # 
  function(i, n_t) {
    #
    env_eff <- yr_ef[i,1] + 
      (gamma + yr_ef[i,2]) * expect_f(n_t, theta) / x^theta
    # 
    inv_logit(env_eff + fixed)
  }
}

# 
# reproduction function (female only by age)
# 
build_tw <- function(int_nod, mod_par) {
  # 
  assign_required(int_nod, mod_par)
  #
  fixed <- b_0_f + b_z1_f * log(x) + b_a1_f * a + b_az_f * a * log(x)
  # 
  function(i, n_t) {
    #
    env_eff <- yr_ef[i,1] + 
      (gamma + yr_ef[i,2]) * expect_f(n_t, theta) / x^theta
    # 
    inv_logit(env_eff + fixed) * (a > 0) # <- lambs can't have twins
  }
}

# 
# lamb size function (female only)
# 
build_bw_lmb <- function(int_nod, mod_par) {
  # 
  assign_required(int_nod, mod_par)
  # calculate the terms that don't vary
  fixed <- b_0_f + is_m * b_0_m + is_tw * b_0_tw + 
           b_a1_f * a + b_z1_f * log(x) + b_az_f * a * log(x) + b_a2_f * a^2
  #
  function(i, n_t) {
    env_eff <- yr_ef[i,1] + 
      (gamma + yr_ef[i,2]) * expect_f(n_t, theta) / x^theta
    dnorm(log(x1), fixed + env_eff, sigma) / x1
  }
}

build_su_lmb <- function(int_nod, mod_par) {
  # 
  assign_required(int_nod, mod_par)
  # FIXME!
  a <- 2
  # calculate the terms that don't vary
  fixed <- b_0_f + is_m * b_0_m + is_tw * b_0_tw + 
           b_z1_lm * log(x) + b_a1_f * a + b_a2_f * a^2
           
  #
  function(i, n_t) {
    env_eff <- yr_ef[i,1]
    inv_logit(fixed + env_eff)
  }
}

build_gr_lmb <- function(int_nod, mod_par) {
  # 
  assign_required(int_nod, mod_par)
  # FIXME!
  a <- 2
  # calculate the terms that don't vary
  fixed <- b_0_f + is_m * b_0_m + is_tw * b_0_tw + 
           b_z1_lm * log(x) + b_a1_f * a + b_a2_f * a^2
    
  #
  function(i, n_t) {
    env_eff <- yr_ef[i,1]
    dnorm(log(x1), fixed + env_eff, sigma) / x1
  }
}


build_get_age <- function(n_node) {
  force(n_node)
  function(n_t, a) nt[seq.int(1, n_node)+n_node*a] 
}

n_node <- 50
n_ages <- 15 
a_set <- seq.int(0, n_ages-1)
x_vals <- seq(5, 40,  length = n_node)
x_lmb_vals <- seq(0.1, 4, length = n_node)

get_age <- build_get_age(n_node)

dummy_yr <- list(yr_ef = matrix(0, ncol = 2))
mp_su <- c(par_fix$survival, dummy_yr)
mp_gr <- c(par_fix$growth,   dummy_yr)

node_sg <- list(x1 = x_vals, x = x_vals, a = a_set)
su <- build_su(node_sg, mp_su)
gr <- build_su(node_sg, mp_gr)

su_t <- su(1, nt)
gr_t <- gr(1, nt)

a <- 1
{su_t("f", a) * gr_t("f", a)} %*% get_age(n_t$f, a)




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Testing
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

x_vals <- seq(5, 40, length = 200)
x_lmb_vals <- seq(0.1, 4, length = 200)

nt <- 200 * dnorm(x_vals, mean = 12, sd = 1.5) +
      300 * dnorm(x_vals, mean = 22, sd = 3.0)
plot(x_vals, nt)

i_age <- function(a) {
  seq_along(x_vals)+length(x_vals)*a
}

# survival
mp <- c(par_fix$survival, list(yr_ef = matrix(0, ncol = 2)))
nd <- list(x = x_vals, a = 0:20)
f1 <- build_su(nd, mp)
f2 <- f1(1, nt)
plot(x_vals, f2("m", 0), ylim = c(0,1))

# reproduction
mp <- c(par_fix$reproduction, list(yr_ef = matrix(0, ncol = 2)))
nd <- list(x = x_vals, a = 0:20)
f1 <- build_rp(nd, mp)
f2 <- f1(1, nt)
plot(x_vals, f2[i_age(1)], ylim = c(0,1))

# twinning
mp <- c(par_fix$twinning, list(yr_ef = matrix(0, ncol = 2)))
nd <- list(x = x_vals, a = 0:20)
f1 <- build_tw(nd, mp)
f2 <- f1(1, nt)
plot(x_vals, f2[i_age(1)], ylim = c(0,1))

# growth
mp <- c(par_fix$growth, list(yr_ef = matrix(0, ncol = 2)))
nd <- list(x1 = x_vals, x = x_vals, a = 0:20)
f1 <- build_gr(nd, mp)
f2 <- f1(1, nt)
gk <- f2("m", 1)
dim(gk) <- rep(length(x_vals), 2)
image(x_vals, x_vals, gk, col = heat.colors(512))

# lamb birth weight
mp <- c(par_fix$lamb.birth, list(yr_ef = matrix(0, ncol = 2)))
nd <- list(is_m = c(0,1), is_tw = c(0,1), 
           x1 = x_lmb_vals, x = x_vals, a = 0:20)
f1 <- build_bw_lmb(nd, mp)
f2 <- f1(1, nt)
dim(f2) <- sapply(nd, length)
gk <- f2[1,1,,,5]
image(x_vals, x_lmb_vals, t(gk), col = heat.colors(512))

# lamb survival
mp <- c(par_fix$lamb.survival, list(yr_ef = matrix(0, ncol = 2)))
nd <- list(is_m = c(0,1), is_tw = c(0,1), x = x_lmb_vals)
f1 <- build_su_lmb(nd, mp)
f2 <- f1(1, nt)
dim(f2) <- sapply(nd, length)
plot(x_lmb_vals, f2[1,1,], ylim = c(0,1))

# lamb spring growth
mp <- c(par_fix$lamb.growth, list(yr_ef = matrix(0, ncol = 2)))
nd <- list(is_m = c(0,1), is_tw = c(0,1), x1 = x_vals, x = x_lmb_vals)
f1 <- build_gr_lmb(nd, mp)
f2 <- f1(1, nt)
dim(f2) <- sapply(nd, length)
gk <- f2[1,1,,]
image(x_lmb_vals, x_vals, t(gk), col = heat.colors(512))












