library(ggplot2)
library(dplyr)
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
## Relabelling
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# lookup table for demographic function names
c(survival      = "su",
  growth        = "gr",
  reproduction  = "rp",
  twinning      = "tw",
  lamb.birth    = "bw_lmb", 
  lamb.growth   = "gr_lmb", 
  lamb.survival = "su_lmb") -> func_labels

# lookup table for parameter names
c(intercept     = "b_0_f",
  Intercept     = "b_0_f",
  sexM          = "b_0_m",
  sex.lambM     = "b_0_m",
  z             = "b_z1_f",
  z_sexM        = "b_z1_m",
  age           = "b_a1_f",
  sexM_age      = "b_a1_m",
  IageE2        = "b_a2_f",
  z_age         = "b_az_f",
  sexM_IageE2   = "b_a2_m",
  offaretwins   = "b_0_tw",
  z_birth.lamb  = "b_z1_lm",
  b_equivdens_1 = "gamma",
  theta         = "theta",
  sigma         = "sigma") -> parm_labels

processed_post <- 
  bind_rows(processed_post) %>% 
  rename(p = parameter, f = demog.var) 

par_fix_post <- 
  filter(processed_post, !grepl("r_2_", p)) %>%
  mutate(f = func_labels[f], p = parm_labels[p]) %>%
  group_by(f, p)

par_ran_post <- 
  filter(processed_post,  grepl("r_2_", p)) %>%
  mutate(f = func_labels[f]) %>%
  group_by(f, p)
  
unique(par_fix_post$f)
unique(par_fix_post$p)
unique(par_ran_post$f)
unique(par_ran_post$p)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

dummy_yr <- list(yr_ef = matrix(0, ncol = 2))

par_fix <- 
  par_fix_post %>%
  select(-iteration, -chain) %>%
  summarise_if(is.numeric, mean) %>% 
  split(.$f) %>%
  lapply(function(df) {
    pars <- df$value
    names(pars) <- df$p
    c(as.list(pars), dummy_yr)
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## utlity functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

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
      dnorm(log(x_), mean = env_eff + fixed[["f", a+1]], sd = sigma) / x_
    })
    from["m",] <- lapply(age_set, function(a) {
      dnorm(log(x_), mean = env_eff + fixed[["m", a+1]], sd = sigma) / x_
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
    dnorm(log(x_), fixed + env_eff, sigma) / x_
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
    dnorm(log(x_), fixed + env_eff, sigma) / x_
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## utlity functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

numeric_midpoint <- function(x) {
  weights <- (x[2] - x[1]) / x[3]
  nodes <- seq(x[1] + weights/2, x[2] - weights/2, length.out = x[3])
  rule <- list(nodes, weights, x[3])
  names(rule) <- c("val", "wgt", "num")
  attr(rule, "state") <- "numeric"
  attr(rule, "rule") <- "midpoint"
  rule
}

numeric_integer <- function() {
  rule <- list(seq.int(x[1], x[2]), 1, x[2]-x[1]+1)
  names(rule) <- c("val", "wgt", "num")
  attr(rule, "state") <- "numeric"
  attr(rule, "rule") <- "integer"
  rule
}

categorical_integer <- function(x) {
  rule <- list(seq.int(x[1], x[2]), 1, x[2]-x[1]+1)
  names(rule) <- c("val", "wgt", "num")
  attr(rule, "state") <- "categorical"
  attr(rule, "rule") <- "integer"
  rule
}

categorical_nominal <- function(x) {
  rule <- list(x, 1, length(x))
  names(rule) <- c("val", "wgt", "num")
  attr(rule, "state") <- "categorical"
  attr(rule, "rule") <- "nominal"
  rule
}

which_list_part <- function (states) {
  type_info <- sapply(states, function(x) attributes(x)$state)
  type_info == "categorical"
}

new_n_t <- function(...) {
  #
  dots <- list(...) 
  # work out which states are in list form
  is_list_part <- which_list_part(dots)
  is_list_part
  # build the list part of the state vector
  list_part <- dots[is_list_part]
  n_t <- mk_list_array(lapply(list_part, '[[', "val"))
  attrib_n_t <- attributes(n_t)
  # insert the numeric state vector
  vect_part <- dots[!is_list_part]
  size <- prod(sapply(vect_part, '[[', "num"))
  zeros <- vector(mode = "double", size)
  n_t <- lapply(n_t, function(dummy) zeros)
  # make sure the numeric states are 1st in the attributes
  ord <- c(which(!is_list_part), which(is_list_part))
  dots <- dots[ord]
  # 
  state_names <- list(state_names = names(dots))
  #
  state_rules <- sapply(dots, function(x) attributes(x)[c("state","rule")])
  state_rules <- list(state_rules = state_rules)
  # 
  details <- sapply(dots, '[')
  details <- list(details = details)
  #
  attributes(n_t) <- c(attrib_n_t, state_names, state_rules, details)
  n_t
}

marginal_density <- function(n_t, margin = "x") {
  #
  attrib_n_t <- attributes(n_t)
  #
  n_t <- unlist(n_t)
  dim(n_t) <- unlist(attrib_n_t$details["num",])
  m <- which(attrib_n_t$state_names %in% margin)
  n_t <- apply(n_t, m, sum)
  # GENERALISE THIS
  dlta <- unlist(attrib_n_t$details["wgt",]) # <- only works if one weight
  if (!("x" %in% margin)) n_t <- n_t * dlta["x"]
  n_t
}

kernel_info <- function(state, what = "val", ...) {
  
  dots <- list(...)
  
  n <- length(dots)
  values <- list()
  suffix <- ""
  
  for (i in seq_along(dots)) {
    attribs <- attributes(dots[[n-i+1]])
    is_type <- attribs$state_rules["state",] == state
    keep <- attribs$details[what, is_type]
    names(keep) <- paste0(attribs$state_names[is_type], suffix)
    values[[n-i+1]] <- keep
    suffix <- paste0(suffix ,"_")
  }
  
  do.call(c, values)
}

assign_elements <- function(vec) {
  list2env(as.list(vec), envir = parent.frame())
}

ipm_fun <- function(mod_par, f_fix, f_env, f_dmg, ...) {
  #
  dots <- list(...)
  #
  cat_states <- do.call(kernel_info, c("categorical", what = "val", dots))
  num_states <- do.call(kernel_info, c("numeric",     what = "val", dots))
  #
  dims <- do.call(kernel_info, c("numeric", what = "num", dots))
  #
  assign_elements(do.call(expand.grid, num_states))
  # 
  assign_elements(mod_par)
  #
  expect_f <- list(expect_f = build_expect_f(num_states$x))
  assign_elements(expect_f)
  #
  a_set <- cat_states$a
  a_num <- length(a_set)
  
  f_fix <- f_fix[[2]]
  f_env <- f_env[[2]]
  f_dmg <- f_dmg[[2]]
  
  f_fix_eval <- mk_list_array(cat_states)
  for (a in head(a_set, -1)) {
    f_fix_eval[[a + 2, a + 1]] <- eval(f_fix)
  }
  f_fix_eval[[a_num, a_num]] <- eval(f_fix)
  
  #
  fun <- mk_list_array(cat_states)
  # 
  subset <- function(a_, a) fun[[a_, a]]
  # 
  function(i, n_t) {
    # 
    f_env_eval <- eval(f_env)
    # 
    for (a in head(a_set, -1)) {
      . <- f_fix_eval[[a + 2, a + 1]] + f_env_eval
      K <- eval(f_dmg)
      dim(K) <- dims
      fun[[a + 2, a + 1]] <<- K
    }
    . <- f_fix_eval[[a_num, a_num]] + f_env_eval
    K <- eval(f_dmg)
    dim(K) <- dims
    fun[[a_num, a_num]] <<- K
    #
    return(subset)
  }
}

age_transitions <- function(a_set) { 
  i <- as.character(c(tail(a_set, -1), tail(a_set, 1)))
  j <- as.character(a_set)
  list(i = i, j = j)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## simulation code
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

x_rule     <- numeric_midpoint(c(5, 40, 40))
a_rule     <- categorical_integer(c(0, 15))
x_lmb_rule <- numeric_midpoint(c(0.1, 4, 20))
s_rule     <- categorical_nominal(c("f", "m"))
t_rule     <- numeric_integer(c(1, 2))

n1 <- new_n_t(x = x_rule, a = a_rule)
n2 <- new_n_t(x = x_rule, a = a_rule)

su_f <- ipm_fun(par_fix$su,
                f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2,
                f_env = ~ (gamma + yr_ef[i,2]) * expect_f(n_t, theta) / x^theta,
                f_dmg = ~ inv_logit(.),
                n2, n1)

su_m <- ipm_fun(par_fix$su,
                f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2 +
                          b_0_m + b_z1_m * log(x) + b_a1_m * a + b_a2_m * a^2,
                f_env = ~ (gamma + yr_ef[i,2]) * expect_f(n_t, theta) / x^theta,
                f_dmg = ~ inv_logit(.),
                n2, n1)

gr_f <- ipm_fun(par_fix$gr,
                f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2,
                f_env = ~ (gamma + yr_ef[i,2]) * expect_f(n_t, theta) / x^theta,
                f_dmg = ~ dnorm(log(x_), mean = ., sd = sigma) / x_,
                n2, n1)

gr_m <- ipm_fun(par_fix$gr,
                f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2 +
                          b_0_m + b_z1_m * log(x) + b_a1_m * a + b_a2_m * a^2,
                f_env = ~ (gamma + yr_ef[i,2]) * expect_f(n_t, theta) / x^theta,
                f_dmg = ~ dnorm(log(x_), mean = ., sd = sigma) / x_,
                n2, n1)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Testing -- simulation code
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

set_lambs <- function(n_t, size, mean, sd) { 
  val <- attr(n_t, "details")[["val", "x"]]
  val
  n_t[["f", "0"]] <- size * dnorm(val, mean = mean, sd = sd)
  n_t[["m", "0"]] <- size * dnorm(val, mean = mean, sd = sd)
  n_t
}

# initial 
n_t <- new_n_t(x = x_rule, s = s_rule, a = a_rule)
n_t <- set_lambs(n_t, 300, 12, 1.5)


a_set <- attributes(n_t)$details[["val", "a"]]

e <- 1
n_t_1 <- new_n_t(x = x_rule, s = s_rule, a = a_rule)

system.time(for (i in 1:1e3) {
  # 
  n_z <- marginal_density(n_t, margin = "x")
  #
  set_j <- age_transitions(a_set)$j
  set_i <- age_transitions(a_set)$i
  # 
  su_f_t <- su_f(e, n_z)
  su_m_t <- su_m(e, n_z)
  gr_f_t <- gr_f(e, n_z)
  gr_m_t <- gr_m(e, n_z)
  #
  for (k in 1:nrow()) {
    # 
    j <- set_j[k]
    i <- set_i[k]
    #
    n_t_1[["f", i]] <- 
      ( su_f_t(i, j) * gr_f_t(i, j) ) %*% n_t[["f", j]]
    n_t_1[["m", i]] <- 
      ( su_m_t(i, j) * gr_m_t(i, j) ) %*% n_t[["m", j]]
  }
  # 
  n_t_1[["f", a_old]] <- n_t_1[["f", a_old]] +
    ( su_f_t(a_old, a_old) * gr_f_t(a_old, a_old) ) %*% n_t[["f", a_old]]
  n_t_1[["m", a_old]] <- n_t_1[["m", a_old]] +
    ( su_f_t(a_old, a_old) * gr_f_t(a_old, a_old) ) %*% n_t[["f", a_old]]
  
})





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Testing -- to be updated...
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












