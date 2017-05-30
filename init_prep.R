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

numeric_integer <- function(x) {
  rule <- list(seq.int(x[1], x[2]), 1, x[2]-x[1]+1)
  names(rule) <- c("val", "wgt", "num")
  attr(rule, "state") <- "numeric"
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

new_n_t <- function(...) {
  #
  dots <- list(...) 
  # work out where to 
  type_info <- sapply(dots, function(x) attributes(x)[c("state","rule")])
  case1 <- type_info["state",] == "categorical"
  case2 <- type_info["state",] == "numeric" & type_info["rule",] == "integer"
  
  list_part <- dots[case1 | case2]
  n_t <- mk_list_array(lapply(list_part, '[[', "val"))
  attrib_n_t <- attributes(n_t)
  
  vect_part <- dots[!(case1 | case2)]
  size <- prod(sapply(vect_part, '[[', "num"))
  zeros <- rep(0, size)
  n_t <- lapply(n_t, function(dummy) zeros)
  
  attributes(n_t) <- c(attrib_n_t, dots)
  n_t
}

# careful -- order of states must match, here and in new_n_t
marginal_density <- function(n_t, margin = "x") {
  # grab the bits we need
  dims <- c(attr(n_t, "x")$all$num, attr(n_t, "s")$num, attr(n_t, "a")$num)
  dx <- attr(n_t, "x")$all$wgt
  #
  n_t <- unlist(n_t)
  dim(n_t) <- dims
  m <- which(c("x", "s", "a") %in% margin)
  n_t <- apply(n_t, m, sum)
  # 
  if (!("x" %in% margin)) n_t <- n_t * dx
  n_t
}


extract_rules <- function(rules, match_str) {

  splits <- strsplit(match_str, "_", fixed = TRUE)
  depth <- sapply(splits, length)
  pref <- sapply(splits, `[`, 1)
  suff <- sapply(splits, `[`, 2)
  out <- list()
  for (i in seq_along(match_str)) {
    if (depth[i]==1) {
      out[[i]] <- rules[[pref[i]]] $ val
    } 
    if (depth[i]==2) {
      out[[i]] <- rules[[ pref[i] ]][[ suff[i] ]] $ val
    }
  }
  names(out) <- pref
  out
}

kern_nodes <- function(n_t, codomain, domain) {
  rules <- attributes(n_t)
  cd <- extract_rules(rules, codomain)
  dm <- extract_rules(rules, domain)
  names(cd) <- paste0(names(cd), "1")
  c(cd, dm)
}


set_lambs <- function(n_t, size, mean, sd) { 
  val <- attr(n_t, "x")$val
  n_t[["f", '0']] <- size * dnorm(val, mean = mean, sd = sd) / 2
  n_t[["m", '0']] <- size * dnorm(val, mean = mean, sd = sd) / 2
  n_t
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## simulation code
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

x <- numeric_midpoint(c(5,  40, 40))
a <- numeric_integer(c(0, 15))
s <- categorical_nominal(c("f", "m"))

n_t_0 <- new_n_t(x = x, s = s, a = a)
n_t_0 <- set_lambs(n_t_0, 200, 12, 1.5)
n_t_0[["f","0"]]

for (i in 1:10000) marginal_density(n_t_0, c("s", "x"))

node_sg <- kern_nodes(n_t_0, c("x_all"), c("x_all", "a"))

su <- build_su(node_sg, par_fix$su)
gr <- build_gr(node_sg, par_fix$gr)


su_t <- su(1, n_t)
gr_t <- gr(1, n_t)

a <- 0
P <- su_t("f", a) * gr_t("f", a)
dim(P) <- 
P %*% get_age(n_t$f, 0)



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












