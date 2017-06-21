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
  function(n, theta) sum(n * x^theta) * dx
}

# 
# lamb size function (female only)


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

mapply_expr <- function(expr, variable_args, fixed_args = list()) {
  #
  mapply_args <- list(
    FUN = function(...) eval(expr, envir = list(...)),
    SIMPLIFY = FALSE, MoreArgs = fixed_args
  ) %>% c(variable_args)
  # 
  do.call(mapply, mapply_args)
}

# never used this...
make_alist <- function(arg_names) {
  if ( any(arg_names == "") ) stop("missing argument name")
  args <- replicate(length(arg_names), substitute())
  names(args) <- arg_names
  args
}

ipm_fun <- function(mod_par, f_fix, f_env, f_dmg, cat_trans, ...) {
  #
  dots <- list(n2, n1)
  # extract the information we require from the state atrtributes
  cat_states <- do.call(kernel_info, c("categorical", what = "val", dots))
  num_states <- do.call(kernel_info, c("numeric",     what = "val", dots))
  array_dims <- do.call(kernel_info, c("numeric",     what = "num", dots))
  # list containing objects that is invariant across function evaluations
  fixed_args <- c(do.call(expand.grid, num_states), mod_par)
  # list containing objects that vary across 
  variable_args <- cat_trans
  # FIXME -- hack to make the numeric age states a numeric vector
  variable_args$a_ <- as.numeric(cat_trans$a_)
  variable_args$a  <- as.numeric(cat_trans$a )
  # FIXME -- function to calculate the expectation wrt to density function
  expect_f <- build_expect_f(num_states$x)
  # list-array to store the function evaluations
  eval_func <- mk_list_array(cat_states)
  # subset function for the above list-array
  subset <- function(...) {
    do.call(`[[`, c(x = quote(eval_func), list(...)))
  }
  # matrix to index into list of evaluated functions 
  index <- as.matrix(cat_trans)
  # evaluate the time invariant part of the additive model
  f_fix_eval <- mapply_expr(f_fix[[2]], variable_args, fixed_args)
  # 
  function(...) {
    # 
    f_args <- c(list(...), fixed_args, expect_f = expect_f)
    .additive <- 
      mapply_expr(f_env[[2]], list(0), f_args) %>% 
      mapply(`+`, f_fix_eval, ., SIMPLIFY = FALSE) 
    #
    v_args <- c(.additive = list(.additive), variable_args)
    eval_func[index] <<- 
      mapply_expr(f_dmg[[2]], v_args, fixed_args) %>%
      lapply(function(x) {
        dim(x) <- array_dims
        x
      })
    # FIXME -- the above `add dim` stpe is not very efficient
    return(subset)
  }
}


ipm_kern <- function(kern, funs, cat_trans, ...) {
  # evaluate set of demographic functions
  fixed_args <- lapply(funs, do.call, list(...))
  # evaluate kernel over categorical states
  mapply_expr(kern[[2]], cat_trans, fixed_args)  
}

# FIXME
iter_kern <- function(kern, funs, cat_trans) {
  # 
  force(funs)
  #
  trans_from <- as.character(cat_trans$a)
  to_states <- as.character(unique(cat_trans$a_))
  trans_to <- factor(cat_trans$a_, levels = to_states)
  # 
  iter_expr <- as.call(list(quote(`%*%`), kern[[2]], quote(n_t)))
  #
  mapply_args <- list(
    FUN = function(...) eval(iter_expr, envir = list(...)),
    SIMPLIFY = FALSE, MoreArgs = NA
  ) %>% c(cat_trans, n_t = NA)
  # 
  function(.nt, ...) {
    # add the evaluated demographic functions to the list of mapply arguments
    mapply_args$MoreArgs <<- lapply(funs, do.call, list(...))
    # pad out the state list so that it matches number of transitions + add
    mapply_args$n_t <<- .nt[trans_from]
    # 
    do.call(mapply, mapply_args) %>% 
      split(trans_to) %>%
      lapply(function(x) Reduce(`+`, x))
  }
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## simulation code
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

calculate_P_transitions <- function(a_set) { 
  a_ <- c(tail(a_set, -1), tail(a_set, 1))
  a_ <- as.character(a_)
  a  <- a_set
  a  <- as.character(a)
  data.frame(a_ = a_, a = a)
}

calculate_F_transitions <- function(a_set) { 
  a_  <- as.character(a_set)
  t_ <- c("s", "t")
  s_ <- c("f", "m")
  out <- expand.grid(s_ = s_, t_ = t_, a_ = a_, stringsAsFactors = FALSE)
  cbind(out, a = out$a_)
}


min_a <- 0
max_a <- 15

a_rule     <- categorical_integer(c(min_a, max_a))
x_rule     <- numeric_midpoint(c(5, 40, 40))
s_rule     <- categorical_nominal(c("f", "m"))
t_rule     <- categorical_nominal(c("s", "t"))
x_lmb_rule <- numeric_midpoint(c(0.1, 4, 20))

##
## 1. define up the survival-growth iteration functions
##

P_transitions  <- calculate_P_transitions(min_a:max_a)

n2 <- n1 <- new_n_t(x = x_rule, a = a_rule)

su_f <- ipm_fun(par_fix$su,
                f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2,
                f_env = ~ yr_ef[i,1] + (gamma + yr_ef[i,2]) * expect_f(n, theta) / x^theta,
                f_dmg = ~ inv_logit(.additive),
                P_transitions, n2, n1)

su_m <- ipm_fun(par_fix$su,
                f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2 +
                          b_0_m + b_z1_m * log(x) + b_a1_m * a + b_a2_m * a^2,
                f_env = ~ yr_ef[i,1] + (gamma + yr_ef[i,2]) * expect_f(n, theta) / x^theta,
                f_dmg = ~ inv_logit(.additive),
                P_transitions, n2, n1)

gr_f <- ipm_fun(par_fix$gr,
                f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2,
                f_env = ~ yr_ef[i,1] + (gamma + yr_ef[i,2]) * expect_f(n, theta) / x^theta,
                f_dmg = ~ dnorm(log(x_), mean = .additive, sd = sigma) / x_,
                P_transitions, n2, n1)

gr_m <- ipm_fun(par_fix$gr,
                f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2 +
                          b_0_m + b_z1_m * log(x) + b_a1_m * a + b_a2_m * a^2,
                f_env = ~ yr_ef[i,1] + (gamma + yr_ef[i,2]) * expect_f(n, theta) / x^theta,
                f_dmg = ~ dnorm(log(x_), mean = .additive, sd = sigma) / x_,
                P_transitions, n2, n1)




P_iter_m <- 
  iter_kern(kern = ~ su(a_, a) * gr(a_, a), 
            funs = list(su = su_m, gr = gr_m), 
            P_transitions)

##
## 2. define up the fecundity iteration functions
##

n1 <- new_n_t(x = x_rule, a = a_rule)
n2 <- new_n_t(s = s_rule, t = t_rule, x = x_lmb_rule, a = a_rule)

F_transitions <- calculate_F_transitions(min_a:max_a)

su <- ipm_fun(par_fix$su,
              f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2,
              f_env = ~ yr_ef[i,1] + (gamma + yr_ef[i,2]) * expect_f(n, theta) / x^theta,
              f_dmg = ~ inv_logit(.additive),
              F_transitions, n2, n1)

rp <- ipm_fun(par_fix$rp,
              f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2,
              f_env = ~ yr_ef[i,1] + (gamma + yr_ef[i,2]) * expect_f(n, theta) / x^theta,
              f_dmg = ~ inv_logit(.additive),
              F_transitions, n2, n1)

tw <- ipm_fun(par_fix$tw,
              f_fix = ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_az_f * a * log(x),
              f_env = ~ yr_ef[i,1] + (gamma + yr_ef[i,2]) * expect_f(n, theta) / x^theta,
              f_dmg = ~ inv_logit(.additive),
              F_transitions, n2, n1)

bw <- ipm_fun(par_fix$bw_lmb,
              f_fix = ~ b_0_f + (s_ == "m") * b_0_m + (t_ == "t") * b_0_tw + 
                b_a1_f * a + b_z1_f * log(x) + b_az_f * a * log(x) + b_a2_f * a^2,
              f_env = ~ yr_ef[i,1] + (gamma + yr_ef[i,2]) * expect_f(n, theta) / x^theta,
              f_dmg = ~ dnorm(log(x_), .additive, sigma) / x_,
              F_transitions, n2, n1)

F1_iter <- 
  iter_kern(kern = ~ bw(s_, t_, a_, a) * tw(s_, t_, a_, a) * 
                     rp(s_, t_, a_, a) * su(s_, t_, a_, a), 
            funs = list(su = su, rp = rp, tw = tw, bw = bw), 
            F_transitions)

##
## 3. lamb spring-summer iteration functions
##

su_lmb <- ipm_fun(par_fix$su_lmb,
                  f_fix = ~ b_0_f + is_m * b_0_m + is_tw * b_0_tw + 
                    b_z1_lm * log(x) + b_a1_f * a + b_a2_f * a^2,
                  f_env = ~ yr_ef[i,1],
                  f_dmg = ~ inv_logit(.additive),
                  F_transitions, n2, n1)


gr_lmb <- ipm_fun(par_fix$gr_lmb,
                  f_fix = ~ b_0_f + is_m * b_0_m + is_tw * b_0_tw + 
                    b_z1_lm * log(x) + b_a1_f * a + b_a2_f * a^2,
                  f_env = ~ yr_ef[i,1],
                  f_dmg = ~ dnorm(log(x_), .additive, sigma) / x_,
                  F_transitions, n2, n1)

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

index_matrix <- function(...) {
  x <- list(...)
  x <- do.call(cbind, x)
  which_keep <- grepl("_{1}", names(x))
  x <- unique(x[which_keep])
  as.matrix(x)
}

# initial state 
n_t_0 <- new_n_t(s = s_rule, x = x_rule, a = a_rule)
# intermediate state
n_t_I <- new_n_t(s = s_rule, t = t_rule, x = x_lmb_rule, a = a_rule)
# final state
n_t_1 <- new_n_t(s_rule, x = x_rule, a = a_rule)
# initialise
n_t_0 <- set_lambs(n_t, 300, 12, 1.5)
# indexing vectors
P_index_f <- index_matrix(s_ = "f", P_transitions)
P_index_m <- index_matrix(s_ = "m", P_transitions)
F_index   <- index_matrix(F_transitions)

n_t <- n_t_0
#
n_z <- marginal_density(n_t, margin = "x")
#
P_kern_f <- 
  ipm_kern(kern = ~ su(a_, a) * gr(a_, a), 
           funs = list(su = su_f, gr = gr_f), 
           P_transitions, i = 1, n = n_z)


sapply(P_kern_f, dim)

