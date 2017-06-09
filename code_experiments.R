
library(lazyeval)

x_rule <- numeric_midpoint(c(5, 40, 30))
a_rule <- categorical_integer(c(0, 12))

n1 <- new_n_t(x = x_rule, a = a_rule)
n2 <- new_n_t(x = x_rule, a = a_rule)

cat_states <- do.call(kernel_domain, c("categorical", list(n1, n2)))
num_states <- do.call(kernel_domain, c("numeric",     list(n1, n2)))

#
a_set <- cat_states$a

fun_1 <-  ~ b_0_f + b_z1_f * log(x) + b_a1_f * a + b_a2_f * a^2
fun_2 <-  ~ (gamma)
fun_3 <-  ~ inv_logit(f_all)

fun_1 <- f_rhs(fun_1)
fun_2 <- f_rhs(fun_2)
fun_3 <- f_rhs(fun_3)

env_params <- new.env()
env_domain <- new.env(parent = env_params)
env_dummy  <- new.env(parent = env_domain)

list2env(par_fix$su, envir = env_params)
list2env(do.call(expand.grid, num_states), envir = env_domain)

ls.str(env_params)
ls.str(env_domain)

list_array <- mk_list_array(cat_states)

state_values <- cbind(a_ = c(tail(a_set, -1), max(a_set)), 
                      a  = a_set)

state_index <- as.character(state_values)
attributes(state_index) <- attributes(state_values)

create_demog_func <- function(envir, f1, f2, f3, 
                              list_array, state_index, state_values) {
  
  enclos <- baseenv()
  
  f1_eval <- list_array

  i_seq <- seq_len(nrow(state_index))

  # for (i in i_seq) {
  #   .Internal(assign("a", state_values[i, "a"], envir, FALSE))
  #   f1_eval[ state_index[i, , drop = FALSE] ] <- 
  #     list(.Internal(eval(f1, envir, enclos)))
  # }
  
  f1_eval[ state_index ] <- 
    lapply(i_seq, function(i) {
      .Internal(assign("a", state_values[i, "a"], envir, FALSE))
      .Internal(eval(f1, envir, enclos))
    })
  
  results <- f_all <- list_array
  
  getter <- function(a) results[[a+2, a+1]]
  
  function(i) {
    # 
    f2_eval <- .Internal(eval(fun_2, envir, enclos))
  
    #
    f_all[ state_index ] <<- 
      mapply(`+`, f1_eval[ state_index ], f2_eval, SIMPLIFY = FALSE)

    results[ state_index ] <<- 
      lapply(i_seq, function(i) {
        .Internal(assign("a", state_values[i, "a"], envir, FALSE))
        .Internal(assign("f_all", unlist(f_all[ state_index[i, , drop = FALSE] ]) , envir, FALSE))
        .Internal(eval(f3, envir, enclos))
      })
    
    getter
  }
}

system.time({
  f1 <- create_demog_func(env_domain, fun_1, fun_2, fun_3, 
                          list_array, state_index, state_values)
  for (i in 1:10000) {f2 <- f1(1); for (j in 1:10) f2(j)}
})


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## general messing
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

all.vars()


pars <- list(a = 1, b = 2)
node <- list(x = 1:1000)

env_nodes <- new.env()
env_evals <- new.env(parent = env_nodes)

list2env(node, envir = env_nodes)
list2env(pars, envir = env_evals)

## using base for the evaluation

form <- ~ {y <- a + b * x}
expr <- f_rhs(form)
all.vars(form)

system.time(
  for (i in 1:1e4) { # number of iterations
    for (j in 1:15) { # number of categorical states
      for (k in 1:1) { # number of demographic functions
        eval(expr, env_evals)
      } 
    }
})


env_evals$b








  

