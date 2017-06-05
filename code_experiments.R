
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
fun_2 <-  ~ (gamma + yr_ef[i,2]) * expect_f(n_t, theta) / x^theta
fun_3 <-  ~ inv_logit(f_all)

fun_1 <- f_rhs(fun_1)
fun_2 <- f_rhs(fun_2)
fun_3 <- f_rhs(fun_3)

env_params <- new.env()
env_domain <- new.env(parent = env_params)

list2env(par_fix$su, envir = env_params)
list2env(do.call(expand.grid, num_states), envir = env_domain)

ls.str(env_params)
ls.str(env_domain)

list_array <- mk_list_array(cat_states)

create_demog_func <- function(envir, f1, f2, f3, list_array, state_index) {
  
  enclos <- baseenv()
  
  a_num <- max(state_index)
  
  f1_eval <- list_array
  
  for (a in head(state_index, -1)) {
    .Internal(assign("a", a, envir, FALSE))
    f1_eval[[a + 2, a + 1]] <- .Internal(eval(f1, envir, enclos))
  }
  .Internal(assign("a", a_num, envir, FALSE))
  f1_eval[[a_num, a_num]] <- .Internal(eval(f1, envir, enclos))
  
  f1_eval
}

system.time(
for (i in 1:10000) 
  create_demog_func(env_domain, fun_1, fun_2, fun_3, list_array, a_set)
)


