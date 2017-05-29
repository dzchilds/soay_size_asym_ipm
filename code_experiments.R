### fastest way to store discrete states

n1 <- 25 # n mesh
n2 <- 12 # n states
n_sim <- 1e4
x_sml <- 1:n1
x_big <- 1:(n1*n2)

M_sml <- matrix(rnorm(n1^2), nrow = n1, ncol = n1)
M_big <- matrix(rnorm(n1 * n1 * n2), nrow = n1, ncol = n1 * n2)

system.time(
  for (i in 1:n_sim) {
    for (j in 1:n2) {
      (M_sml * M_sml * M_sml * M_sml) %*% x_big[x_sml+2*n1]
    } 
    (M_big * M_big * M_big * M_big) %*% x_big
  })

system.time(
  for (i in 1:n_sim) {
    for (j in 1:n2) {
      (M_sml * M_sml * M_sml * M_sml) %*% x_sml
      (M_sml * M_sml * M_sml * M_sml) %*% x_sml
    } 
  })



x <- 1:(n1*n2)
n_calc <- n_sim
system.time(for (i in 1:n_calc) (M * M * M) %*% x)

###

x <- list()
x$f <- vector(mode = "list", length = 20)
names(x$f) <- letters[1:20]
x$m <- vector(mode = "list", length = 20)
names(x$m) <- letters[1:20]
system.time(for (i in 1:1e7) x[["m"]][["d"]] <- addto)


x <- mk_list_array(list(c("f","m"), letters[1:20]))
addto <- 1:1000
system.time(for (i in 1:1e7) x[["m","a"]] <- addto)


x <- mk_list_array(list(c("f","m"), letters[1:20]))
dum_list <- vector(mode = "list", length = 21)
names(dum_list) <- letters[1:21]
x["f",] <- dum_list
x
