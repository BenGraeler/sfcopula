# adopted from http://freakonometrics.hypotheses.org/2435, 

lower_emp_biv_joint_dep_fun <- function(u) {
  stopifnot(ncol(u) == 2)
  empFun <- function(x) sum((u[,1]<=x)&(u[,2]<=x))/sum(u[,1]<=x)
  function(x) sapply(x,empFun)
}

upper_emp_biv_joint_dep_fun <- function(u) {
  stopifnot(ncol(u) == 2)
  empFun <- function(x) sum((u[,1]>=x)&(u[,2]>=x))/sum(u[,1]>=x)
  function(x) sapply(x,empFun)
}

emp_biv_joint_dep_fun <- function(u) {
  stopifnot(ncol(u) == 2)
  
  function(z) {
    res <- z
    res[z>0.5] <- upper_emp_biv_joint_dep_fun(u)(z[z>0.5])
    res[z<=0.5] <- lower_emp_biv_joint_dep_fun(u)(z[z<=0.5])
    return(res)
  }
}

##

lower_biv_joint_dep_fun <- function(copula) {
  stopifnot(copula@dimension == 2)
  function(z) pCopula(cbind(z,z),copula)/z
}

upper_biv_joint_dep_fun <- function(copula) {
  stopifnot(copula@dimension == 2)
  function(z) (1-2*z+pCopula(cbind(z,z),copula))/(1-z)
}

biv_joint_dep_fun <- function(copula) {
  function(z) {
    res <- z
    res[z>0.5] <- upper_biv_joint_dep_fun(copula)(z[z>0.5])
    res[z<=0.5] <- lower_biv_joint_dep_fun(copula)(z[z<=0.5])
    return(res)
  }
}