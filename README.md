# DCal
 DCal package is designed for treatment effect estimation with high-dimensional nuisance parameters. It uses a new Double-Calibration strategy that corrects the estimation bias of the nuisance parameter estimates computated by regularized high-dimensional techniques. Details of the methodology can be found in Lin Liu, Xinbo Wang and Yuhao Wang (2023) "Root-n consistent semiparametric learning with high-dimensional nuisance functions under minimal sparsity"[arXiv link](https://arxiv.org/abs/2305.04174). The package relies on the optimisation software ['MOSEK'](https://www.mosek.com/) which must be installed separately; see the documentation for 'Rmosek'. 

## Installation
```R
devtools::install_github("Cinbo-Wang/DCal")
```

## Examples
We use the following examples to illustrate the basic usage of the DCal package. Here, the OR model is sparse linear model, while the propensity model is dense nonlinear.

```R
p = 400; s_or = 10; n = 200;rho=0.9;rd_num = 1
Sigma_X <- matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    Sigma_X[i,j] <- rho ** abs(i-j)
  }
}
X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigma_X)

set.seed(rd_num)
# Sparse PS

gamma_true <- rep(0,p)
act_loc <- 1:s_or # sample(1:p,size=d_gamma,replace = F)
gamma_true[act_loc] <- runif(length(act_loc),1,2)#*sample(c(1,-1),size=d_gamma,replace=T)
gamma_true <- gamma_true / norm(gamma_true,type='2')

U <- X %*% gamma_true
pi_W <- 1 / (1+exp(-U))
W <- rbinom(n,size=1,p=pi_W)

# nonlinear dense OR
require(splines)
beta_true <- rep(0,p)
act_loc <- 1:p 
for(j in act_loc){
  beta_true[j] <- 1/j 
}
beta_true<-  beta_true / norm(beta_true,type='2')

X2 <- 2 / (1+exp(-X[,2]))
X3 <- exp(X[,3]/2)
X4 <- X[,4]/(1+exp(X[,3]))
X5 <- X[,4] * X[,5] / 10 

non_lin <- scale(cbind(X2,X3,X4,X5)) %*% c(1,-1/2,1/3,-1/4) + bs(X[, 1], df = 100) %*% c(1/(1:100))
term_cm <- (abs(non_lin)+0.05)**(-1) + X %*% beta_true 
potential_outcome_treat <-  term_cm - 1
potential_outcome_control <- term_cm
Y <- potential_outcome_treat * W + potential_outcome_control * (1-W) + rnorm(n,0,1)
tau_treat_g <- mean(potential_outcome_treat)

# Slow version
mean_treat_dcal_ls <- DCal.mean_treat(X,Y,W,Y.family = 'gaussian',B=2,is.scale = FALSE,alpha = 0.9, is.parallel=TRUE,core_num = 2)
ate_vec <- mean_treat_dcal_ls$ATE_vec
sd_vec <- sqrt(mean_treat_dcal_ls$var_ATE_vec)
CI_lb_vec <- ate_vec - qnorm(0.975) * sd_vec
CI_ub_vec <- ate_vec + qnorm(0.975) * sd_vec
ate_vec
CI_lb_vec
CI_ub_vec

# Faster version
mean_treat_dcal_ls <- DCal_star.mean_treat(X,Y,W,Y.family = 'gaussian',B=6,is.scale = FALSE,alpha = 0.9, is.parallel=FALSE)
ate_vec <- mean_treat_dcal_ls$ATE_vec
sd_vec <- sqrt(mean_treat_dcal_ls$var_ATE_vec)
CI_lb_vec <- ate_vec - qnorm(0.975) * sd_vec
CI_ub_vec <- ate_vec + qnorm(0.975) * sd_vec
ate_vec
CI_lb_vec
CI_ub_vec


```
