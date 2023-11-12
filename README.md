# DCal
 DCal package is designed for treatment effect estimation with high-dimensional nuisance parameters. It uses a new Double-Calibration strategy that corrects the estimation bias of the nuisance parameter estimates computated by regularized high-dimensional techniques. Details of the methodology can be found in Lin Liu and Yuhao Wang (2023) "Root-n consistent semiparametric learning with high-dimensional nuisance functions under minimal sparsity" <arXiv:2305.04174>. The package relies on the optimisation software 'MOSEK' <https://www.mosek.com/> which must be installed separately; see the documentation for 'Rmosek'. 

## Installation
```R
devtools::install_github("Cinbo-Wang/DCal")
```

## Examples
We use the following examples to illustrate the basic usage of the DCal package. Here, the OR model is sparse linear model, while the propensity model is dense nonlinear.

```R
p = 300; s_or = 10; n = 100;rho=0.9;rd_num = 1
d_beta = s_or
Sigma_X <- matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    Sigma_X[i,j] <- rho ** abs(i-j)
  }
}

X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigma_X)
# dense propensity model
Xf <- X[,1:4]
Xf[,1] <- exp(0.5*X[,1])
Xf[,2] <- 10 + X[,2]/(1+exp(X[,1]))
Xf[,3] <- (0.05*X[,1]*X[,3]+0.6)**2
Xf[,4] <- (X[,2]+X[,4]+10)**2

gamma_true <- rep(0,p)
for(j in 1:p){
  gamma_true[j] <-  1/j
}
gamma_true <- gamma_true / norm(gamma_true,type='2')

lp <- scale(Xf[,1:4]) %*% c(1,-1/2,1/4,-1/8) + X %*% gamma_true
pi_W <- 1/(1+exp(-lp))
pi_W <- pmin(pmax(pi_W, 0.05), 0.95)
W <- rbinom(n = n,size=1,p=pi_W)

# sparse linear OR
beta_true <- rep(0,p)
act_loc <- 1:d_beta # Confounder
beta_true[act_loc] <- runif(d_beta,1,2)*sample(c(1,-1),size=d_beta,replace = T)
beta_true <-  beta_true / norm(beta_true,type='2')
potential_outcome_treat <- X %*% beta_true + 1
potential_outcome_control <- X %*% beta_true
Y <- potential_outcome_treat*W + potential_outcome_control*(1-W) + rnorm(n,0,1)
tau_treat <- mean(potential_outcome_treat)
mean_treat_cf_dcal_ls <- DCal.mean_treat(X,Y,W,B=3,r1 = NULL,pi = NULL,
                                            is.scale = T,Y.family = 'gaussian',alpha = 0.9, is.parallel=F,core_num=3)

mean_treat_cf_dcal_ls

```