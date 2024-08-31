
#' Estimation of E[Y(1)] from observational data using Double Calibration Estimator.
#'
#' @param X The n by p input data matrix.
#' @param Y The n dimensional observed response.
#' @param W The n dimensional binary vector indicating treatment assignment.
#' @param Y.family Response type. A character string representing one of the built-in families, including 'gaussian' and 'binomial'.
#' @param B The number of iterations for random splits for cross-fitting.
#' @param is.scale Whether or not to scale the input data matrix.
#' @param r1_init Optional n dimensional vector of an initial estimate of E[Y_i |
#'   X_i, T_i=1] for i = 1, ..., n. The default is NULL.
#' @param pi_init Optional n dimensional vector of an initial estimate of E[W_i |
#'   X_i] for i = 1, ..., n. The default is NULL.
#' @param alpha The elastic net mixing parameter, with \(0 \eqn{\leq} alpha \eqn{\leq} 1\).
#' @param is.parallel Whether to perform parallel computation. Default is False.
#' @param core_num Number of cores used for parallel computation.
#' @param ratio_violate_dcal The maximum allowable ratio of violations in inequations during the Double Calibrated step. The default value is 0.01.
#' @return
#' \item{ATE_sc}{Single calibation estimation of E[Y(1)] averaged over B random splits.}
#' \item{ATE_sc_var}{Estimated variance of the ATE_sc averaged over B random splits.}
#' \item{ATE_dc}{Double calibation estimation of E[Y(1)] averaged over B random splits.}
#' \item{ATE_dc_var}{Estimated variance of the ATE_dc averaged over B random splits.}
#' \item{ATE_mat}{A matrix containing ATE_sc, ATE_sc_var, ATE_dc and ATE_dc_var in each random split.}
#'
#' @author Xinbo Wang
#' @references Lin Liu, and Yuhao Wang. (2023) \emph{Root-n consistent semiparametric learning with
#' high-dimensional nuisance functions under minimal sparsity.} \url{https://doi.org/10.48550/arXiv.2305.04174}


#' @examples
#' \dontrun{
#' # Sparse OR, dense nonlinear PS--------
#' p = 400; s_or = 10; n = 200;rho=0.9;rd_num = 1
#' Sigma_X <- matrix(0,p,p)
#' for(i in 1:p){
#'   for(j in 1:p){
#'     Sigma_X[i,j] <- rho ** abs(i-j)
#'   }
#' }
#'
#' X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigma_X)
#' # dense propensity model
#' Xf <- X[,1:4]
#' Xf[,1] <- exp(0.5 * X[,1])
#' Xf[,2] <- 10 + X[,2] / (1 + exp(X[,1]))
#' Xf[,3] <- (0.05 * X[,1] * X[,3] + 0.6)**2
#' Xf[,4] <- (X[,2] + X[,4] + 10)**2
#'
#' gamma_true <- rep(0,p)
#' for(j in 1:p){
#'   gamma_true[j] <-  1/j
#' }
#' gamma_true <- gamma_true / norm(gamma_true,type='2')
#'
#' lp <- scale(Xf[,1:4]) %*% c(1,-1/2,1/4,-1/8) + X %*% gamma_true
#' pi_W <- 1/(1+exp(-lp))
#' pi_W <- pmin(pmax(pi_W, 0.05), 0.95)
#' W <- rbinom(n = n,size=1,p=pi_W)
#'
#' # sparse linear OR
#' beta_true <- rep(0,p)
#' act_loc <- 1:s_or # Confounder
#' beta_true[act_loc] <- runif(s_or,1,2)
#' beta_true <-  beta_true / norm(beta_true,type='2')
#' potential_outcome_treat <- X %*% beta_true + 1
#' potential_outcome_control <- X %*% beta_true - 1
#' Y <- potential_outcome_treat * W + potential_outcome_control * (1-W) + rnorm(n,0,1)
#' tau_treat <- mean(potential_outcome_treat)
#'
#' mean_treat_dcal_ls <- DCal.mean_treat(X,Y,W,B=6,r1_init = NULL,pi_init = NULL,
#'                      is.scale = FALSE,Y.family = 'gaussian',alpha = 0.9, is.parallel=FALSE)
#'
#' mean_treat_dcal_ls
#' CI_dcal <- c(
#'       mean_treat_cf_dcal_ls$ATE_dc - qnorm(0.975)*sqrt(mean_treat_cf_dcal_ls$ATE_dc_var),
#'       mean_treat_cf_dcal_ls$ATE_dc + qnorm(0.975)*sqrt(mean_treat_cf_dcal_ls$ATE_dc_var)
#' )
#' tau_treat;
#' CI_dcal
#' }

#' @export DCal.mean_treat
#' @import glmnet
#' @import foreach
#' @import doParallel
#' @import stats
DCal.mean_treat <-
  function(X,
           Y,
           W,
           Y.family = c('gaussian', 'binomial')[1],
           B = 6,
           is.scale = TRUE,
           r1_init = NULL,
           pi_init = NULL,
           alpha = 0.9,
           is.parallel = FALSE,
           core_num = NULL,
           ratio_violate_dcal=0.01) {
    # require(glmnet)
    if (is.scale) {
      scl <- apply(X, 2, sd, na.rm = TRUE)
      is.binary <-
        apply(X, 2, function(xx)
          sum(xx == 0) + sum(xx == 1) == length(xx))
      scl[is.binary] <- 1
      X <- scale(X, center = FALSE, scale = scl)
    }


    if (is.null(r1_init)) {
      fit.out.treated <-
        cv.glmnet(X[W == 1,],
                  Y[W == 1],
                  family = Y.family,
                  alpha = alpha,
                  nfolds = 5)
      r1_out <-
        predict(fit.out.treated, newx = X, type = 'response')
    } else{
      r1_out <- r1_init
    }

    if (is.null(pi_init)) {
      fit.prop <-
        cv.glmnet(X,
                  W,
                  family = "binomial",
                  alpha = alpha,
                  nfolds = 5)
      pi_hat <- predict(fit.prop, newx = X, type = 'response')
      pi_hat <- pmax(pmin(pi_hat, 0.99), 0.01)
    } else{
      pi_hat <-  pmax(pmin(pi_init, 0.99), 0.01)
    }


    loc0 <- which(W == 0)
    loc1 <- which(W == 1)
    if (is.parallel) {
      type <-
        ifelse(.Platform$OS.type == 'windows', 'PSOCK', 'FORK')
      core_num <-
        ifelse(!is.null(core_num),
               core_num,
               ifelse(.Platform$OS.type == 'windows', 4, min(25, B)))
      cl <- parallel::makeCluster(core_num, type = type)
      doParallel::registerDoParallel(cl)
      ATE_mat <-
        foreach::foreach(b = 1:B, .combine = 'rbind',.export = c('quad.prog.Lagr','double_cali.Lagr')) %dopar% {
          set.seed(b)

          # estimate mu using cross-fitting
          mu_hat <- rep(0, nrow(X))
          loc0_2_group <- sample(1:2, size = length(loc0), replace = TRUE)
          loc1_2_group <- sample(1:2, size = length(loc1), replace = TRUE)

          for (foldid in 1:2) {
            loc_main <- c(loc0[loc0_2_group == foldid], loc1[loc1_2_group == foldid])
            loc_aux <- c(loc0[loc0_2_group != foldid], loc1[loc1_2_group != foldid])

            loc_main <- sample(loc_main, size = length(loc_main))
            loc_aux <- sample(loc_aux, size = length(loc_aux))

            X_main <- cbind(1, scale(X[loc_main,], scale = FALSE))
            X_aux <- cbind(1, scale(X[loc_aux,], scale = FALSE))

            W_main <- W[loc_main]
            Y_main <- Y[loc_main]
            W_aux <- W[loc_aux]
            Y_aux <- Y[loc_aux]

            n_main <- nrow(X_main)
            p <- ncol(X_main)
            n_aux <- nrow(X_aux)

            r1_aux <- r1_out[loc_aux]

            pi_hat_aux <- pi_hat[loc_aux]
            pi_hat_main <- pi_hat[loc_main]

            Pi_aux <- diag(1 - pi_hat_aux)
            R_aux <- W_aux * (Y_aux - r1_aux) / pi_hat_aux
            Pi_main <- diag(1 - pi_hat_main)

            target <-  t(X_aux) %*% Pi_aux %*% R_aux / n_aux


            ratio_violate <- sum(abs(target)>0.5*sqrt(log(p) / n_main)) / length(target)
            if (ratio_violate< 0.01) {
              # If constraint is satisfied.
              mu_hat_main <- rep(0, n_main)
            } else{
              M <- t(X_main) %*% Pi_main / n_main
              weight_term <- ifelse(ratio_violate>0.1,2,3/2)
              mu_hat_main_tmp_ls <-
                quad.prog.Lagr(
                  M = M,
                  target = target,
                  Mr = max(abs(Y)),
                  kappa = 0.5,
                  verbose = FALSE,
                  weight_term = weight_term
                )
              if (mu_hat_main_tmp_ls$is_converge) {
                mu_hat_main <- mu_hat_main_tmp_ls$mu_hat
              } else{
                mu_hat_main <- rep(0, n_main)
              }
            }

            mu_hat[loc_main] <- mu_hat_main
          }

          mean_treat_init <-
            sum(W * (Y - mu_hat - r1_out) / pi_hat) / sum(W / pi_hat) + sum(mu_hat + r1_out) / nrow(X)
          mean_treat_init_vec <- W * (Y - mu_hat - r1_out) / pi_hat + mu_hat + r1_out
          mean_treat_init_var <- sum((mean_treat_init_vec - mean_treat_init) ** 2) / (nrow(X) ** 2)

          # Double Calibration for pi_tilde (weight_tilde)
          weight_init <- 1 / pi_hat
          if (Y.family == 'binomial') {
            psi_deri <- as.numeric(r1_out  * (1 - r1_out))
          } else{
            psi_deri <- rep(1, nrow(X))
          }
          X_tilde <- cbind(diag(psi_deri) %*% X, mu_hat)
          n <- nrow(X_tilde)
          p <- ncol(X_tilde)
          lambda_vec <- sqrt(log(p) / n) * sqrt(colSums(X_tilde ** 2)) / sqrt(n)

          constr_left <- abs(colSums(diag(as.numeric(W / pi_hat - 1)) %*% X_tilde)) / n
          ratio_violate <- sum(constr_left > lambda_vec)/length(constr_left)

          if (ratio_violate < ratio_violate_dcal) {
            weight_tilde <- weight_init
          } else{
            weight_ls <-
              double_cali.Lagr(X,
                               W,
                               mu_hat,
                               psi_deri,
                               kappa = 0.5,
                               verbose = FALSE)
            weight_tilde <- weight_init
            if (weight_ls$is_converge) {
              weight_tilde[W == 1] <- weight_ls$weight_tilde
            }
          }

          mean_treat_ped <- sum(W * (Y - mu_hat - r1_out) * weight_tilde) / sum(weight_tilde * W) + sum(mu_hat + r1_out) / nrow(X)
          mean_treat_ped_vec <- W * (Y - mu_hat - r1_out) * weight_tilde + mu_hat + r1_out

          mean_treat_ped_var <- sum((mean_treat_ped_vec - mean_treat_ped) ** 2) / (nrow(X) ** 2)

          c(mean_treat_init,
            mean_treat_init_var,
            mean_treat_ped,
            mean_treat_ped_var)

        }
      doParallel::stopImplicitCluster()
      parallel::stopCluster(cl)
    } else{
      ATE_mat <- matrix(0, nrow = B, ncol = 4)
      for (b in 1:B) {
        # estimate mu using cross-fitting
        loc0_2_group <- sample(1:2, size = length(loc0), replace = TRUE)
        loc1_2_group <- sample(1:2, size = length(loc1), replace = TRUE)
        mu_hat <- rep(0, nrow(X))
        for (foldid in 1:2) {
          loc_main <- c(loc0[loc0_2_group == foldid], loc1[loc1_2_group == foldid])
          loc_aux <- c(loc0[loc0_2_group != foldid], loc1[loc1_2_group != foldid])

          loc_main <- sample(loc_main, size = length(loc_main))
          loc_aux <- sample(loc_aux, size = length(loc_aux))

          X_main <- cbind(rep(1, length(loc_main)), scale(X[loc_main,], scale = FALSE))
          X_aux <- cbind(rep(1, length(loc_aux)), scale(X[loc_aux,], scale = FALSE))

          W_main <- W[loc_main]
          Y_main <- Y[loc_main]
          W_aux <- W[loc_aux]
          Y_aux <- Y[loc_aux]

          n_main <- nrow(X_main)
          p <- ncol(X_main)
          n_aux <- nrow(X_aux)

          r1_aux <- r1_out[loc_aux]

          pi_hat_aux <- pi_hat[loc_aux]
          pi_hat_main <- pi_hat[loc_main]

          Pi_aux <- diag(1 - pi_hat_aux)
          R_aux <- W_aux * (Y_aux - r1_aux) / pi_hat_aux
          Pi_main <- diag(1 - pi_hat_main)
          target <-  t(X_aux) %*% Pi_aux %*% R_aux / n_aux

          ratio_violate <- sum(abs(target)>0.5*sqrt(log(p) / n_main)) / length(target)

          if (ratio_violate < 0.01) {
            mu_hat_main <- rep(0, n_main)
          } else{
            M <- t(X_main) %*% Pi_main / n_main
            weight_term <- ifelse(ratio_violate>0.1,2,3/2)
            mu_hat_main_tmp_ls <-
              quad.prog.Lagr(
                M = M,
                target = target,
                Mr = max(abs(Y)),
                kappa = 0.5,
                verbose = FALSE,
                weight_term = weight_term
              )
            if (mu_hat_main_tmp_ls$is_converge) {
              mu_hat_main <- mu_hat_main_tmp_ls$mu_hat
            } else{
              mu_hat_main <- rep(0, n_main)
            }
          }
          mu_hat[loc_main] <- mu_hat_main
        }
        mean_treat_init <- sum(W * (Y - mu_hat - r1_out) / pi_hat) / sum(W / pi_hat) + sum(mu_hat + r1_out) / nrow(X)
        mean_treat_init_vec <- W * (Y - mu_hat - r1_out) / pi_hat + mu_hat + r1_out
        mean_treat_init_var <- sum((mean_treat_init_vec - mean_treat_init) ** 2) / (nrow(X) ** 2)

        # Double Calibration for pi_tilde (weight_tilde)
        weight_init <- 1 / pi_hat
        if (Y.family == 'binomial') {
          psi_deri <- as.numeric(r1_out  * (1 - r1_out))
        } else{
          psi_deri <- rep(1, nrow(X))
        }
        X_tilde <- cbind(diag(psi_deri) %*% X, mu_hat)
        n <- nrow(X_tilde)
        p <- ncol(X_tilde)
        lambda_vec <- sqrt(log(p) / n) * sqrt(colSums(X_tilde ** 2)) / sqrt(n)

        constr_left <- abs(colSums(diag(as.numeric(W / pi_hat - 1)) %*% X_tilde)) / n
        ratio_violate <- sum(constr_left > lambda_vec)/length(constr_left)

        if (ratio_violate < ratio_violate_dcal) {
          weight_tilde <- weight_init
        } else{
          weight_ls <-
            double_cali.Lagr(X,
                             W,
                             mu_hat,
                             psi_deri,
                             kappa = 0.5,
                             verbose = FALSE)
          weight_tilde <- weight_init
          if (weight_ls$is_converge) {
            weight_tilde[W == 1] <- weight_ls$weight_tilde
          }
        }

        mean_treat_ped <- sum(W * (Y - mu_hat - r1_out) * weight_tilde) / sum(W * weight_tilde) + sum(mu_hat + r1_out) / nrow(X)
        mean_treat_ped_vec <- W * (Y - mu_hat - r1_out) * weight_tilde + mu_hat + r1_out
        mean_treat_ped_var <- sum((mean_treat_ped_vec - mean_treat_ped) ** 2) / (nrow(X) ** 2)

        ATE_mat[b,] <- c(mean_treat_init,
                         mean_treat_init_var,
                         mean_treat_ped,
                         mean_treat_ped_var)

      }
    }

    colnames(ATE_mat) <- c('ATE_sc', 'ATE_sc_var', 'ATE_dc', 'ATE_dc_var')
    return(
      list(
        ATE_sc = mean(ATE_mat[, 1]),
        ATE_sc_var = mean(ATE_mat[, 2]),
        ATE_dc = mean(ATE_mat[, 3]),
        ATE_dc_var = mean(ATE_mat[, 4]),
        ATE_mat = ATE_mat
      )
    )

  }

# Additional functions: -------------------------
## Solve for mu: Quadratic program--------------------
#' @keywords internal
quad.prog.Lagr <- function(M,
                           target,
                           Mr,
                           kappa = 0.5,
                           verbose = FALSE,
                           weight_term = NULL) {
  # min (eta, mu): 1+n dim
  # require(Rmosek)
  prob <- list(sense = 'min')
  prob$c <- rep(0, ncol(M) + 1)

  prob$qobj$i <- 1:(ncol(M) + 1)
  prob$qobj$j <- 1:(ncol(M) + 1)
  prob$qobj$v <- 2 * c(kappa, rep((1 - kappa) / (ncol(M)) ** weight_term, ncol(M)))

  prob$A <-  rbind(cbind(rep(-1, nrow(M)), M),
                   cbind(rep(-1, nrow(M)), -M))
  bvec = c(target,-target)

  prob$bc <- rbind(blc = c(rep(-Inf, 2 * nrow(M))),
                   buc = c(target,-target))
  prob$bx <- rbind(blx = c(0, rep(-Mr, ncol(M))),
                   bux = c(Inf, rep(Mr, ncol(M))))
  prob$iparam <- list(INTPNT_MAX_ITERATIONS = 400)

  prob$dparam <- list(
    INTPNT_QO_TOL_DFEAS = 1.0e-12,
    INTPNT_QO_TOL_PFEAS = 1.0e-12,
    INTPNT_QO_TOL_INFEAS = 1.0e-12,
    INTPNT_QO_TOL_REL_GAP = 1.0e-12
  )

  if (verbose) {
    mosek.out <- Rmosek::mosek(prob)
  } else{
    mosek.out <- Rmosek::mosek(prob, opts = list(verbose = 1))
  }
  mu_hat <- mosek.out$sol$itr$xx[1 + 1:ncol(M)]
  eta_est <- mosek.out$sol$itr$xx[1]

  prob_status <- mosek.out$sol$itr$prosta
  sol_status <- mosek.out$sol$itr$solsta

  if (prob_status == 'PRIMAL_AND_DUAL_FEASIBLE' &
      sol_status == 'OPTIMAL') {
    is_converge <- TRUE
  } else{
    is_converge <- FALSE
    mu_hat = rep(0, ncol(M))
  }

  return(list(
    mu_hat = mu_hat,
    eta_est = eta_est,
    is_converge = is_converge
  ))
}

## Sovle for pi_tilde -----------------------

#' @keywords internal
double_cali.Lagr <-
  function(X,
           W,
           mu_hat,
           psi_deri = NULL,
           kappa = 0.5,
           verbose = FALSE) {
    # require(Rmosek)

    p <- ncol(X)
    n <- nrow(X)
    if (p > n) {
      X_aug <- matrix(runif(n * (p - n),-1, 1), nrow = n)
    } else{
      X_aug <- NULL
    }

    X_tilde <- cbind(diag(psi_deri) %*% X, mu_hat, X_aug)
    X1_tilde <- X_tilde[W == 1,]

    n <- nrow(X_tilde)
    p <- ncol(X_tilde)
    n1 <- nrow(X1_tilde)

    lambda0_vec <- sqrt(colSums(X_tilde ** 2)) / sqrt(n)
    X_tilde_mean <- colMeans(X_tilde) # bar
    prob <- list(sense = 'min')
    prob$c <- c(rep(0, n1), 0)

    prob$qobj$i <- 1:(n1 + 1)
    prob$qobj$j <- 1:(n1 + 1)
    prob$qobj$v <- 2 * c(rep(kappa / n ** 2, n1), 1 - kappa)

    prob$A <- rbind(cbind(t(X1_tilde / n),-lambda0_vec),
                    cbind(-t(X1_tilde / n),-lambda0_vec),
                    c(rep(1, n1), 0))

    prob$bc <- rbind(blc = c(rep(-Inf, 2 * p), n),
                     buc = c(X_tilde_mean,-X_tilde_mean, n)) # constraint condition

    prob$bx <- rbind(blx = c(rep(1e-5, n1), 0),
                     bux = c(rep(1 / 0.05, n1), Inf)) # variable condition

    if (verbose) {
      mosek.out <- Rmosek::mosek(prob)
    } else{
      mosek.out <- Rmosek::mosek(prob, opts = list(verbose = 1))
    }
    weight_tilde <- mosek.out$sol$itr$xx[1:n1]
    eta_pi <- mosek.out$sol$itr$xx[1 + n1]

    prob_status <- mosek.out$sol$itr$prosta
    sol_status <- mosek.out$sol$itr$solsta

    if (prob_status == 'PRIMAL_AND_DUAL_FEASIBLE' &
        sol_status == 'OPTIMAL') {
      is_converge <- TRUE
    } else{
      is_converge <- FALSE
    }

    return(list(
      weight_tilde = weight_tilde,
      eta_pi = eta_pi,
      is_converge = is_converge
    ))

  }
