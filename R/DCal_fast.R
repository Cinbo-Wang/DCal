#' Efficient Doubly-Calibrated Estimation of E[Y(1)] for Observational Data
#'
#' Implements a computationally efficient version of Double Calibration for estimating expected outcomes under treatment.
#' By minimizing L2-norm of inverse propensity weights (1/proensity score) during calibration instead of L1-norm of coefficients,
#' this function achieves faster computation than \code{DCal.mean_treat} while maintaining root-n consistency under
#' minimal sparsity conditions. Uses cross-fitting to avoid overfitting.
#'
#' @param X An \eqn{n \times p} numeric matrix of high-dimensional covariates.
#' @param Y Numeric vector of length \eqn{n} containing observed outcomes.
#' @param W Numeric binary vector of length \eqn{n} (0 = control, 1 = treatment) for treatment assignments.
#' @param Y.family Character specifying outcome type: \code{"gaussian"} for continuous or \code{"binomial"} for binary outcomes.
#' @param B Number of random sample splits for cross-fitting. Default: \code{3}.
#' @param K Number of folds per split for cross-fitting. Default: \code{2} (recommended).
#' @param is.scale Logical indicating whether to standardize \code{X}. Default: \code{FALSE}.
#' @param r1_init Optional initial estimates for \eqn{E[Y|X,W=1]}. Numeric vector of length \eqn{n}.
#' @param pi_init Optional initial propensity score estimates \eqn{\hat\pi(X) = P(W=1|X)}. Numeric vector of length \eqn{n}.
#' @param gamma_init Optional initial vector for calibration coefficients. Default: \code{NULL}.
#' @param alpha Elastic net mixing parameter (0 = ridge, 1 = lasso) for initial nuisance function estimation.
#' @param is.parallel Logical indicating whether to parallelize across sample splits. Default: \code{FALSE}.
#' @param core_num Number of cores for parallel computation when \code{is.parallel = TRUE}.
#'
#' @return A list containing:
#' \item{ATE_vec}{Numeric vector of point estimates averaged over splits: [Single-Calibration, Double-Calibration] for E[Y(1)]}
#' \item{var_ATE_vec}{Numeric vector of variance estimates averaged over splits: [Single-Calibration, Double-Calibration]}
#' \item{ATE_mat}{\eqn{B \times 2} matrix of point estimates per split (Column 1: Single-Calibration, Column 2: Double-Calibration)}
#' \item{var_ATE_mat}{\eqn{B \times 2} matrix of variance estimates per split}
#'
#' @author Xinbo Wang
#' @references
#' Liu, L., Wang, X., Liu, L. and Wang, Y. (2023) \emph{Root-n consistent semiparametric learning with high-dimensional nuisance functions under minimal sparsity. arXiv preprint, arXiv:2305.04174}, \doi{10.48550/arXiv.2305.04174}.\cr
#' Wang, Y. and Shah, R. D. (2025) \emph{Debiased inverse propensity score weighting for estimation of average treatment effects with high-dimensional confounders. The Annals of Statistics, Vol. 52(5), 1978-2003}, \doi{10.1214/24-AOS2409}.


#' @examples
#' p = 400
#' s_or = 10
#' n = 200
#' rho = 0.9
#' rd_num = 1
#' Sigma_X <- matrix(0, p, p)
#' for (i in 1:p) {
#'   for (j in 1:p) {
#'     Sigma_X[i, j] <- rho ** abs(i - j)
#'   }
#' }
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma_X)
#'
#' set.seed(rd_num)
#'
#' # Sparse PS
#' gamma_true <- rep(0, p)
#' act_loc <- 1:s_or
#' gamma_true[act_loc] <- runif(length(act_loc), 1, 2)
#' gamma_true <- gamma_true / norm(gamma_true, type = '2')
#'
#' U <- X %*% gamma_true
#' pi_W <- 1 / (1 + exp(-U))
#' W <- rbinom(n, size = 1, p = pi_W)
#'
#' # nonlinear dense OR
#' require(splines)
#' beta_true <- rep(0, p)
#' act_loc <- 1:p
#' for (j in act_loc) {
#'   beta_true[j] <- 1 / j
#' }
#' beta_true <-  beta_true / norm(beta_true, type = '2')
#'
#' X2 <- 2 / (1 + exp(-X[, 2]))
#' X3 <- exp(X[, 3] / 2)
#' X4 <- X[, 4] / (1 + exp(X[, 3]))
#' X5 <- X[, 4] * X[, 5] / 10
#'
#' non_lin <- scale(cbind(X2, X3, X4, X5)) %*% c(1, -1 / 2, 1 / 3, -1 / 4) +
#'     bs(X[, 1], df = 100) %*% c(1 / (1:100))
#' term_cm <- (abs(non_lin) + 0.05) ** (-1) + X %*% beta_true
#' potential_outcome_treat <-  term_cm - 1
#' potential_outcome_control <- term_cm
#' Y <- potential_outcome_treat * W + potential_outcome_control * (1 - W) + rnorm(n, 0, 1)
#' tau_treat_g <- mean(potential_outcome_treat)
#'
#'
#' mean_treat_dcal_ls <- DCal_fast.mean_treat(
#'   X,
#'   Y,
#'   W,
#'   Y.family = 'gaussian',
#'   B = 3,
#'   K = 2,
#'   is.scale = FALSE,
#'   alpha = 0.9,
#'   is.parallel = TRUE,
#'   core_num = 2
#' )
#' ate_vec <- mean_treat_dcal_ls$ATE_vec
#' sd_vec <- sqrt(mean_treat_dcal_ls$var_ATE_vec)
#' CI_lb_vec <- ate_vec - qnorm(0.975) * sd_vec
#' CI_ub_vec <- ate_vec + qnorm(0.975) * sd_vec
#' ate_vec
#' CI_lb_vec
#' CI_ub_vec

#' @export DCal_fast.mean_treat
#' @import parallel
#' @import glmnet
#' @import foreach
#' @import doParallel
#' @import stats
DCal_fast.mean_treat <- function(X,
                                 Y,
                                 W,
                                 Y.family = 'gaussian',
                                 B = 3,
                                 K = 2,
                                 is.scale = FALSE,
                                 r1_init = NULL,
                                 pi_init = NULL,
                                 gamma_init = NULL,
                                 alpha = 0.9,
                                 is.parallel = FALSE,
                                 core_num = NULL) {
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
      glmnet::cv.glmnet(X[W == 1, ],
                Y[W == 1],
                family = Y.family,
                alpha = alpha)
    r1_out <-
      predict(fit.out.treated, newx = X, type = 'response')
  } else{
    r1_out <- r1_init
  }

  if (is.null(pi_init)) {
    fit.prop <-
      glmnet::cv.glmnet(X,
                W,
                family = "binomial",
                alpha = alpha)
    pi_hat <- predict(fit.prop, newx = X, type = 'response')
    pi_hat <- pmax(pmin(pi_hat, 0.99), 0.01)
    gamma_hat <- as.numeric(coef(fit.prop))

  } else{
    pi_hat <-  pmax(pmin(pi_init, 0.99), 0.01)
    gamma_hat <- gamma_init
  }

  loc0 <- which(W == 0)
  loc1 <- which(W == 1)
  if (is.parallel) {
    # require(doParallel);require(foreach)
    type <- ifelse(.Platform$OS.type == 'windows', 'PSOCK', 'FORK')
    core_num <- ifelse(!is.null(core_num),
                       core_num,
                       ifelse(.Platform$OS.type == 'windows', 4, min(25, B)))
    cl <- parallel::makeCluster(core_num, type = type)
    doParallel::registerDoParallel(cl)
    ATE_mat_full <-
      foreach::foreach(
        b = 1:B,
        .combine = 'rbind',
        .export = c('quad.prog.mu', 'double_cali_pi')
      ) %dopar% {
        set.seed(b)

        # estimate mu using cross-fitting
        mu_hat <- rep(0, nrow(X))
        loc0_2_group <- sample(1:K, size = length(loc0), replace = TRUE)
        loc1_2_group <- sample(1:K, size = length(loc1), replace = TRUE)

        for (foldid in 1:K) {
          loc_main <- c(loc0[loc0_2_group == foldid], loc1[loc1_2_group == foldid])
          loc_aux <- c(loc0[loc0_2_group != foldid], loc1[loc1_2_group != foldid])

          loc_main <- sample(loc_main, size = length(loc_main))
          loc_aux <- sample(loc_aux, size = length(loc_aux))

          X_main <- cbind(1, scale(X[loc_main, ], scale = FALSE))
          X_aux <- cbind(1, scale(X[loc_aux, ], scale = FALSE))

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

          ratio_violate <- sum(abs(target) > 0.1 * sqrt(log(p) / n_main)) / length(target)
          eta_r <- 0.1 * sqrt(log(p) / n_main)
          if (sum(abs(target) > eta_r) == 0) {
            mu_hat_main <- rep(0, n_main)
          } else{
            for (ratio in seq(0.1, 5, 0.1)) {
              eta_r <- ratio * sqrt(log(p) / n_main)
              M <- t(X_main) %*% Pi_main / n_main

              lb <- min((Y - r1_out)[loc_main])
              ub <- max((Y - r1_out)[loc_main])
              mu_hat_main_tmp_ls <-
                quad.prog.mu(
                  M = M,
                  target = target,
                  lb = lb,
                  ub = ub,
                  eta_r = eta_r,
                  verbose = FALSE
                )
              if (mu_hat_main_tmp_ls$is_converge) {
                mu_hat_main <- mu_hat_main_tmp_ls$mu_hat
                break
              } else{
                mu_hat_main <- rep(0, n_main)
              }
            }
          }

          mu_hat[loc_main] <- mu_hat_main
        }
        mean_treat_init <-
          sum(W * (Y - mu_hat - r1_out) / pi_hat) / sum(W / pi_hat) + sum(mu_hat + r1_out) / nrow(X)
        mean_treat_init_vec <- W * (Y - mu_hat - r1_out) / pi_hat + mu_hat + r1_out
        mean_treat_init_var <- sum((mean_treat_init_vec - mean_treat_init) ** 2) / (nrow(X) ** 2)

        # Double Calibration for pi_tilde
        weight_init <- 1 / pi_hat
        if (Y.family == 'binomial') {
          psi_deri <- as.numeric(r1_out  * (1 - r1_out))
        } else{
          psi_deri <- rep(1, nrow(X))
        }
        X_tilde <- cbind(diag(psi_deri) %*% X, mu_hat)
        n <- nrow(X_tilde)

        lambda_vec <- 0.1 * sqrt(log(ncol(X) + 1) / n) * sqrt(colSums(X_tilde ** 2)) / sqrt(n)
        constr_left <- abs(colSums(diag(as.numeric(W / pi_hat - 1)) %*% X_tilde)) / n
        ratio_violate <- sum(constr_left > lambda_vec) / length(constr_left)
        pi_tilde_init <- 1 / pi_hat
        if (ratio_violate > 0) {
          # Step 2: Find the minimal eta_pi with feasible solution by directly minimizing pi
          weight_tilde <- NULL
          for (ratio in seq(0.1, 5, 0.1)) {
            eta_pi <- ratio * sqrt(log(ncol(X) + 1) / nrow(X))
            result_dcal_init_ls <- double_cali_pi(X, W, mu_hat, psi_deri, eta_pi)


            if (result_dcal_init_ls$is_converge) {
              weight_tilde <- result_dcal_init_ls$weight_tilde
              break
            }
          }
          if (!is.null(weight_tilde)) {
            pi_tilde_init[W == 1] <- 1 / weight_tilde
          }

        }

        mean_treat_ped <- sum((Y - mu_hat - r1_out) * W / pi_tilde_init) / sum(W /
                                                                                 pi_tilde_init) + sum(mu_hat + r1_out) / nrow(X)
        mean_treat_ped_vec <- (Y - mu_hat - r1_out) * W / pi_tilde_init + mu_hat + r1_out
        mean_treat_ped_var <- sum((mean_treat_ped_vec - mean_treat_ped) ** 2) / (nrow(X) ** 2)



        c(mean_treat_init,
          mean_treat_ped,
          mean_treat_init_var,
          mean_treat_ped_var)

      }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)


    ATE_mat <- ATE_mat_full[, 1:2]
    var_ATE_mat <- ATE_mat_full[, 3:4]

  } else{
    ATE_mat <- var_ATE_mat <- matrix(0, nrow = B, ncol = 2)
    for (b in 1:B) {
      set.seed(b)
      # estimate mu using cross-fitting
      loc0_2_group <- sample(1:K, size = length(loc0), replace = TRUE)
      loc1_2_group <- sample(1:K, size = length(loc1), replace = TRUE)
      mu_hat <- rep(0, nrow(X))
      for (foldid in 1:K) {
        loc_main <- c(loc0[loc0_2_group == foldid], loc1[loc1_2_group == foldid])
        loc_aux <- c(loc0[loc0_2_group != foldid], loc1[loc1_2_group != foldid])

        loc_main <- sample(loc_main, size = length(loc_main))
        loc_aux <- sample(loc_aux, size = length(loc_aux))

        X_main <- cbind(rep(1, length(loc_main)), scale(X[loc_main, ], scale = FALSE))
        X_aux <- cbind(rep(1, length(loc_aux)), scale(X[loc_aux, ], scale = FALSE))

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

        ratio_violate <- sum(abs(target) > 0.1 * sqrt(log(p) / n_main)) / length(target)
        eta_r <- 0.1 * sqrt(log(p) / n_main)
        if (sum(abs(target) > eta_r) == 0) {
          mu_hat_main <- rep(0, n_main)
        } else{
          for (ratio in seq(0.1, 5, 0.1)) {
            eta_r <- ratio * sqrt(log(p) / n_main)
            M <- t(X_main) %*% Pi_main / n_main

            lb <- min((Y - r1_out)[loc_main])
            ub <- max((Y - r1_out)[loc_main])
            mu_hat_main_tmp_ls <-
              quad.prog.mu(
                M = M,
                target = target,
                lb = lb,
                ub = ub,
                eta_r = eta_r,
                verbose = FALSE
              )
            if (mu_hat_main_tmp_ls$is_converge) {
              mu_hat_main <- mu_hat_main_tmp_ls$mu_hat
              break
            } else{
              mu_hat_main <- rep(0, n_main)
            }
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
      lambda_vec <- 0.1 * sqrt(log(ncol(X) + 1) / n) * sqrt(colSums(X_tilde ** 2)) / sqrt(n)
      constr_left <- abs(colSums(diag(as.numeric(W / pi_hat - 1)) %*% X_tilde)) / n
      ratio_violate <- sum(constr_left > lambda_vec) / length(constr_left)
      pi_tilde_init <- 1 / pi_hat
      if (ratio_violate > 0) {
        # Step 2: Find the minimal eta_pi with feasible solution by directly minimizing pi
        weight_tilde <- NULL
        for (ratio in seq(0.1, 5, 0.1)) {
          eta_pi <- ratio * sqrt(log(ncol(X) + 1) / nrow(X))
          weight_init <- 1 / pi_hat[W == 1]
          result_dcal_init_ls <- double_cali_pi(X, W, mu_hat, psi_deri, eta_pi)

          if (result_dcal_init_ls$is_converge) {
            weight_tilde <- result_dcal_init_ls$weight_tilde
            break
          }
        }
        if (!is.null(weight_tilde)) {
          pi_tilde_init[W == 1] <- 1 / weight_tilde
        }

      }

      mean_treat_ped <- sum((Y - mu_hat - r1_out) * W / pi_tilde_init) / sum(W /
                                                                               pi_tilde_init) + sum(mu_hat + r1_out) / nrow(X)
      mean_treat_ped_vec <- (Y - mu_hat - r1_out) * W / pi_tilde_init + mu_hat + r1_out
      mean_treat_ped_var <- sum((mean_treat_ped_vec - mean_treat_ped) ** 2) / (nrow(X) ** 2)

      ATE_mat[b, ] <- c(mean_treat_init, mean_treat_ped)
      var_ATE_mat[b, ] <- c(mean_treat_init_var, mean_treat_ped_var)

    }
  }

  colnames(ATE_mat) <- colnames(var_ATE_mat) <- c('SCal', 'DCal')

  return(
    list(
      ATE_vec = colMeans(ATE_mat),
      var_ATE_vec = colMeans(var_ATE_mat),
      ATE_mat = ATE_mat,
      var_ATE_mat = var_ATE_mat
    )
  )

}



# Additional functions: -------------------------
## Solve for mu: Quadratic program--------------------
#' @keywords internal
quad.prog.mu <- function(M,
                      target,
                      lb,
                      ub,
                      eta_r = NULL,
                      verbose = FALSE) {
  # min (mu[n]): n dim
  # require(Rmosek)
  prob <- list(sense = 'min')
  prob$c <- rep(0, ncol(M))

  prob$qobj$i <- 1:ncol(M)
  prob$qobj$j <- 1:ncol(M)
  prob$qobj$v <- 2 * rep(1 / ncol(M), ncol(M))

  prob$A <-  rbind(cbind(M), cbind(-M))

  p <- nrow(M)
  n <- ncol(M)
  bvec = c(target + eta_r, -target + eta_r)

  prob$bc <- rbind(blc = c(rep(-Inf, 2 * nrow(M))), buc = bvec)

  prob$bx <- rbind(blx = c(rep(lb, ncol(M))), bux = c(rep(ub, ncol(M))))

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
  mu_hat <- mosek.out$sol$itr$xx[1:ncol(M)]

  prob_status <- mosek.out$sol$itr$prosta
  sol_status <- mosek.out$sol$itr$solsta

  if (prob_status == 'PRIMAL_AND_DUAL_FEASIBLE' &
      sol_status == 'OPTIMAL') {
    is_converge <- TRUE
  } else{
    is_converge <- FALSE
    mu_hat = rep(0, ncol(M))
  }

  return(list(mu_hat = mu_hat, is_converge = is_converge))
}

## Sovle for pi_tilde -----------------------

#' @keywords internal
double_cali_pi <- function(X,
                           W,
                           mu_hat,
                           psi_deri = NULL,
                           eta_pi = NULL,
                           verbose = FALSE) {
  # require(Rmosek)
  Xaug <- cbind(1, X)
  p <- ncol(Xaug)
  n <- nrow(Xaug)
  if(n <= p){
    X_tilde <- Xaug
  }else{
    X_tilde <- cbind(Xaug, matrix(runif(n * (n - p),-1, 1), nrow = n))
  }

  if(all(unique(psi_deri) == 1)){
    # gaussian
    if(n <= p){
      X_cmb <- cbind(Xaug, mu_hat)
    }else{
      X_cmb <- cbind(Xaug, mu_hat, X_tilde[, -c(1:p)])
    }
  }else{
    X_cmb <- cbind(diag(psi_deri) %*% Xaug, mu_hat, X_tilde)
  }

  num_cond <- ncol(X_cmb)
  X1_cmb <- X_cmb[W == 1,]
  n1 <- nrow(X1_cmb)

  lambda0_vec <- sqrt(colSums(X_cmb ** 2)) / sqrt(n) * eta_pi
  X_cmb_mean <- colMeans(X_cmb)


  prob <- list(sense = 'min')
  prob$c <- c(rep(0, n1))
  prob$qobj$i <- 1:(n1)
  prob$qobj$j <- 1:(n1)
  prob$qobj$v <- 2 * rep(1 / n1, n1)
  prob$A <- rbind(cbind(t(X1_cmb / n)), cbind(-t(X1_cmb / n)), rep(1, n1))
  prob$bc <- rbind(
    blc = c(rep(-Inf, 2 * num_cond), n),
    buc = c(
      X_cmb_mean + lambda0_vec - 1e-4,
      -X_cmb_mean + lambda0_vec  - 1e-4,
      n
    )
  )
  prob$bx <- rbind(blx = c(rep(1e-5, n1)),
                   bux = c(rep(1 / 0.01, n1)))

  if (verbose) {
    mosek.out <- Rmosek::mosek(prob)
  } else{
    mosek.out <- Rmosek::mosek(prob, opts = list(verbose = 1))
  }
  weight_tilde <- mosek.out$sol$itr$xx[1:n1]

  prob_status <- mosek.out$sol$itr$prosta
  sol_status <- mosek.out$sol$itr$solsta
  if (prob_status == 'PRIMAL_AND_DUAL_FEASIBLE' &
      sol_status == 'OPTIMAL') {
    is_converge <- TRUE
  } else{
    is_converge <- FALSE
  }
  return(list(weight_tilde = weight_tilde, is_converge = is_converge))
}
