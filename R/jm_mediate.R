#' Description of jm_mediate()
#' @author Lina Kramer \email{linakramer015@gmail.com}
#' @export
#' @import JMbayes2
#' @import splines
#' @import statmod
#' @importFrom stats formula quantile
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom MASS mvnorm
#' @param jointfit A \code{\link[JMbayes2]{jm}} joint model object.
#' @param ds_surv Survival data frame that was used to fit the relative risk
#' submodel.
#' @param ds_long Long format data frame that was used to fit the
#' longitudinal submodel.
#' @param time_eval Time points at which to compute effects.
#' @param trt_name Name of the treatment variable.
#' @param time_var Time variable in the longitudinal model.
#' @param L_inner Number of draws for the random effects integration.
#' @param n_mcmc How many of the mcmc draws are used.
#' @param seed Used for random draws.
#' @returns A list containing:
#' summary: long format data frame of total, direct, and indirect effects over
#' time;
#' pop_surv: matrix of survival estimates for
#' each posterior draw under
#' s11 = P(do(A^D=1, A^M=1)),
#' s10 = P(do(A^D=0, A^M=1)),
#' s01 = P(do(A^D=1, A^M=0)),
#' s00 = P(do(A^D=0, A^M=0));
#' pop_cif = 1-pop_surv;
#' draws: list of effects for each posterior draw;
#' settings: list of settings used for computation.
#' @examples
#' # First fit a joint model:
#' library(JMbayes2)
#' data("prothro")
#' data("prothros")
#'
#' # longitudinal submodel
#' lmefit <- lme(pro~ ns(time, 3):treat + ns(time,3),
#'              random = ~time | id,
#'              data = prothro,
#'              control = lmeControl(opt = "optim"))
#' # survival submodel
#' coxfit <- coxph(Surv(Time, death)~ treat,
#'             data = prothros, x = TRUE, na.action = na.omit)
#'
#' # jointmodel
#' jointfit <- JMbayes2::jm(coxfit, lmefit, time_var = "time", n_chains =2L,
#'                           n_iter = 40000L, n_burnin = 5000L, n_thin = 3L,
#'                           data_Surv = prothros)
#'
#' # Mediation analysis:
#' med_prothro <- jm_mediate(
#'                   jointfit = jointfit,
#'                   ds_surv = prothros,
#'                   ds_long = prothro,
#'                   time_eval = 0:12,
#'                   trt_name = "treat",
#'                   time_var = "time",
#'                   n_mcmc = 500
#'                   )

# Main function ################################################################

jm_mediate <- function(
    jointfit,
    ds_surv,
    ds_long,
    time_eval,
    trt_name,
    time_var,
    L_inner = 50,
    n_mcmc = NULL,
    seed = 15
) {

  set.seed(seed)

  # Extract MCMC draws
  draws <- extract_mcmc_draws(jointfit, n_mcmc)
  K <- nrow(draws$betas)
  n_re <- (-1 + sqrt(1 + 8 * ncol(draws$D))) / 2 # automatically calculates n of random effects

  #Survival design matrices
  W_list     <- build_scenario_W(jointfit, ds_surv, trt_name)
  n_subjects <- nrow(ds_surv)

  # Quadrature
  quad_list <- precompute_quadrature(
    time_eval,
    knots_bs  = jointfit$control$knots$`1`,
    degree_bs = jointfit$control$Bsplines_degree
  )

  # Precompute all design matrices
  id_var      <- jointfit$model_info$var_names$idVar
  template_df <- ds_long[!duplicated(ds_long[[id_var]]), ]
  template_df <- template_df[match(ds_surv[[id_var]], template_df[[id_var]]), ]


  # get fixed and random effects terms from jointfit
  terms_FE <- jointfit[["model_info"]][["terms"]][["terms_FE_noResp"]][[1]]
  terms_RE <- jointfit[["model_info"]][["terms"]][["terms_RE"]][[1]]

  # precompute X and Z matrices
  dm_list <- precompute_design_matrices(
    quad_list, time_eval, terms_FE, terms_RE,
    template_df, time_var, trt_name, n_subjects
  )

  # Storage
  pop_surv <- array(
    NA,
    dim = c(length(time_eval), K, 4),
    dimnames = list(NULL, NULL, c("s11", "s10", "s00", "s01"))
  )

  message(sprintf("Running: %d posterior draws x %d RE draws x %d subjects x %d time points",
                  K, L_inner, n_subjects, length(time_eval)))

  pb <- txtProgressBar(min = 0, max = K, style = 3)

  for (k in seq_len(K)) { # for each posterior draw

    setTxtProgressBar(pb, k)

    betas_k     <- draws$betas[k, ]
    alphas_k    <- draws$alphas[k, ]
    gammas_k    <- draws$gammas[k, ]
    Bs_gammas_k <- draws$Bs_gammas[k, ]
    D_k         <- reconstruct_D(draws$D[k, ], n_re)

    # Constant across j and l
    Wgamma_trt  <- as.numeric(W_list$W_trt  %*% gammas_k)
    Wgamma_ctrl <- as.numeric(W_list$W_ctrl %*% gammas_k)

    for (j in seq_along(time_eval)) { # for posterior draw for each time point

      if (time_eval[j] == 0) {
        pop_surv[j, k, ] <- 1
        next # needs to skip t0
      }

      qj <- quad_list[[j]]
      dm <- dm_list[[j]]

      # Get baseline hazard and fixed mediator parts
      log_h0   <- as.numeric(qj$B %*% Bs_gammas_k)         # [Q]
      fixed_z1 <- as.numeric(dm$X_z1 %*% betas_k)           # [n_sub * Q]
      fixed_z0 <- as.numeric(dm$X_z0 %*% betas_k)           # [n_sub * Q]

      acc_s11 <- 0
      acc_s10 <- 0
      acc_s00 <- 0
      acc_s01 <- 0

      for (l in seq_len(L_inner)) { # for each posterior draw for each time point for each random effects draw

        a_sampled <- MASS::mvrnorm(n = n_subjects, mu = rep(0, n_re), Sigma = D_k) # sample from multivariate normal distribution with covariance D_k


        # (note: this is just a faster version of obtaining potential mediator values than jm.predict() )
        # Random effects contribution: Z * a, per subject
        random <- rowSums(dm$Z * a_sampled[dm$sub_idx, , drop = FALSE])  # [n_sub * Q]

        # Mediator matrices: [n_subjects x Q]
        M1_mat <- matrix(fixed_z1 + random, nrow = n_subjects, byrow = TRUE)
        M0_mat <- matrix(fixed_z0 + random, nrow = n_subjects, byrow = TRUE)

        # Survival for 4 scenarios
        lp_11 <- sweep(alphas_k * M1_mat, 1, Wgamma_trt, "+")
        lp_11 <- sweep(lp_11, 2, log_h0, "+")

        lp_00 <- sweep(alphas_k * M0_mat, 1, Wgamma_ctrl, "+")
        lp_00 <- sweep(lp_00, 2, log_h0, "+")

        lp_10 <- sweep(alphas_k * M1_mat, 1, Wgamma_ctrl, "+")
        lp_10 <- sweep(lp_10, 2, log_h0, "+")

        lp_01 <- sweep(alphas_k * M0_mat, 1, Wgamma_trt, "+")
        lp_01 <- sweep(lp_01, 2, log_h0, "+")

        w <- qj$w
        acc_s11 <- acc_s11 + mean(exp(-(exp(lp_11) %*% w)))
        acc_s10 <- acc_s10 + mean(exp(-(exp(lp_10) %*% w)))
        acc_s00 <- acc_s00 + mean(exp(-(exp(lp_00) %*% w)))
        acc_s01 <- acc_s01 + mean(exp(-(exp(lp_01) %*% w)))
      }

      pop_surv[j, k, "s11"] <- acc_s11 / L_inner # means over inner RE draws
      pop_surv[j, k, "s10"] <- acc_s10 / L_inner
      pop_surv[j, k, "s00"] <- acc_s00 / L_inner
      pop_surv[j, k, "s01"] <- acc_s01 / L_inner
    }
  }

  close(pb)

  results <- summarise_effects(pop_surv, time_eval)

  list(
    summary  = results$summary,
    pop_surv = pop_surv,
    pop_cif  = results$pop_cif,
    draws    = results$draws,
    settings = list(K = K, L_inner = L_inner, n_subjects = n_subjects,
                    n_times = length(time_eval), n_re = n_re)
  )
}




# Helper functions #############################################################

# reconstructs the covariance matrix of random effects
# inputs:
# d_vec: vector of random effects covariances D, extracted from joint model fit in JMbayes2
# n_re: number of random effects
reconstruct_D <- function(d_vec, n_re) {
  D <- matrix(0, n_re, n_re)
  D[lower.tri(D, diag = TRUE)] <- d_vec
  D[upper.tri(D)] <- t(D)[upper.tri(D)]
  D # returns D matrix
}

# builds the survival design matrices for intervention and no intervention on surival path
# inputs:
# jointfit: joint model fir in JMbayes2
# ds_surv: survival data frame
# trt_name: treatment variable name in survival submodel
build_scenario_W <- function(jointfit, ds_surv, trt_name) {
  form_s <- formula(jointfit$model_info$terms$terms_Surv_noResp)
  W <- model.matrix(form_s, data = ds_surv)
  W <- W[, -1, drop = FALSE]

  trt_col <- which(colnames(W) == trt_name)
  if (length(trt_col) == 0) {
    trt_col <- grep(paste0("^", trt_name), colnames(W)) # find correct col
    if (length(trt_col) > 1) {
      stop("Multiple columns match '", trt_name, "': ",
           paste(colnames(W)[trt_col], collapse = ", "),
           ". Treatment factor must have exactly 2 levels.")
    }
    if (length(trt_col) == 0) {
      stop("Treatment variable '", trt_name, "' not found in survival design matrix. ",
           "Available columns: ", paste(colnames(W), collapse = ", "))
    }
    message(sprintf("Using column '%s' for treatment.", colnames(W)[trt_col]))
  }
  W_trt <- W_ctrl <- W
  W_trt[, trt_col]  <- 1
  W_ctrl[, trt_col] <- 0
  list(W_trt = W_trt, W_ctrl = W_ctrl) # returns list of two design matrices
}

# extracts all mcmc_draws (when n_mcmc is NULL) or randomly samples n_mcmc draws
extract_mcmc_draws <- function(jointfit, n_mcmc = NULL) {
  draws <- list(
    betas     = do.call(rbind, jointfit$mcmc$betas1), # do.call joins all mcmc chains
    alphas    = do.call(rbind, jointfit$mcmc$alphas),
    gammas    = do.call(rbind, jointfit$mcmc$gammas),
    Bs_gammas = do.call(rbind, jointfit$mcmc$bs_gammas),
    D         = do.call(rbind, jointfit$mcmc$D)
  )
  n_total <- nrow(draws$betas) # automatically counts how many total mcmc draws
  if (!is.null(n_mcmc) && n_mcmc < n_total) {
    idx <- sample(n_total, n_mcmc)
    draws <- lapply(draws, function(x) x[idx, , drop = FALSE])
  }
  draws # returns mcmc draws of betas, alphas, gammas, BS_gammas, and D
}

# precompute the 15 point Gauss-Legendre quadrature that is used for the approximation of the survival function
# inputs:
# time_eval: timepoints at which to evaluate (input by user)
# knots_bs: B-spline knots, extracted from joint model object
# degree_bs: B-spline degree, extracted from joint model object
precompute_quadrature <- function(time_eval, knots_bs, degree_bs) {
  gl <- statmod::gauss.quad(15, kind = "legendre")
  lapply(time_eval, function(t_j) {
    if (t_j == 0) return(NULL)
    u <- (t_j / 2) * (gl$nodes + 1)
    w <- (t_j / 2) * gl$weights
    B <- splines::splineDesign(knots = knots_bs, x = u,
                               ord = degree_bs + 1, outer.ok = TRUE)
    list(u = u, w = w, B = B) # returns list of length(time_eval)-1 including 15 quad points u, 15  weights w, and the B-spline design matrix
  })
}

# precompute fixed and random design matrices at quadrature points
# for each time point and both treatment levels.
# inputs:
# quad_list: output of precompute_quadrature()
# time_eval
# terms_FE & terms_RE: fixed effects and random effects terms extracted from joint model object
# template_df: dataset including one row per subject extracted from ds_long
# time_var: time variable used in joint model
# trt_var: treatment variable used in longitudinal submodel
# n_subjects: number of individual subjects in data set

precompute_design_matrices <- function(quad_list, time_eval, terms_FE, terms_RE,
                                       template_df, time_var, trt_var, n_subjects) {

  n_quad <- 15L
  sub_idx <- rep(seq_len(n_subjects), each = n_quad)

  # Detect treatment variable type and figure out the two levels
  trt_col_orig <- template_df[[trt_var]]
  if (is.null(trt_col_orig)) {
    stop("Treatment variable '", trt_var, "' not found in ds_long")
  }

  if (is.factor(trt_col_orig)) {
    trt_levels <- levels(trt_col_orig)
    if (length(trt_levels) != 2) {
      stop("Treatment factor must have exactly 2 levels, found: ",
           paste(trt_levels, collapse = ", "))
    }
    trt_value_0 <- factor(trt_levels[1], levels = trt_levels)
    trt_value_1 <- factor(trt_levels[2], levels = trt_levels)
    message(sprintf("Treatment '%s' is a factor: control = '%s', treatment = '%s'",
                    trt_var, trt_levels[1], trt_levels[2]))
  } else if (is.numeric(trt_col_orig) || is.logical(trt_col_orig)) {
    trt_value_0 <- 0
    trt_value_1 <- 1
  } else if (is.character(trt_col_orig)) {
    trt_levels <- sort(unique(trt_col_orig))
    if (length(trt_levels) != 2) {
      stop("Treatment character must have exactly 2 unique values")
    }
    trt_value_0 <- trt_levels[1]
    trt_value_1 <- trt_levels[2]
  } else {
    stop("Unsupported treatment variable type: ", class(trt_col_orig))
  }

  lapply(seq_along(time_eval), function(j) {

    if (time_eval[j] == 0) return(NULL)

    qj <- quad_list[[j]]

    df_quad <- template_df[rep(seq_len(n_subjects), each = n_quad), ]
    df_quad[[time_var]] <- rep(qj$u, times = n_subjects)

    Z_quad <- model.matrix(terms_RE, data = df_quad)

    df_quad[[trt_var]] <- trt_value_1
    X_z1 <- model.matrix(terms_FE, data = df_quad)

    df_quad[[trt_var]] <- trt_value_0
    X_z0 <- model.matrix(terms_FE, data = df_quad)

    list(X_z1 = X_z1, X_z0 = X_z0, Z = Z_quad, sub_idx = sub_idx)
  })
}

# calculates NIE, NDE, and TE based on cumulative incidence functions pop_surv at each time_eval
summarise_effects <- function(pop_surv, time_eval) {
  pop_cif <- 1 - pop_surv

  nie_draws <- pop_cif[, , "s11"] - pop_cif[, , "s01"] # checked
  nde_draws <- pop_cif[, , "s01"] - pop_cif[, , "s00"]
  te_draws  <- pop_cif[, , "s11"] - pop_cif[, , "s00"]

  make_summary <- function(draws, label) {
    data.frame(
      time     = time_eval,
      effect   = label,
      estimate = rowMeans(draws),
      lower    = apply(draws, 1, quantile, 0.025), # Bayesian CIs
      upper    = apply(draws, 1, quantile, 0.975)
    )
  }

  summary_df <- rbind(
    make_summary(te_draws, "Total"),
    make_summary(nde_draws, "Direct"),
    make_summary(nie_draws, "Indirect")
  )
  summary_df$effect <- factor(summary_df$effect, levels = c("Total", "Direct", "Indirect"))

  list(summary = summary_df, pop_cif = pop_cif,
       draws = list(NIE = nie_draws, NDE = nde_draws, TE = te_draws))
}

