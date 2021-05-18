library(lavaan)
install.packages("MASS")
library(MASS)

# Data generation process --------------------------------------------------------
# nobs: number of observations
# nind: number of indicators per latent variable
dgp <- function(beta = c(0, 1),
                nobs = 100, 
                xi_var = 1,
                xi_mean = 0,
                zeta_var = 1,
                zeta_mean = 0, 
                delta_mean = 0,
                dv_lower = 0.14,
                dv_upper = 4.14,
                epsilon_mean = 0,
                ev_lower = 0.14,
                ev_upper = 4.14,
                lambda_x_lower = 0.5,
                lambda_x_upper = 1.5,
                lambda_y_lower = 0.5,
                lambda_y_upper = 1.5,
                nind = 5) {
  # Latent variables
  xi <- rnorm(nobs,
              xi_mean,
              sqrt(xi_var)) # true x (exogenous latent variable)
  zeta <- rnorm(nobs,
                zeta_mean,
                sqrt(zeta_var))  # residual of true y
  eta <- beta[1] + beta[2] * xi + zeta  # true y (endogenous latent variable)
  
  # Generate error variances and factor loadings
  delta_var <- matrix(NA, nrow = nobs, ncol = nind) 
  epsilon_var <- matrix(NA, nrow = nobs, ncol = nind) 
  lambda_x <- matrix(NA, nrow = nobs, ncol = nind) 
  lambda_y <- matrix(NA, nrow = nobs, ncol = nind) 
  
  for (i in 1:nind) {
    delta_var[, i] <- runif(1, dv_lower, dv_upper)
    epsilon_var[,i] <- runif(1, ev_lower, ev_upper)
    lambda_x[, i] <- runif(1, lambda_x_lower, lambda_x_upper)
    lambda_y[, i] <- runif(1, lambda_y_lower, lambda_y_upper)
  }
  
  # Indicators and their measurement errors
  x <- matrix(NA, nrow = nobs, ncol = nind) 
  y <- matrix(NA, nrow = nobs, ncol = nind)
  delta <- matrix(NA, nrow = nobs, ncol = nind)
  epsilon <- matrix(NA, nrow = nobs, ncol = nind)
  for (i in 1:nind) {
    delta[, i] <-rnorm(nobs,
                       delta_mean,
                       sqrt(delta_var[1, i]))  # measurement error for x_i
    epsilon[, i] <- rnorm(nobs,
                          epsilon_mean, 
                          sqrt(epsilon_var[1, i]))  # measurement error for y_i
    x[, i] <- xi * lambda_x[, i] + delta[, i]  # measured x_i
    y[, i] <- eta * lambda_y[, i] + epsilon[, i]  # measured y_i
  }
  
  # Calculate rowsums used for McDonald's Omega
  x_sum <- apply(x, 1, sum)  # rowsum of all indicators for xi
  y_sum <- apply(y, 1, sum)  # rowsum of all indicators for eta
  lambda_x_sum <- apply(lambda_x, 1, sum)
  lambda_y_sum <- apply(lambda_y, 1, sum)
  var_delta_sum <- apply(delta_var, 1, sum)
  var_epsilon_sum <- apply(epsilon_var, 1, sum)
  
  # Calculate reliability (McDonald's Omega) for x- and y-variables
  rel_x <- ((lambda_x_sum)^2)/(((lambda_x_sum)^2) + var_delta_sum)
  rel_y <- ((lambda_y_sum)^2)/(((lambda_y_sum)^2) + var_epsilon_sum)
  
  # Output of dgp() function with column names
  gen_dat <- data.frame(xi, lambda_x, x, x_sum, rel_x,
                        delta, delta_var,
                        zeta, eta, 
                        lambda_y, y, y_sum, rel_y,
                        epsilon, epsilon_var)  # store generated data as data frame
  delta_colnames <- paste("delta", 1:nind, sep="")
  x_colnames <- paste("x", 1:nind, sep="")
  epsilon_colnames <- paste("epsilon", 1:nind, sep="")
  y_colnames <- paste("y", 1:nind, sep="")
  delta_var_colnames <- paste("delta_var", 1:nind, sep="")
  epsilon_var_colnames <- paste("epsilon_var", 1:nind, sep="")
  lambda_x_colnames <- paste("lambda_x", 1:nind, sep="")
  lambda_y_colnames <- paste("lambda_y", 1:nind, sep="")
  names(gen_dat) <- c("xi", 
                      lambda_x_colnames, x_colnames, "x_sum", "rel_x",
                      delta_colnames, delta_var_colnames, 
                      "zeta", "eta", 
                      lambda_y_colnames, y_colnames, "y_sum", "rel_y",
                      epsilon_colnames, epsilon_var_colnames)
  
  # return data frame with all generated data
  return(gen_dat)
}


# Single iteration ----------------------------------------------------
iterate <- function(dv_lower = 0.14,
                    dv_upper = 4.14,
                    ev_lower = 0.14,
                    ev_upper = 4.14,
                    lambda_x_lower = 0.5,
                    lambda_x_upper = 1.5,
                    lambda_y_lower = 0.5,
                    lambda_y_upper = 1.5,
                    nind = 5,
                    cfa_model =
                      "# latent variable definitions (5 indicators)
                    xi =~ x1 + x2 + x3 + x4 + x5
                    eta =~ y1 + y2 + y3 + y4 + y5
                    # regression (5 indicators)
                    eta ~ xi",  # WARNING: make sure to adjust cfa_model if nind != 5
                    ...) {
  # make sure lavaan is installed
  require(lavaan)
  
  # generate data using dgp()
  dat <- dgp(dv_lower = dv_lower, 
             dv_upper = dv_upper,
             ev_lower = ev_lower,
             ev_upper = ev_upper,
             lambda_x_lower = lambda_x_lower,
             lambda_x_upper = lambda_x_upper,
             lambda_y_lower = lambda_y_lower,
             lambda_y_upper = lambda_y_upper,
             nind = nind,
             ...)
  
  # create empty variable used for potential warnings
  warn <- NA
  
  
  # fit data with linear regression
  lm_fit <- summary(lm(y_sum ~ x_sum, data = dat))  # regression model based on rowsums of all indicators
  
  # create data frame with observed variables only
  dat_xy <- cbind(dat[, (nind + 2):(2 * nind + 1)], 
                  dat[, (5 * nind + 6):(6 * nind + 5)])
  
  # create empty data frame for CFA coefficients (used when delta_var & epsilon_var == 0)
  cfa_coefs <- data.frame(
    lhs = "eta",
    op = "~",
    rhs = "xi",
    est = NA,
    se = NA,
    z = NA,
    pvalue = NA,
    ci.lower = NA,
    ci.upper = NA,
    stringsAsFactors = FALSE
  )
  
  # create empty list to store CFA parameters in
  cfa_res <- list(param = cfa_coefs, warn = NA, w = FALSE)
  
  # fit data with SEM and capture output
  cfa_res <- tryCatch({
    list(param = parameterEstimates(cfa(model = cfa_model, data = dat_xy)), warn = NA, w = FALSE)
  }, warning = function(w) {
    list(param = parameterEstimates(cfa(model = cfa_model, data = dat_xy)), warn = w, w = TRUE)
  }, error = function(e) {
    list(param = cfa_coefs, warn = e, w = TRUE)
  }
  )
  
  # extract relevant information from linear regression model
  lm_coefs <- lm_fit$coefficients  # extract info about regression coefficient only
  
  # extract info about regression estimate from SEM
  cfa_coefs <- cfa_res$param[(cfa_res$param$lhs == "eta")
                             & (cfa_res$param$rhs == "xi"), ]
  
  # return info about both models
  return(list(LM = lm_coefs, SEM = cfa_coefs, 
              warn = cfa_res$warn, w = cfa_res$w, 
              dat = dat))
}


# Simulation function -------------------------------------------------------
# Creates nind indicators of two latent variables, adds random measurement errors delta and epsilon
# to the indicators, fits linear regression model (LM) and structural equation model (SEM) to data
# and iterates niter times through all steps. Allows for comparison of LM and SEM with a varying
# number of observations, measurement error variances and number of indicators.
simulate <- function(niter = 10,  # number of iterations per treatment combination
                     beta = c(0, 1),
                     nobs = 500,
                     xi_var = 1,
                     xi_mean = 0,
                     zeta_var = 1,
                     zeta_mean = 0,
                     dv_lower = c(1.14, 0.14),
                     dv_upper = c(3.14, 4.14),
                     delta_mean = 0,
                     ev_lower = c(1.14, 0.14),
                     ev_upper = c(3.14, 4.14),
                     epsilon_mean = 0,
                     lambda_x_lower = c(1, 0.5),
                     lambda_x_upper = c(1, 1.5),
                     lambda_y_lower = c(1, 0.5),
                     lambda_y_upper = c(1, 1.5),
                     nind = 5,  # number of indicators per latent variable
                     cfa_model = 
                       "# latent variable definitions (5 indicators)
                     xi =~ x1 + x2 + x3 + x4 + x5
                     eta =~ y1 + y2 + y3 + y4 + y5
                     # regression (5 indicators)
                     eta ~ xi",
                     progress = TRUE) {  # show elapsed time and status of completion
  
  # make sure input has correct format (!!!WIP!!!)
  stopifnot()
  
  # show warnings as they occur (remove #-symbol)
  #options(warn=1)
  
  # start timer if progress = TRUE
  if (isTRUE(progress)) {
    start <- proc.time()
  }
  
  # print generic information about simulation
  ntreat <- 
    32 # store number of treatments
  t0 <- as.character(Sys.time())
  sim_info <- data.frame(c(t0,
                           niter, 
                           ntreat, 
                           niter * ntreat, 
                           " "))
  colnames(sim_info) <- NULL
  rownames(sim_info) <- c("Started at: ",
                          "Iterations per treatment: ", 
                          "Number of treatments: ", 
                          "Total iterations: ", 
                          " ")
  print(sim_info)
  
  # create data frame with all treatment combinations, code each treatment combination as an integer
  template <- data.frame(niter = rep(1:niter, times = 32),
                         nobs = rep(nobs, each = 32 * niter),
                         nind = rep(nind, each = 32 * niter),
                         dv_lower = rep(dv_lower, each = niter, times = 16),
                         dv_upper = rep(dv_upper, each = niter, times = 16),
                         ev_lower = rep(ev_lower, each = 2 * niter, times = 8),
                         ev_upper = rep(ev_upper, each = 2 * niter, times = 8),
                         lambda_x_lower = rep(lambda_x_lower, each = 4 * niter, times = 4),
                         lambda_x_upper = rep(lambda_x_upper, each = 4 * niter, times = 4),
                         lambda_y_lower = rep(lambda_y_lower, each = 8 * niter, times = 2),
                         lambda_y_upper = rep(lambda_y_upper, each = 8 * niter, times = 2),
                         treatment = as.factor(rep(1:ntreat, each = niter)),
                         dv = rep(c("small", "big"), each = niter, times = 16),
                         ev = rep(c("small", "big"), each = 2 * niter, times = 8),
                         lambda_x = rep(c("small", "big"), each = 4 * niter, times = 4),
                         lambda_y = rep(c("small", "big"), each = 8 * niter, times = 2),
                         method = rep(c("LM", "SEM"), each = 16 * niter))
  template$rel_x <- NA
  template$rel_y <- NA
  
  # add additional columns used to store simulation results
  col_names <- c("b1_est", "b1_SE", "p_value", "bias_b1", "abs_bias_b1")
  template[, col_names] <- NA
  
  # add column for storing warnings and error messages
  template$warnings <- NA
  template$w <- NA
  
  # split template in half and create separate templates for LM and SEM
  LM <- template[template$method == "LM", ]
  SEM <- template[template$method == "SEM", ]
  
  # create progress bar and print information about simulation
  if (isTRUE(progress)) {
    pb <- txtProgressBar(min = 0, max = (nrow(template) / 2), style = 3)
  }
  
  # initiate iterate(), store results in data frames, iterate niter times per treatment
  for (i in 1:(nrow(template) / 2)) {
    # initiate one iteration and store its results
    iteration_res <- iterate(beta = beta, 
                             nobs = template$nobs[i],
                             xi_var = xi_var,
                             xi_mean = xi_mean,
                             zeta_var = zeta_var,
                             zeta_mean = zeta_mean,
                             dv_lower = template$dv_lower[i],
                             dv_upper = template$dv_upper[i],
                             delta_mean = delta_mean,
                             ev_lower = template$ev_lower[i],
                             ev_upper = template$ev_upper[i],
                             epsilon_mean = epsilon_mean,
                             lambda_x_lower = template$lambda_x_lower[i],
                             lambda_x_upper = template$lambda_x_upper[i],
                             lambda_y_lower = template$lambda_y_lower[i],
                             lambda_y_upper = template$lambda_y_upper[i],
                             nind = template$nind[i],
                             cfa_model = cfa_model)
    
    # extract info from results relevant for the linear model (regression) and store it
    LM$b1_est[i] <- iteration_res$LM[2, 1]
    LM$b1_SE[i] <- iteration_res$LM[2, 2]
    LM$p_value[i] <- iteration_res$LM[2, 4]
    LM$bias_b1[i] <- LM$b1_est[i] - beta[2]
    LM$abs_bias_b1[i] <- abs(LM$bias_b1[i])
    LM$w[i] <- FALSE
    LM$rel_x[i] <- iteration_res$dat[1, "rel_x"]
    LM$rel_y[i] <- iteration_res$dat[1, "rel_y"]
    
    # extract info from results relevant for the SEM and store it
    SEM$b1_est[i] <- iteration_res$SEM[1, 4]
    SEM$b1_SE[i] <- iteration_res$SEM[1, 5]
    SEM$p_value[i] <- iteration_res$SEM[1, 7]
    if (!is.na(iteration_res$SEM[1, 4])) {  # only calculate bias if b1_est is not NA
      SEM$bias_b1[i] <- SEM$b1_est[i] - beta[2]
      SEM$abs_bias_b1[i] <- abs(SEM$bias_b1[i])
    }
    SEM$warnings[i] <- iteration_res$warn["message"]
    SEM$w[i] <- iteration_res$w
    SEM$rel_x[i] <- iteration_res$dat[1, "rel_x"]
    SEM$rel_y[i] <- iteration_res$dat[1, "rel_y"]
    
    # show progress bar
    if (isTRUE(progress)) {
      setTxtProgressBar(pb, i)  # update progress bar
    }
  }
  
  # show total time elapsed
  if (isTRUE(progress)) {
    close(pb)  # close progress bar
    end <- proc.time()
    elapsed <- end - start
    print(elapsed)
  }
  
  # return parts of generated data (template data frame) and simulation results for LM and SEM
  return(list(LM = LM, SEM = SEM))
}

# Simulate and save results -----------------------------------------------
raw_data <- simulate(niter = 1000,  # number of iterations per treatment combination
                     beta = c(0, 1),
                     nobs = 500,
                     xi_var = 1,
                     xi_mean = 0,
                     zeta_var = 1,
                     zeta_mean = 0,
                     dv_lower = c(1.14, 0.14),
                     dv_upper = c(3.14, 4.14),
                     delta_mean = 0,
                     ev_lower = c(1.14, 0.14),
                     ev_upper = c(3.14, 4.14),
                     epsilon_mean = 0,
                     lambda_x_lower = c(1, 0.5),
                     lambda_x_upper = c(1, 1.5),
                     lambda_y_lower = c(1, 0.5),
                     lambda_y_upper = c(1, 1.5),
                     nind = 5,  # number of indicators per latent variable
                     cfa_model = 
                       "# latent variable definitions (5 indicators)
                     xi =~ x1 + x2 + x3 + x4 + x5
                     eta =~ y1 + y2 + y3 + y4 + y5
                     # regression (5 indicators)
                     eta ~ xi",
                     progress = TRUE)  # show elapsed time and status of completion

save(raw_data,file="raw_data_study2.Rda")

sim_res <- raw_data

# store LM and SEM as separate data frames (useful for looking at them)
LM <- sim_res$LM
SEM <- sim_res$SEM

# bind both data frames back together (used for analysis and graphical representation)
sim_res <- rbind(LM, SEM)
sim_res$method <- as.factor(sim_res$method)
sim_res$nind <- factor(sim_res$nind, labels = c("5 indicators"))

# add TRUE to significant results
sim_res$significant <- (sim_res$p_value <= 0.05)

# calculate various summarising statistics for each condition
#min <- rep(NA, nlevels(sim_res$treatment))
#max <- rep(NA, nlevels(sim_res$treatment))

mean_b1_est <- rep(NA, nlevels(sim_res$treatment))
sim_res$mean_b1_est <- mean_b1_est

mean_se <- NA
sim_res$mean_se <- NA

sd_b1_est <- NA
sim_res$sd_b1_est <- NA

sim_res$mean_bias_b1 <- NA

sim_res$sd_bias_b1 <- NA

sim_res$mean_abs_bias_b1 <- NA

sim_res$sd_abs_bias_b1 <- NA

sim_res$se_acc_treatment <- NA

sim_res$power <- NA

for (i in unique(sim_res$treatment)) {
  # create data frame for a single treatment
  temp_sr <- sim_res[sim_res$treatment == i, ]
  temp_sr <- temp_sr[((temp_sr$w == FALSE) &
                        (complete.cases(temp_sr[, c(1:24)]) == TRUE)), ]
  
  # calculate the mean and SD of the SE(B1) per treatment
  mean_se <- mean(temp_sr$b1_SE)
  
  sd_b1_est <- sd(temp_sr$b1_est)
  
  # calculate things needed for graphical representation and outcomes
  sim_res[sim_res$treatment == i, ]$mean_se <- mean_se
  
  sim_res[sim_res$treatment == i, ]$mean_b1_est <- mean(temp_sr$b1_est)
  
  sim_res[sim_res$treatment == i, ]$sd_b1_est <- sd_b1_est
  
  sim_res[sim_res$treatment == i, ]$mean_bias_b1 <- mean(temp_sr$bias_b1)
  
  sim_res[sim_res$treatment == i, ]$sd_bias_b1 <- sd(temp_sr$bias_b1)
  
  sim_res[sim_res$treatment == i, ]$mean_abs_bias_b1 <- mean(temp_sr$abs_bias_b1)
  
  sim_res[sim_res$treatment == i, ]$sd_abs_bias_b1 <- sd(temp_sr$abs_bias_b1)
  
  sim_res[sim_res$treatment == i, ]$se_acc_treatment <- (mean_se - sd_b1_est) / sd_b1_est
}

# save data frame as file
save(sim_res,file="sim_res_study2.Rda")


# PREPARE DATA FOR ANALYSIS -----------------------------------------------
# Calculate power
for (i in unique(sim_res$treatment)) {
  temp <- sim_res[((sim_res$treatment == i) &
                     (sim_res$significant == TRUE) &
                     (sim_res$w == FALSE) &
                     (complete.cases(sim_res[, -c(25, 36)]) == TRUE)), ]
  nsig <- nrow(temp)
  sim_res[(sim_res$treatment == i), ]$power <- nsig / max(sim_res$niter)
}

# save data frame as file
save(sim_res,file="sim_res_study2.Rda")