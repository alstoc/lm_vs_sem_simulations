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
                delta_var = 1, 
                delta_mean = 0,
                epsilon_var = 1,
                epsilon_mean = 0,
                lambda_x = c(1, 1, 1),  # length(lambda_x) == max(nind)
                lambda_y = c(1, 1, 1),  # length(lambda_y) == max(nind)
                nind = 3) {
  # Latent variables
  xi <- rnorm(nobs,
              xi_mean,
              sqrt(xi_var)) # true x (exogenous latent variable)
  zeta <- rnorm(nobs,
                zeta_mean,
                sqrt(zeta_var))  # residual of true y
  eta <- beta[1] + beta[2] * xi + zeta  # true y (endogenous latent variable)
  
  # Indicators and their measurement errors
  x <- matrix(NA, nrow = nobs, ncol = nind) 
  y <- matrix(NA, nrow = nobs, ncol = nind)
  delta <- matrix(NA, nrow = nobs, ncol = nind)
  epsilon <- matrix(NA, nrow = nobs, ncol = nind)
  for (i in 1:nind) {
    delta[, i] <-rnorm(nobs,
                       delta_mean,
                       sqrt(delta_var))  # measurement error for x_i
    epsilon[, i] <- rnorm(nobs,
                          epsilon_mean, 
                          sqrt(epsilon_var))  # measurement error for y_i
    x[, i] <- xi * lambda_x[i] + delta[, i]  # measured x_i
    y[, i] <- eta * lambda_y[i] + epsilon[, i]  # measured y_i
  }
  
  # Output of dgp() function with column names
  gen_dat <- data.frame(xi, delta, x, zeta, eta, epsilon, y)  # store generated data as data frame  
  delta_colnames <- paste("delta", 1:nind, sep="")
  x_colnames <- paste("x", 1:nind, sep="")
  epsilon_colnames <- paste("epsilon", 1:nind, sep="")
  y_colnames <- paste("y", 1:nind, sep="")
  names(gen_dat) <- c("xi", delta_colnames, x_colnames, "zeta", "eta", epsilon_colnames, y_colnames)
  
  # return data frame with all generated data
  return(gen_dat)
}


# Single iteration ----------------------------------------------------
iterate <- function(beta = c(0, 1), 
                    delta_var = 1,
                    epsilon_var = 1,
                    nind = 3,
                    cfa_model = 
                      "# latent variable definitions
                    xi =~ x1 + x2 + x3
                    eta =~ y1 + y2 + y3
                    # regression
                    eta ~ xi",  # WARNING: make sure to adjust cfa_model if nind != 3
                    ...) {
  # make sure lavaan is installed
  require(lavaan)
  
  # generate data using dgp()
  dat <- dgp(beta = beta, 
             nind = nind, 
             delta_var = delta_var, 
             epsilon_var = epsilon_var, ...)
  
  # create empty variable used for potential warnings
  warn <- NA
  
  # fit data with linear regression
  x_sum <- apply(dat[, (nind + 2):(2 * nind + 1)], 1, sum)  # rowsum of all indicators for xi
  y_sum <- apply(dat[, (3 * nind + 4):(4 * nind + 3)], 1, sum)  # rowsum of all indicators for eta
  lm_fit <- summary(lm(y_sum ~ x_sum))  # regression model based on rowsums of all indicators
  
  # create data frame with observed variables only
  dat_xy <- cbind(dat[, (nind + 2):(2 * nind + 1)], 
                  dat[, (3 * nind + 4):(4 * nind + 3)])
  
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
  
  # fit data with structural equation model only if data contains measurement error
  if ((delta_var != 0) & (epsilon_var != 0)) {
    # fit data with SEM and capture output
    cfa_res <- tryCatch({
      list(param = parameterEstimates(cfa(model = cfa_model, data = dat_xy)), warn = NA, w = FALSE)
    }, warning = function(w) {
      list(param = parameterEstimates(cfa(model = cfa_model, data = dat_xy)), warn = w, w = TRUE)
    }, error = function(e) {
      list(param = cfa_coefs, warn = e, w = TRUE)      
    })
  }
  
  # extract relevant information from linear regression model
  lm_coefs <- lm_fit$coefficients  # extract info about regression coefficient only
  
  # extract info about regression estimate from SEM
  if ((delta_var != 0) & (epsilon_var != 0)) {
    cfa_coefs <- cfa_res$param[(cfa_res$param$lhs == "eta")
                               & (cfa_res$param$rhs == "xi"), ] 
  }
  
  # return info about both models
  return(list(LM = lm_coefs, SEM = cfa_coefs, warn = cfa_res$warn, w = cfa_res$w, dat = dat))
}


# Simulation function -------------------------------------------------------
# Creates nind indicators of two latent variables, adds random measurement errors delta and epsilon
# to the indicators, fits linear regression model (LM) and structural equation model (SEM) to data
# and iterates niter times through all steps. Allows for comparison of LM and SEM with a varying
# number of observations, measurement error variances and number of indicators.
simulate <- function(niter = 1000,  # number of iterations per treatment combination
                     beta = c(0, 1),
                     nobs = c(50, 100, 250, 500, 1000),
                     xi_var = 1,
                     xi_mean = 0,
                     zeta_var = 1,
                     zeta_mean = 0,
                     delta_var = list(c(0, 0.157895, 1/3, 0.75, 1, 3),
                                      c(0, 0.263158, 5/9, 1.25, 5/3, 5)),
                     delta_mean = 0,
                     epsilon_var = list(c(0, 0.157895, 1/3, 0.75, 1, 3),
                                        c(0, 0.263158, 5/9, 1.25, 5/3, 5)),
                     epsilon_mean = 0,
                     lambda_x = c(1, 1, 1, 1, 1),  # length(lambda_x) == max(nind)
                     lambda_y = c(1, 1, 1, 1, 1),  # length(lambda_y) == max(nind)
                     nind = c(3, 5),  # number of indicators per latent variable
                     cfa_model = list(  # WARNING: make sure to adjust cfa_model if nind != c(3, 5)
                       model_1 =
                         "# latent variable definitions (3 indicators)
                       xi =~ x1 + x2 + x3
                       eta =~ y1 + y2 + y3
                       # regression (3 indicators)
                       eta ~ xi",
                       model_2 =
                         "# latent variable definitions (5 indicators)
                       xi =~ x1 + x2 + x3 + x4 + x5
                       eta =~ y1 + y2 + y3 + y4 + y5
                       # regression (5 indicators)
                       eta ~ xi"),
                     progress = FALSE) {  # show elapsed time and status of completion
  
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
    length(nobs) * length(nind) * length(delta_var[[1]]) * length(epsilon_var[[1]]) * 2 # store number of treatments
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
  template_1 <- expand.grid(niter = 1:niter,
                            nobs = nobs,
                            nind = nind[1],
                            delta_var = delta_var[[1]],
                            epsilon_var = epsilon_var[[1]])
  template_2 <- expand.grid(niter = 1:niter,
                            nobs = nobs,
                            nind = nind[2],
                            delta_var = delta_var[[2]],
                            epsilon_var = epsilon_var[[2]])
  template <- rbind(template_1, template_2)
  
  # add additional columns used to store simulation results
  col_names <- c("b1_est", "b1_SE", "p_value", "bias_b1", "abs_bias_b1")
  template[, col_names] <- NA
  
  # add empty column for treatment
  template$treatment <- NA
  
  # add column for storing warnings and error messages
  template$warnings <- NA
  template$w <- NA
  
  # create separate data frame for LM and SEM (used later to store simulation results)
  LM <- template
  SEM <- template
  
  # fill column with treatment information
  nlvl <- 
    length(nobs) * length(nind) * length(delta_var[[1]]) * length(epsilon_var[[1]])
  
  LM$treatment <- 
    factor(rep(1:nlvl, each = niter))
  
  SEM$treatment <- 
    factor(rep((nlvl + 1):(nlvl * 2), each = niter))
  
  # convert LM$treatment into factor with proper labels
  label <- rep(NA, nlevels(LM$treatment))
  for (i in 1:nlevels(LM$treatment)) {
    label[i] <- paste("n=", LM$nobs[i * niter], " | ",
                      "nind=", LM$nind[i * niter], " | ",
                      "delta_var=", round(LM$delta_var[i * niter], digits = 2), " | ",
                      "epsilon_var=", round(LM$epsilon_var[i * niter], digits = 2), " | ",
                      "method=LM")
  }
  LM$treatment <- factor(LM$treatment,
                         labels = label)
  
  level <- (nlevels(SEM$treatment) + 1):(nlevels(SEM$treatment) * 2)
  label <- rep(NA, nlevels(SEM$treatment))
  for (i in 1:nlevels(SEM$treatment)) {
    label[i] <- paste("n=", SEM$nobs[i * niter], " | ",
                      "nind=", SEM$nind[i * niter], " | ",
                      "delta_var=", round(SEM$delta_var[i * niter], digits = 2), " | ",
                      "epsilon_var=", round(SEM$epsilon_var[i * niter], digits = 2), " | ",
                      "method=SEM")
  }
  SEM$treatment <- factor(SEM$treatment,
                          labels = label)
  
  # create progress bar and print information about simulation
  if (isTRUE(progress)) {
    pb <- txtProgressBar(min = 0, max = nrow(template), style = 3)
  }
  
  # initiate iterate(), store results in data frames, iterate niter times per treatment
  for (i in 1:nrow(template)) {
    # create variable for choosing which parameters to use as iterate() argument
    # (delta_var, epsilon_var and cfa_model depend on number of indicators)
    if (template$nind[i] == nind[1]) {
      j <- 1
    } else if (template$nind[i] == nind[2]) {
      j <- 2
    }
    
    
    # initiate one iteration and store its results
    iteration_res <- iterate(beta = beta, 
                             nobs = template$nobs[i],
                             xi_var = xi_var,
                             xi_mean = xi_mean,
                             zeta_var = zeta_var,
                             zeta_mean = zeta_mean,
                             delta_var = template$delta_var[i],
                             delta_mean = delta_mean,
                             epsilon_var = template$epsilon_var[i],
                             epsilon_mean = epsilon_mean,
                             lambda_x = lambda_x,
                             lambda_y = lambda_y,
                             nind = template$nind[i],
                             cfa_model = cfa_model[[j]])
    
    # extract info from results relevant for the linear model (regression) and store it
    LM$b1_est[i] <- iteration_res$LM[2, 1]
    LM$b1_SE[i] <- iteration_res$LM[2, 2]
    LM$p_value[i] <- iteration_res$LM[2, 4]
    LM$bias_b1[i] <- LM$b1_est[i] - beta[2]
    LM$abs_bias_b1[i] <- abs(LM$bias_b1[i])
    LM$w[i] <- FALSE
    
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

# Simulate and save results for v1.3 -----------------------------------------------
raw_data <- simulate(niter = 1000,  # number of iterations per treatment combination
                     beta = c(0, 1),
                     nobs = c(50, 100, 250, 500, 1000),
                     xi_var = 1,
                     xi_mean = 0,
                     zeta_var = 4,
                     zeta_mean = 0,
                     delta_var = list(c(0, 0.157895, 1/3, 0.75, 1.285714, 3),
                                      c(0, 0.263158, 5/9, 1.25, 2.142857, 5)),
                     delta_mean = 0,
                     epsilon_var = list(c(0, 0.157895, 1/3, 0.75, 1.285714, 3),
                                        c(0, 0.263158, 5/9, 1.25, 2.142857, 5)),
                     epsilon_mean = 0,
                     lambda_x = c(1, 1, 1, 1, 1),  # length(lambda_x) == max(nind)
                     lambda_y = c(1, 1, 1, 1, 1),  # length(lambda_y) == max(nind)
                     nind = c(3, 5),  # number of indicators per latent variable
                     cfa_model = list(  # WARNING: make sure to adjust cfa_model if nind != c(3, 5)
                       model_1 =
                         "# latent variable definitions (3 indicators)
                       xi =~ x1 + x2 + x3
                       eta =~ y1 + y2 + y3
                       # regression (3 indicators)
                       eta ~ xi",
                       model_2 =
                         "# latent variable definitions (5 indicators)
                       xi =~ x1 + x2 + x3 + x4 + x5
                       eta =~ y1 + y2 + y3 + y4 + y5
                       # regression (5 indicators)
                       eta ~ xi"),
                     progress = TRUE)  # show elapsed time and status of completion

save(raw_data,file="raw_data_study1_3.Rda")

sim_res <- raw_data

# store LM and SEM as separate data frames (useful for looking at them)
LM <- sim_res$LM
SEM <- sim_res$SEM

# add method column
LM$method <- "LM"
SEM$method <- "SEM"

# add column for reliability
rel_x <- rep(c(1, 0.95, 0.9, 0.8, 0.7, 0.5), each = 5000, times = 12)
rel_y <- rep(c(1, 0.95, 0.9, 0.8, 0.7, 0.5), each = 30000, times = 2)
LM <- cbind(LM[, c(1:3)], rel_x, rel_y, LM[, -c(1:3)])
SEM <- cbind(SEM[, c(1:3)], rel_x, rel_y, SEM[, -c(1:3)])

# bind both dataframes back together (used for analysis and graphical representation)
sim_res <- rbind(LM, SEM)
sim_res$method <- as.factor(sim_res$method)
sim_res$nind <- factor(sim_res$nind, labels = c("3 indicators", "5 indicators"))

# add TRUE to significant results
sim_res$significant <- (sim_res$p_value <= 0.05)

# calculate SD and mean of B1 est. and SE(B1) in each condition for ggplot
min <- rep(NA, nlevels(sim_res$treatment))
max <- rep(NA, nlevels(sim_res$treatment))
mean_se <- rep(NA, nlevels(sim_res$treatment))
sd_b1_est <- rep(NA, nlevels(sim_res$treatment))
mean_b1_est <- rep(NA, nlevels(sim_res$treatment))
se_acc <- rep(NA, nlevels(sim_res$treatment))

for (i in 1:(nlevels(sim_res$treatment))) {
  # define lower and upper limits of each condition (which each have 1000 iterations)
  min[i] <- (1 + (1000 * ((i - 1))))
  max[i] <- (1000 * i)
  
  temp_sr <- sim_res[c(min[i]:max[i]), ]
  temp_sr <- temp_sr[((temp_sr$w == FALSE) &
                        (complete.cases(temp_sr[, c(1:12)]) == TRUE)), ]
  
  # calculate the mean of the SE(B1) per condition
  mean_se[i] <- mean(temp_sr$b1_SE)
  sim_res$mean_se[min[i]:max[i]] <- mean_se[i]
  
  # calculate things needed for graphical representation and outcomes
  sd_b1_est[i] <- sd(temp_sr$b1_est)
  sim_res$sd_b1_est[min[i]:max[i]] <- sd_b1_est[i]
  
  mean_b1_est[i] <-mean(temp_sr$b1_est)
  sim_res$mean_b1_est[min[i]:max[i]] <- mean_b1_est[i]
  
  se_acc[i] <- (mean_se[i] - sd_b1_est[i]) / sd_b1_est[i]
  sim_res$se_acc[min[i]:max[i]] <- se_acc[i]
}

# PREPARE DATA FOR ANALYSIS
# Create factor for each level of reliability for x and y (used in simple effects analysis)
# Do this by splitting sim_res by method and number of indicators
sim_res_lm_3 <- sim_res[(sim_res$nind == "3 indicators") &
                          (sim_res$method == "LM"), ]
sim_res_lm_3 <- sim_res_lm_3[order(sim_res_lm_3$rel_x), ]  # sort by rel_x
sim_res_lm_3$simple_method_relx <- rep(1:6, each = 30000)  # create factor for simple effects
sim_res_lm_3$simple_method_relx <- factor(sim_res_lm_3$simple_method_relx,
                                          levels = c(1:6),
                                          labels = c(
                                            "0.5_relx_3ind_LM",
                                            "0.7_relx_3ind_LM",
                                            "0.8_relx_3ind_LM",
                                            "0.9_relx_3ind_LM",
                                            "0.95_relx_3ind_LM",
                                            "1.0_relx_3ind_LM"
                                          ))
sim_res_lm_3 <- sim_res_lm_3[order(sim_res_lm_3$rel_y), ]  # sort by rel_y
sim_res_lm_3$simple_method_rely <- rep(1:6, each = 30000)  # create factor for simple effects
sim_res_lm_3$simple_method_rely <- factor(sim_res_lm_3$simple_method_rely,
                                          levels = c(1:6),
                                          labels = c(
                                            "0.5_rely_3ind_LM",
                                            "0.7_rely_3ind_LM",
                                            "0.8_rely_3ind_LM",
                                            "0.9_rely_3ind_LM",
                                            "0.95_rely_3ind_LM",
                                            "1.0_rely_3ind_LM"
                                          ))


sim_res_sem_3 <- sim_res[(sim_res$nind == "3 indicators") &
                           (sim_res$method == "SEM"), ]
sim_res_sem_3 <- sim_res_sem_3[order(sim_res_sem_3$rel_x), ]
sim_res_sem_3$simple_method_relx <- rep(7:12, each = 30000)
sim_res_sem_3$simple_method_relx <- factor(sim_res_sem_3$simple_method_relx,
                                           levels = c(7:12),
                                           labels = c(
                                             "0.5_relx_3ind_SEM",
                                             "0.7_relx_3ind_SEM",
                                             "0.8_relx_3ind_SEM",
                                             "0.9_relx_3ind_SEM",
                                             "0.95_relx_3ind_SEM",
                                             "1.0_relx_3ind_SEM"
                                           ))
sim_res_sem_3 <- sim_res_sem_3[order(sim_res_sem_3$rel_y), ]
sim_res_sem_3$simple_method_rely <- rep(7:12, each = 30000) 
sim_res_sem_3$simple_method_rely <- factor(sim_res_sem_3$simple_method_rely,
                                           levels = c(7:12),
                                           labels = c(
                                             "0.5_rely_3ind_SEM",
                                             "0.7_rely_3ind_SEM",
                                             "0.8_rely_3ind_SEM",
                                             "0.9_rely_3ind_SEM",
                                             "0.95_rely_3ind_SEM",
                                             "1.0_rely_3ind_SEM"
                                           ))


sim_res_lm_5 <- sim_res[(sim_res$nind == "5 indicators") &
                          (sim_res$method == "LM"), ]
sim_res_lm_5 <- sim_res_lm_5[order(sim_res_lm_5$rel_x), ]
sim_res_lm_5$simple_method_relx <- rep(13:18, each = 30000)
sim_res_lm_5$simple_method_relx <- factor(sim_res_lm_5$simple_method_relx,
                                          levels = c(13:18),
                                          labels = c(
                                            "0.5_relx_5ind_LM",
                                            "0.7_relx_5ind_LM",
                                            "0.8_relx_5ind_LM",
                                            "0.9_relx_5ind_LM",
                                            "0.95_relx_5ind_LM",
                                            "1.0_relx_5ind_LM"
                                          ))
sim_res_lm_5 <- sim_res_lm_5[order(sim_res_lm_5$rel_y), ]
sim_res_lm_5$simple_method_rely <- rep(13:18, each = 30000)
sim_res_lm_5$simple_method_rely <- factor(sim_res_lm_5$simple_method_rely,
                                          levels = c(13:18),
                                          labels = c(
                                            "0.5_rely_5ind_LM",
                                            "0.7_rely_5ind_LM",
                                            "0.8_rely_5ind_LM",
                                            "0.9_rely_5ind_LM",
                                            "0.95_rely_5ind_LM",
                                            "1.0_rely_5ind_LM"
                                          ))


sim_res_sem_5 <- sim_res[(sim_res$nind == "5 indicators") &
                           (sim_res$method == "SEM"), ]
sim_res_sem_5 <- sim_res_sem_5[order(sim_res_sem_5$rel_x), ]
sim_res_sem_5$simple_method_relx <- rep(19:24, each = 30000)
sim_res_sem_5$simple_method_relx <- factor(sim_res_sem_5$simple_method_relx,
                                           levels = c(19:24),
                                           labels = c(
                                             "0.5_relx_5ind_SEM",
                                             "0.7_relx_5ind_SEM",
                                             "0.8_relx_5ind_SEM",
                                             "0.9_relx_5ind_SEM",
                                             "0.95_relx_5ind_SEM",
                                             "1.0_relx_5ind_SEM"
                                           ))
sim_res_sem_5 <- sim_res_sem_5[order(sim_res_sem_5$rel_y), ]
sim_res_sem_5$simple_method_rely <- rep(19:24, each = 30000) 
sim_res_sem_5$simple_method_rely <- factor(sim_res_sem_5$simple_method_rely,
                                           levels = c(19:24),
                                           labels = c(
                                             "0.5_rely_5ind_SEM",
                                             "0.7_rely_5ind_SEM",
                                             "0.8_rely_5ind_SEM",
                                             "0.9_rely_5ind_SEM",
                                             "0.95_rely_5ind_SEM",
                                             "1.0_rely_5ind_SEM"
                                           ))


# save simple effects factors in sim_res and turn them into proper factors
sim_res <- rbind(sim_res_lm_3,
                 sim_res_sem_3,
                 sim_res_lm_5,
                 sim_res_sem_5)

# save data frame as file
save(sim_res,file="sim_res_study1_3.Rda")