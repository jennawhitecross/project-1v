################################################################################
# TRUE VUS / METRIC FUNCTIONS
################################################################################

calc.bias <- function(vus.est, vus.true, n.sim) {
  bias <- mean(vus.est - vus.true, na.rm = TRUE)  
  se <- sd(vus.est - vus.true, na.rm = TRUE) / sqrt(sum(!is.na(vus.est)))
  return(list(bias = bias, se = se))
}

calc.rmse <- function(vus.est, vus.true, n.sim) {
  mse <- mean((vus.est - vus.true)^2, na.rm = TRUE)
  rmse <- sqrt(mse)
  se <- sd((vus.est - vus.true)^2, na.rm = TRUE) / sqrt(sum(!is.na(vus.est)))
  return(list(rmse = rmse, se = se))
}

calc.cdf <- function(x, distbn, theta){
  if (identical(distbn, rnorm)){
    cdf <- pnorm(x, theta[1], theta[2])
  } else if (identical(distbn, rgamma)) {
    cdf <- pgamma(x, theta[1], theta[2])
  } else if (identical(distbn, rlnorm)) {
    cdf <- plnorm(x, theta[1], theta[2])
  }
  return(cdf)
}

calc.inverse.cdf <- function(p, distbn, theta){
  if (identical(distbn, rnorm)){
    inverse.cdf <- qnorm(p, theta[1], theta[2])
  } else if (identical(distbn, rgamma)) {
    inverse.cdf <- qgamma(p, theta[1], theta[2])
  } else if (identical(distbn, rlnorm)) {
    inverse.cdf <- qlnorm(p, theta[1], theta[2])
  }
  return(inverse.cdf)
}

calc.true.vus <- function(theta.X, theta.Y, theta.Z, distbn.X, distbn.Y, distbn.Z) {
  
  roc <- function(p1, p3) {
    lower <- calc.inverse.cdf(p1, distbn.X, theta.X)
    upper <- calc.inverse.cdf(1 - p3, distbn.Z, theta.Z)
    ifelse(lower <= upper,
           calc.cdf(upper, distbn.Y, theta.Y) - calc.cdf(lower, distbn.Y, theta.Y),
           0)
  }
  integrand <- function(p3) integrate(function(p1) roc(p1, p3), 0.001, 0.99)$value
  vus <- integrate(Vectorize(integrand), 0.001, 0.99)$value
  
  return(vus)
}

calc.sample <- function(n, theta, distbn){
  if (identical(distbn, rnorm)){
    sample <- rnorm(n, theta[1], theta[2])
  } else if (identical(distbn, rgamma)) {
    sample <- rgamma(n, theta[1], theta[2])
  } else if (identical(distbn, rlnorm)) {
    sample <- rlnorm(n, theta[1], theta[2])
  }
  return(sample)
}

calc.coverage <- function(ci.list, vus.true) {
  if(length(ci.list) == 0) return(NA)
  
  count <- 0
  n <- length(ci.list)
  
  for (i in 1:n){
    ci <- ci.list[[i]]
    if (!is.null(ci) && length(ci) == 2 && !any(is.na(ci))) {
      if (vus.true >= ci[1] && vus.true <= ci[2]) {
        count <- count + 1
      }
    }
  }
  coverage <- count / n
  return(coverage)
}

calc.avg.length <- function(ci.list) {
  if(length(ci.list) == 0) return(NA)
  
  n <- length(ci.list)
  lengths <- numeric(0)
  
  for (i in 1:n){
    ci <- ci.list[[i]]
    if (!is.null(ci) && length(ci) == 2 && !any(is.na(ci))) {
      lengths <- c(lengths, ci[2] - ci[1])
    }
  }
  
  if (length(lengths) == 0) return(NA)
  avg.length <- mean(lengths)
  return(avg.length)
}

calc.cdtional.length <- function(ci.list, vus.true) {
  if(length(ci.list) == 0) return(NA)
  
  n <- length(ci.list)
  covering.lengths <- numeric(0)
  
  for (i in 1:n) {
    ci <- ci.list[[i]]
    if (!is.null(ci) && length(ci) == 2 && !any(is.na(ci))) {
      if (vus.true >= ci[1] && vus.true <= ci[2]) {
        covering.lengths <- c(covering.lengths, ci[2] - ci[1])
      }
    }
  }
  
  if (length(covering.lengths) == 0) return(NA)
  return(mean(covering.lengths))
}

################################################################################
# SIMULATION FUNCTIONS FOR SCENARIOS 1 TO 6
################################################################################

library(boot)
library(trinROC)
library(nleqslv)
library(numDeriv)
library(survival)
library(MASS)

repetition <- function(X, Y, Z, tune.X, tune.Y, tune.Z, 
                       u.var.inc, bc.inc, n.boot = 1000, alpha = 0.05) {
  
  # point estimates
  ep.vus <- emp.vus(X, Y, Z, dat = NULL, old.version = TRUE)
  k1.vus <- k1.method(X,Y,Z)
  mle.out <- trinVUS.test(X, Y, Z)
  mle.vus <- mle.out$estimate
  if (bc.inc == TRUE) {
    bc.XYZ <- boxcoxROC(X, Y, Z, verbose = FALSE)
    lambda.opt <- bc.XYZ$lambda
    bc.out <- trinVUS.test(bc.XYZ$xbc, bc.XYZ$ybc, bc.XYZ$zbc)
    bc.vus <- bc.out$estimate
  }else{
    bc.vus <- NA
  }
  leh.out <- leh.method(X, Y, Z)
  leh.vus <- leh.out$vus
  mfpt.out <- mfpt.method(X, Y, Z, iterations=1000, burn.in=100, thin=5, 
                                      tuneX=tune.X, tuneY=tune.Y, tuneZ=tune.Z, c=1)
  mfpt.vus <- mfpt.out$vus
  mfpt.acc1 <- mfpt.out$accept1
  mfpt.acc2 <- mfpt.out$accept2
  mfpt.acc3 <- mfpt.out$accept3
  Sigma <- diag(c(0.01, 0.01, 0.01, 0.005, 0.005, 0.005))
  mh.out <- mh.method(X, Y, Z, iterations=10000, burn.in=1000, thin=10, Sigma=Sigma)
  mh.vus <- mh.out$vus
  mh.acc <- mh.out$acc.prob
  
  # asymtpotic
  if (u.var.inc == TRUE) {
    u.ci <- u.var(X,Y,Z,ep.vus)$ci
  }else{
    u.ci <- c(NA,NA)
  }
  mle.ci <- mle.out$conf.int
  if (bc.inc == TRUE) {
    bc.ci <- bc.out$conf.int
  }else{
    bc.ci <- c(NA,NA)
  }
  leh.ci <- leh.out$ci
  
  # z-transformation
  # mle.theta <- z.transformation(mle.vus)
  # var.theta <- variance(mle.vus, mle.out$Sigma[1,1])
  # mle.z.ci <- z.inverse(c(mle.theta - 1.96*sqrt(var.theta), mle.theta + 1.96*sqrt(var.theta)))
  
  # bootstrap
  data.boot <- data.frame(
    marker = c(X, Y, Z),
    class = c(rep(1, length(X)), rep(2, length(Y)), rep(3, length(Z)))
  )
  ep.out.boot <- boot(data.boot, ep.boot, R = n.boot)
  ep.ci.boot <- boot.ci(ep.out.boot, type = c("perc", "bca"), conf = 0.95)
  k1.out.boot <- boot(data.boot, k1.boot, R = n.boot)
  k1.ci.boot <- boot.ci(k1.out.boot, type = c("perc", "bca"), conf = 0.95)
  mle.out.boot <- boot(data.boot, mle.boot, R = n.boot)
  mle.ci.boot <- boot.ci(mle.out.boot, type = c("perc", "bca", "stud"), conf = 0.95)
  if (bc.inc == TRUE) {
    bc.out.boot <- boot(data.boot, bc.boot, R = n.boot)
    bc.ci.boot <- boot.ci(bc.out.boot, type = c("perc", "bca", "stud"), conf = 0.95)
  }else{
    bc.ci.boot <- list(percent = c(NA, NA, NA, NA, NA),
                       bca = c(NA, NA, NA, NA, NA),
                       student = c(NA, NA, NA, NA, NA))
  }
  leh.out.boot <- boot(data.boot, leh.boot, R = n.boot)
  leh.ci.boot <- boot.ci(leh.out.boot, type = c("perc", "bca", "stud"), conf = 0.95)
 
  # bayesian
  mfpt.ci <- bayes.cred(mfpt.out$vus.chain)
  mh.ci <- bayes.cred(mh.out$vus.chain)
  
  # results
  vus.est <- list(ep = ep.vus, k1 = k1.vus, mle = mle.vus, bc = bc.vus, leh = leh.vus, mfpt = mfpt.vus, mh = mh.vus)
  # vus.est <- list(ep = ep.vus, k1 = k1.vus, mle = mle.vus, mle.z = z.inverse(mle.theta), bc = bc.vus, leh = leh.vus)
  asymptotic.est <- list(u = u.ci, mle = mle.ci, bc = bc.ci, leh = leh.ci)
  # asymptotic.est <- list(u = u.ci, mle = mle.ci, mle.z = mle.z.ci, bc = bc.ci, leh = leh.ci)
  percent.est <- list(ep = ep.ci.boot$percent[4:5], k1 = k1.ci.boot$percent[4:5], mle = mle.ci.boot$percent[4:5], bc = bc.ci.boot$percent[4:5], leh = leh.ci.boot$percent[4:5])
  bca.est <- list(ep = ep.ci.boot$bca[4:5], k1 = k1.ci.boot$bca[4:5], mle = mle.ci.boot$bca[4:5], bc = bc.ci.boot$bca[4:5], leh = leh.ci.boot$bca[4:5])
  student.est <- list(mle = mle.ci.boot$student[4:5], bc = bc.ci.boot$student[4:5], leh = leh.ci.boot$student[4:5])
  bayes.est <- list(mfpt.equi = mfpt.ci$equi, mfpt.hpd = mfpt.ci$hpd, mh.equi = mh.ci$equi, mh.hpd = mh.ci$hpd)
  acc.rate <- list(mfpt.mu1 = mfpt.acc1[1], mfpt.sigma1 = mfpt.acc1[2],
                   mfpt.mu2 = mfpt.acc2[1], mfpt.sigma2 = mfpt.acc2[2],
                   mfpt.mu3 = mfpt.acc3[1], mfpt.sigma3 = mfpt.acc3[2], 
                   mh = mh.acc)
  
  return(list(
    vus = vus.est,
    asymptotic = asymptotic.est,
    percent = percent.est,
    bca = bca.est,
    student = student.est,
    bayes = bayes.est,
    acc = acc.rate
  ))
}

simulation <- function(n, theta.X, theta.Y, theta.Z, distbn.X, distbn.Y, distbn.Z, tune.X, tune.Y, tune.Z, u.var.inc, bc.inc, n.sim = 50) {
  
  vus.true <- calc.true.vus(theta.X, theta.Y, theta.Z, distbn.X, distbn.Y, distbn.Z)
  
  if (bc.inc == TRUE){
    method.names <- c("EP", "K1", "MLE", "BC", "LEH", "MFPT", "MH")
    ci.method.names <- c("U", "MLE", "BC", "LEH",
                         "EP.Percent", "EP.BCa",
                         "K1.Percent", "K1.BCa",
                         "MLE.Percent", "MLE.BCa", "MLE.Student",
                         "BC.Percent", "BC.BCa", "BC.Student",
                         "LEH.Percent", "LEH.BCa", "LEH.Student",
                         "MFPT.Equi", "MFPT.HPD", "MH.Equi", "MH.HPD")
    acc.rate.names <- c("MFPT.mu1", "MFPT.mu2", "MFPT.mu3",
                        "MFPT.sigma1", "MFPT.sigma2", "MFPT.sigma3",
                        "MH.vus")
  }else{
    method.names <- c("EP", "K1", "MLE", "LEH", "MFPT", "MH")
    # method.names <- c("EP", "K1", "MLE", "MLE.Z", "LEH")
    ci.method.names <- c("U", "MLE", "LEH",
                         "EP.Percent", "EP.BCa",
                         "K1.Percent", "K1.BCa",
                         "MLE.Percent", "MLE.BCa", "MLE.Student",
                         "LEH.Percent", "LEH.BCa", "LEH.Student",
                         "MFPT.Equi", "MFPT.HPD", "MH.Equi", "MH.HPD")
    # ci.method.names <- c("U", "MLE", "MLE.Z", "LEH",
    #                      "EP.Percent", "EP.BCa",
    #                      "K1.Percent", "K1.BCa",
    #                      "MLE.Percent", "MLE.BCa", "MLE.Student",
    #                      "LEH.Percent", "LEH.BCa", "LEH.Student")
    acc.rate.names <- c("MFPT.mu1", "MFPT.mu2", "MFPT.mu3",
                        "MFPT.sigma1", "MFPT.sigma2", "MFPT.sigma3",
                        "MH.vus")
  }
  
  results <- list()
  for (i in 1:n.sim){
    X <- calc.sample(n[1], theta.X, distbn.X)
    Y <- calc.sample(n[2], theta.Y, distbn.Y)
    Z <- calc.sample(n[3], theta.Z, distbn.Z)
    results[[i]] <- repetition(X, Y, Z, tune.X, tune.Y, tune.Z, u.var.inc, bc.inc)
    print(i)
  }
  
  vus.estimates <- matrix(NA, n.sim, length(method.names))
  colnames(vus.estimates) <- method.names
  
  acc.rates <- matrix(NA, n.sim, length(acc.rate.names))
  colnames(acc.rates) <- acc.rate.names
  
  if (bc.inc == TRUE){
    
    ci.estimates <- list(
      U = list(),
      MLE = list(),
      BC = list(),
      LEH = list(),
      EP.Percent = list(),
      EP.BCa = list(),
      K1.Percent = list(),
      K1.BCa = list(),
      MLE.Percent = list(),
      MLE.BCa = list(),
      MLE.Student = list(),
      BC.Percent = list(),
      BC.BCa = list(),
      BC.Student = list(),
      LEH.Percent = list(),
      LEH.BCa = list(),
      LEH.Student = list(),
      MFPT.Equi = list(),
      MFPT.HPD = list(),
      MH.Equi = list(),
      MH.HPD = list()
    )
    
    for(i in 1:n.sim) {
      vus.estimates[i, "EP"] <- results[[i]]$vus$ep
      vus.estimates[i, "K1"] <- results[[i]]$vus$k1
      vus.estimates[i, "MLE"] <- results[[i]]$vus$mle
      vus.estimates[i, "BC"] <- results[[i]]$vus$bc
      vus.estimates[i, "LEH"] <- results[[i]]$vus$leh
      vus.estimates[i, "MFPT"] <- results[[i]]$vus$mfpt
      vus.estimates[i, "MH"] <- results[[i]]$vus$mh
      ci.estimates$U[[i]] <- results[[i]]$asymptotic$u
      ci.estimates$MLE[[i]] <- results[[i]]$asymptotic$mle
      ci.estimates$BC[[i]] <- results[[i]]$asymptotic$bc
      ci.estimates$LEH[[i]] <- results[[i]]$asymptotic$leh
      ci.estimates$EP.Percent[[i]] <- results[[i]]$percent$ep
      ci.estimates$EP.BCa[[i]] <- results[[i]]$bca$ep
      ci.estimates$K1.Percent[[i]] <- results[[i]]$percent$k1
      ci.estimates$K1.BCa[[i]] <- results[[i]]$bca$k1
      ci.estimates$MLE.Percent[[i]] <- results[[i]]$percent$mle
      ci.estimates$MLE.BCa[[i]] <- results[[i]]$bca$mle
      ci.estimates$MLE.Student[[i]] <- results[[i]]$student$mle
      ci.estimates$BC.Percent[[i]] <- results[[i]]$percent$bc
      ci.estimates$BC.BCa[[i]] <- results[[i]]$bca$bc
      ci.estimates$BC.Student[[i]] <- results[[i]]$student$bc
      ci.estimates$LEH.Percent[[i]] <- results[[i]]$percent$leh
      ci.estimates$LEH.BCa[[i]] <- results[[i]]$bca$leh
      ci.estimates$LEH.Student[[i]] <- results[[i]]$student$leh
      ci.estimates$MFPT.Equi[[i]] <- results[[i]]$bayes$mfpt.equi
      ci.estimates$MFPT.HPD[[i]] <- results[[i]]$bayes$mfpt.hpd
      ci.estimates$MH.Equi[[i]] <- results[[i]]$bayes$mh.equi
      ci.estimates$MH.HPD[[i]] <- results[[i]]$bayes$mh.hpd
      acc.rates[i, "MFPT.mu1"] <- results[[i]]$acc$mfpt.mu1
      acc.rates[i, "MFPT.mu2"] <- results[[i]]$acc$mfpt.mu2
      acc.rates[i, "MFPT.mu3"] <- results[[i]]$acc$mfpt.mu3
      acc.rates[i, "MFPT.sigma1"] <- results[[i]]$acc$mfpt.sigma1
      acc.rates[i, "MFPT.sigma2"] <- results[[i]]$acc$mfpt.sigma2
      acc.rates[i, "MFPT.sigma3"] <- results[[i]]$acc$mfpt.sigma3
      acc.rates[i, "MH.vus"] <- results[[i]]$acc$mh
    }
  }else{
    
    ci.estimates <- list(
      U = list(),
      MLE = list(),
      # MLE.Z = list(),
      LEH = list(),
      EP.Percent = list(),
      EP.BCa = list(),
      K1.Percent = list(),
      K1.BCa = list(),
      MLE.Percent = list(),
      MLE.BCa = list(),
      MLE.Student = list(),
      LEH.Percent = list(),
      LEH.BCa = list(),
      LEH.Student = list(),
      MFPT.Equi = list(),
      MFPT.HPD = list(),
      MH.Equi = list(),
      MH.HPD = list()
    )
    
    for(i in 1:n.sim) {
      vus.estimates[i, "EP"] <- results[[i]]$vus$ep
      vus.estimates[i, "K1"] <- results[[i]]$vus$k1
      vus.estimates[i, "MLE"] <- results[[i]]$vus$mle
      # vus.estimates[i, "MLE.Z"] <- results[[i]]$vus$mle.z
      vus.estimates[i, "LEH"] <- results[[i]]$vus$leh
      ci.estimates$U[[i]] <- results[[i]]$asymptotic$u
      ci.estimates$MLE[[i]] <- results[[i]]$asymptotic$mle
      # ci.estimates$MLE.Z[[i]] <- results[[i]]$asymptotic$mle.z
      ci.estimates$LEH[[i]] <- results[[i]]$asymptotic$leh
      vus.estimates[i, "MFPT"] <- results[[i]]$vus$mfpt
      vus.estimates[i, "MH"] <- results[[i]]$vus$mh
      ci.estimates$EP.Percent[[i]] <- results[[i]]$percent$ep
      ci.estimates$EP.BCa[[i]] <- results[[i]]$bca$ep
      ci.estimates$K1.Percent[[i]] <- results[[i]]$percent$k1
      ci.estimates$K1.BCa[[i]] <- results[[i]]$bca$k1
      ci.estimates$MLE.Percent[[i]] <- results[[i]]$percent$mle
      ci.estimates$MLE.BCa[[i]] <- results[[i]]$bca$mle
      ci.estimates$MLE.Student[[i]] <- results[[i]]$student$mle
      ci.estimates$LEH.Percent[[i]] <- results[[i]]$percent$leh
      ci.estimates$LEH.BCa[[i]] <- results[[i]]$bca$leh
      ci.estimates$LEH.Student[[i]] <- results[[i]]$student$leh
      ci.estimates$MFPT.Equi[[i]] <- results[[i]]$bayes$mfpt.equi
      ci.estimates$MFPT.HPD[[i]] <- results[[i]]$bayes$mfpt.hpd
      ci.estimates$MH.Equi[[i]] <- results[[i]]$bayes$mh.equi
      ci.estimates$MH.HPD[[i]] <- results[[i]]$bayes$mh.hpd
      acc.rates[i, "MFPT.mu1"] <- results[[i]]$acc$mfpt.mu1
      acc.rates[i, "MFPT.mu2"] <- results[[i]]$acc$mfpt.mu2
      acc.rates[i, "MFPT.mu3"] <- results[[i]]$acc$mfpt.mu3
      acc.rates[i, "MFPT.sigma1"] <- results[[i]]$acc$mfpt.sigma1
      acc.rates[i, "MFPT.sigma2"] <- results[[i]]$acc$mfpt.sigma2
      acc.rates[i, "MFPT.sigma3"] <- results[[i]]$acc$mfpt.sigma3
      acc.rates[i, "MH.vus"] <- results[[i]]$acc$mh
    }
  }
  
  bias.result <- matrix(NA, 2, length(method.names))
  rmse.result <- matrix(NA, 2, length(method.names))
  colnames(bias.result) <- method.names
  colnames(rmse.result) <- method.names
  rownames(bias.result) <- c("Bias", "SE")
  rownames(rmse.result) <- c("RMSE", "SE")
  
  acc.rate.result <- numeric(length(acc.rate.names))
  names(acc.rate.result) <- acc.rate.names
  
  for (i in 1:length(method.names)) {
    vus.est <- vus.estimates[, i]
    
    bias.result[1, i] <- calc.bias(vus.est, vus.true, n.sim)$bias * 1000
    bias.result[2, i] <- calc.bias(vus.est, vus.true, n.sim)$se * 1000
    rmse.result[1, i] <- calc.rmse(vus.est, vus.true, n.sim)$rmse * 1000
    rmse.result[2, i] <- calc.rmse(vus.est, vus.true, n.sim)$se * 1000
  }
  
  for (i in 1:length(acc.rate.names)){
    acc.rate.result[i] <- mean(acc.rates[, i])
  }
  
  coverage.result <- matrix(NA, 1, length(ci.method.names))
  avglength.result <- matrix(NA, 1, length(ci.method.names))
  condlength.result <- matrix(NA, 1, length(ci.method.names))
  
  colnames(coverage.result) <- ci.method.names
  colnames(avglength.result) <- ci.method.names
  colnames(condlength.result) <- ci.method.names
  
  for (i in 1:length(ci.method.names)) {
    ci.list <- ci.estimates[[ci.method.names[i]]]
    coverage.result[1, i] <- calc.coverage(ci.list, vus.true)
    avglength.result[1, i] <- calc.avg.length(ci.list)
    condlength.result[1, i] <- calc.cdtional.length(ci.list, vus.true)
  }
  
  return(list(
    vus.true = vus.true,
    bias = bias.result,
    rmse = rmse.result,
    coverage = coverage.result,
    avg.length = avglength.result,
    cond.length = condlength.result,
    acc.rate = acc.rate.result
  ))
}

set.seed(1)

results.case1.1 <- simulation(
  n = c(20, 20, 20),
  theta.X = c(0, 1),
  theta.Y = c(5, 1),
  theta.Z = c(10, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.5,6),
  tune.Y = c(0.5,6),
  tune.Z = c(0.5,6),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case1.2 <- simulation(
  n = c(20, 20, 30),
  theta.X = c(0, 1),
  theta.Y = c(5, 1),
  theta.Z = c(10, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.5,6),
  tune.Y = c(0.5,6),
  tune.Z = c(0.5,6),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case1.3 <- simulation(
  n = c(20, 30, 50),
  theta.X = c(0, 1),
  theta.Y = c(5, 1),
  theta.Z = c(10, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = (c(0.2,12) + c(0.5,6)) / 2,
  tune.Y = (c(0.2,12) + c(0.5,6)) / 2,
  tune.Z = (c(0.2,12) + c(0.5,6)) / 2,
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case1.4 <- simulation(
  n = c(30, 50, 50),
  theta.X = c(0, 1),
  theta.Y = c(5, 1),
  theta.Z = c(10, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.2,12),
  tune.Y = c(0.2,12),
  tune.Z = c(0.2,12),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case1.5 <- simulation(
  n = c(50, 50, 50),
  theta.X = c(0, 1),
  theta.Y = c(5, 1),
  theta.Z = c(10, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.2,12),
  tune.Y = c(0.2,12),
  tune.Z = c(0.2,12),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)


results.case2.1 <- simulation(
  n = c(20, 20, 20),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.3,3),
  tune.Y = c(0.3,3),
  tune.Z = c(0.3,3),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case2.2 <- simulation(
  n = c(20, 20, 30),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.3,3),
  tune.Y = c(0.3,3),
  tune.Z = c(0.3,3),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case2.3 <- simulation(
  n = c(20, 30, 50),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = (c(0.3,3) + c(0.15,6)) / 2,
  tune.Y = (c(0.3,3) + c(0.15,6)) / 2,
  tune.Z = (c(0.3,3) + c(0.15,6)) / 2,
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case2.4 <- simulation(
  n = c(30, 50, 50),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.15,6),
  tune.Y = c(0.15,6),
  tune.Z = c(0.15,6),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case2.5 <- simulation(
  n = c(50, 50, 50),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.15,6),
  tune.Y = c(0.15,6),
  tune.Z = c(0.15,6),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case3.1 <- simulation(
  n = c(20, 20, 20),
  theta.X = c(0, 1),
  theta.Y = c(0, 1),
  theta.Z = c(0, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,2),
  tune.Z = c(0.2,2),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case3.2 <- simulation(
  n = c(20, 20, 30),
  theta.X = c(0, 1),
  theta.Y = c(0, 1),
  theta.Z = c(0, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,2),
  tune.Z = c(0.2,2),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case3.3 <- simulation(
  n = c(20, 30, 50),
  theta.X = c(0, 1),
  theta.Y = c(0, 1),
  theta.Z = c(0, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = (c(0.1,12) + c(0.2,2)) / 2,
  tune.Y = (c(0.1,12) + c(0.1,2)) / 2,
  tune.Z = (c(0.1,5) + c(0.2,2)) / 2,
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case3.4 <- simulation(
  n = c(30, 50, 50),
  theta.X = c(0, 1),
  theta.Y = c(0, 1),
  theta.Z = c(0, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.1,12),
  tune.Y = c(0.1,12),
  tune.Z = c(0.1,5),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case3.5 <- simulation(
  n = c(50, 50, 50),
  theta.X = c(0, 1),
  theta.Y = c(0, 1),
  theta.Z = c(0, 1),
  distbn.X = rnorm,
  distbn.Y = rnorm,
  distbn.Z = rnorm,
  tune.X = c(0.1,12),
  tune.Y = c(0.1,12),
  tune.Z = c(0.1,5),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case4.1 <- simulation(
  n = c(20, 20, 20),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rlnorm,
  distbn.Y = rlnorm,
  distbn.Z = rlnorm,
  tune.X = c(0.2,4),
  tune.Y = c(1.5,1),
  tune.Z = c(10,0.8),
  u.var.inc = TRUE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case4.2 <- simulation(
  n = c(20, 20, 30),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rlnorm,
  distbn.Y = rlnorm,
  distbn.Z = rlnorm,
  tune.X = c(0.2,4),
  tune.Y = c(1.5,1),
  tune.Z = c(10,0.8),
  u.var.inc = TRUE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case4.3 <- simulation(
  n = c(20, 30, 50),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rlnorm,
  distbn.Y = rlnorm,
  distbn.Z = rlnorm,
  tune.X = (c(0.05,8) + c(0.2,4)) / 2,
  tune.Y = (c(0.8,8) + c(1.5,1)) / 2,
  tune.Z = (c(0.5,2) + c(10,0.8)) / 2,
  u.var.inc = FALSE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case4.4 <- simulation(
  n = c(30, 50, 50),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rlnorm,
  distbn.Y = rlnorm,
  distbn.Z = rlnorm,
  tune.X = c(0.05,8),
  tune.Y = c(0.8,8),
  tune.Z = c(0.5,2),
  u.var.inc = FALSE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case4.5 <- simulation(
  n = c(50, 50, 50),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rlnorm,
  distbn.Y = rlnorm,
  distbn.Z = rlnorm,
  tune.X = c(0.05,8),
  tune.Y = c(0.8,8),
  tune.Z = c(0.5,2),
  u.var.inc = FALSE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case5.1 <- simulation(
  n = c(20, 20, 20),
  theta.X = c(1, 1),
  theta.Y = c(2, 1),
  theta.Z = c(3, 1),
  distbn.X = rgamma,
  distbn.Y = rgamma,
  distbn.Z = rgamma,
  tune.X = c(0.2,3),
  tune.Y = c(0.5,3),
  tune.Z = c(1,1),
  u.var.inc = TRUE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case5.2 <- simulation(
  n = c(20, 20, 30),
  theta.X = c(1, 1),
  theta.Y = c(2, 1),
  theta.Z = c(3, 1),
  distbn.X = rgamma,
  distbn.Y = rgamma,
  distbn.Z = rgamma,
  tune.X = c(0.2,3),
  tune.Y = c(0.5,3),
  tune.Z = c(1,1),
  u.var.inc = TRUE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case5.3 <- simulation(
  n = c(20, 30, 50),
  theta.X = c(1, 1),
  theta.Y = c(2, 1),
  theta.Z = c(3, 1),
  distbn.X = rgamma,
  distbn.Y = rgamma,
  distbn.Z = rgamma,
  tune.X = (c(0.02,15) + c(0.2,3)) / 2,
  tune.Y = (c(0.05,8) + c(0.5,3)) / 2,
  tune.Z = (c(0.05,5) + c(1,1)) / 2,
  u.var.inc = FALSE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case5.4 <- simulation(
  n = c(30, 50, 50),
  theta.X = c(1, 1),
  theta.Y = c(2, 1),
  theta.Z = c(3, 1),
  distbn.X = rgamma,
  distbn.Y = rgamma,
  distbn.Z = rgamma,
  tune.X = c(0.02,15),
  tune.Y = c(0.05,8),
  tune.Z = c(0.05,5),
  u.var.inc = FALSE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case5.5 <- simulation(
  n = c(50, 50, 50),
  theta.X = c(1, 1),
  theta.Y = c(2, 1),
  theta.Z = c(3, 1),
  distbn.X = rgamma,
  distbn.Y = rgamma,
  distbn.Z = rgamma,
  tune.X = c(0.02,15),
  tune.Y = c(0.05,8),
  tune.Z = c(0.05,5),
  u.var.inc = FALSE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case6.1 <- simulation(
  n = c(20, 20, 20),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rnorm,
  distbn.Y = rgamma,
  distbn.Z = rlnorm,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case6.2 <- simulation(
  n = c(20, 20, 30),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rnorm,
  distbn.Y = rgamma,
  distbn.Z = rlnorm,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case6.3 <- simulation(
  n = c(20, 30, 50),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rnorm,
  distbn.Y = rgamma,
  distbn.Z = rlnorm,
  tune.X = (c(0.2,6) + c(0.2,2)) / 2,
  tune.Y = (c(0.01,15) + c(0.1,3)) / 2,
  tune.Z = (c(3,1.5) + c(10,0.5)) / 2,
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case6.4 <- simulation(
  n = c(30, 50, 50),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rnorm,
  distbn.Y = rgamma,
  distbn.Z = rlnorm,
  tune.X = c(0.2,6),
  tune.Y = c(0.01,15),
  tune.Z = c(3,1.5),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case6.5 <- simulation(
  n = c(50, 50, 50),
  theta.X = c(0, 1),
  theta.Y = c(1, 1),
  theta.Z = c(2, 1),
  distbn.X = rnorm,
  distbn.Y = rgamma,
  distbn.Z = rlnorm,
  tune.X = c(0.2,6),
  tune.Y = c(0.01,15),
  tune.Z = c(3,1.5),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

################################################################################
# SCENARIOS 7 TO 8
################################################################################

calc.mixture.vus <- function(theta.X, theta.Y, theta.Z, distbn){
  if (identical(distbn, rnorm)) {
    pdf.Y <- function(y) 0.5 * dnorm(y, theta.Y[1], theta.Y[2]) + 
      0.5 * dnorm(y, theta.Y[3], theta.Y[4])
    cdf.X <- function(x) 0.5 * pnorm(x, theta.X[1], theta.X[2]) + 
      0.5 * pnorm(x, theta.X[3], theta.X[4])
    cdf.Z <- function(z) 0.5 * pnorm(z, theta.Z[1], theta.Z[2]) + 
      0.5 * pnorm(z, theta.Z[3], theta.Z[4])
    
    lower <- -Inf
    upper <- Inf
  } else if (identical(distbn, rgamma)) {
    pdf.Y <- function(y) 0.5 * dgamma(y, theta.Y[1], theta.Y[2]) + 
      0.5 * dgamma(y, theta.Y[3], theta.Y[4])
    cdf.X <- function(x) 0.5 * pgamma(x, theta.X[1], theta.X[2]) + 
      0.5 * pgamma(x, theta.X[3], theta.X[4])
    cdf.Z <- function(z) 0.5 * pgamma(z, theta.Z[1], theta.Z[2]) + 
      0.5 * pgamma(z, theta.Z[3], theta.Z[4])
    
    lower <- 0
    upper <- Inf
  }
  integral <- function(y){
    cdf.X(y)*(1-cdf.Z(y))*pdf.Y(y)
  }
  vus <- integrate(integral, lower = lower, upper = upper)$value
  return(vus)
}

generate.mixture <- function(n, theta, distbn.1, distbn.2) {
  indicator <- rbinom(n, 1, 0.5)  # 0.5 probability for each component
  component.1 <- distbn.1(n, theta[1], theta[2])
  component.2 <- distbn.2(n, theta[3], theta[4])
  ifelse(indicator == 1, component.1, component.2)
}

simulation.mixture <- function(n, theta.X, theta.Y, theta.Z, distbn, tune.X, tune.Y, tune.Z, u.var.inc, bc.inc, n.sim = 50) {
  
  vus.true <- calc.mixture.vus(theta.X, theta.Y, theta.Z, distbn)
  
  if (bc.inc == TRUE){
    method.names <- c("EP", "K1", "MLE", "BC", "LEH", "MFPT", "MH")
    ci.method.names <- c("U", "MLE", "BC", "LEH",
                         "EP.Percent", "EP.BCa",
                         "K1.Percent", "K1.BCa",
                         "MLE.Percent", "MLE.BCa", "MLE.Student",
                         "BC.Percent", "BC.BCa", "BC.Student",
                         "LEH.Percent", "LEH.BCa", "LEH.Student",
                         "MFPT.Equi", "MFPT.HPD", "MH.Equi", "MH.HPD")
    acc.rate.names <- c("MFPT.mu1", "MFPT.mu2", "MFPT.mu3",
                        "MFPT.sigma1", "MFPT.sigma2", "MFPT.sigma3",
                        "MH.vus")
  }else{
    method.names <- c("EP", "K1", "MLE", "LEH", "MFPT", "MH")
    # method.names <- c("EP", "K1", "MLE", "MLE.Z", "LEH")
    ci.method.names <- c("U", "MLE", "LEH",
                         "EP.Percent", "EP.BCa",
                         "K1.Percent", "K1.BCa",
                         "MLE.Percent", "MLE.BCa", "MLE.Student",
                         "LEH.Percent", "LEH.BCa", "LEH.Student",
                         "MFPT.Equi", "MFPT.HPD", "MH.Equi", "MH.HPD")
    # ci.method.names <- c("U", "MLE", "MLE.Z", "LEH",
    #                      "EP.Percent", "EP.BCa",
    #                      "K1.Percent", "K1.BCa",
    #                      "MLE.Percent", "MLE.BCa", "MLE.Student",
    #                      "LEH.Percent", "LEH.BCa", "LEH.Student")
    acc.rate.names <- c("MFPT.mu1", "MFPT.mu2", "MFPT.mu3",
                        "MFPT.sigma1", "MFPT.sigma2", "MFPT.sigma3",
                        "MH.vus")
  }
  
  results <- list()
  for (i in 1:n.sim){
    X <- generate.mixture(n[1], theta.X, distbn, distbn)
    Y <- generate.mixture(n[2], theta.Y, distbn, distbn)
    Z <- generate.mixture(n[3], theta.Z, distbn, distbn)
    
    results[[i]] <- repetition(X, Y, Z, tune.X, tune.Y, tune.Z, u.var.inc, bc.inc)
    print(i)
  }
  
  vus.estimates <- matrix(NA, n.sim, length(method.names))
  colnames(vus.estimates) <- method.names
  
  acc.rates <- matrix(NA, n.sim, length(acc.rate.names))
  colnames(acc.rates) <- acc.rate.names
  
  if (bc.inc == TRUE){
    
    ci.estimates <- list(
      U = list(),
      MLE = list(),
      BC = list(),
      LEH = list(),
      EP.Percent = list(),
      EP.BCa = list(),
      K1.Percent = list(),
      K1.BCa = list(),
      MLE.Percent = list(),
      MLE.BCa = list(),
      MLE.Student = list(),
      BC.Percent = list(),
      BC.BCa = list(),
      BC.Student = list(),
      LEH.Percent = list(),
      LEH.BCa = list(),
      LEH.Student = list(),
      MFPT.Equi = list(),
      MFPT.HPD = list(),
      MH.Equi = list(),
      MH.HPD = list()
    )
    
    for(i in 1:n.sim) {
      vus.estimates[i, "EP"] <- results[[i]]$vus$ep
      vus.estimates[i, "K1"] <- results[[i]]$vus$k1
      vus.estimates[i, "MLE"] <- results[[i]]$vus$mle
      vus.estimates[i, "BC"] <- results[[i]]$vus$bc
      vus.estimates[i, "LEH"] <- results[[i]]$vus$leh
      vus.estimates[i, "MFPT"] <- results[[i]]$vus$mfpt
      vus.estimates[i, "MH"] <- results[[i]]$vus$mh
      ci.estimates$U[[i]] <- results[[i]]$asymptotic$u
      ci.estimates$MLE[[i]] <- results[[i]]$asymptotic$mle
      ci.estimates$BC[[i]] <- results[[i]]$asymptotic$bc
      ci.estimates$LEH[[i]] <- results[[i]]$asymptotic$leh
      ci.estimates$EP.Percent[[i]] <- results[[i]]$percent$ep
      ci.estimates$EP.BCa[[i]] <- results[[i]]$bca$ep
      ci.estimates$K1.Percent[[i]] <- results[[i]]$percent$k1
      ci.estimates$K1.BCa[[i]] <- results[[i]]$bca$k1
      ci.estimates$MLE.Percent[[i]] <- results[[i]]$percent$mle
      ci.estimates$MLE.BCa[[i]] <- results[[i]]$bca$mle
      ci.estimates$MLE.Student[[i]] <- results[[i]]$student$mle
      ci.estimates$BC.Percent[[i]] <- results[[i]]$percent$bc
      ci.estimates$BC.BCa[[i]] <- results[[i]]$bca$bc
      ci.estimates$BC.Student[[i]] <- results[[i]]$student$bc
      ci.estimates$LEH.Percent[[i]] <- results[[i]]$percent$leh
      ci.estimates$LEH.BCa[[i]] <- results[[i]]$bca$leh
      ci.estimates$LEH.Student[[i]] <- results[[i]]$student$leh
      ci.estimates$MFPT.Equi[[i]] <- results[[i]]$bayes$mfpt.equi
      ci.estimates$MFPT.HPD[[i]] <- results[[i]]$bayes$mfpt.hpd
      ci.estimates$MH.Equi[[i]] <- results[[i]]$bayes$mh.equi
      ci.estimates$MH.HPD[[i]] <- results[[i]]$bayes$mh.hpd
      acc.rates[i, "MFPT.mu1"] <- results[[i]]$acc$mfpt.mu1
      acc.rates[i, "MFPT.mu2"] <- results[[i]]$acc$mfpt.mu2
      acc.rates[i, "MFPT.mu3"] <- results[[i]]$acc$mfpt.mu3
      acc.rates[i, "MFPT.sigma1"] <- results[[i]]$acc$mfpt.sigma1
      acc.rates[i, "MFPT.sigma2"] <- results[[i]]$acc$mfpt.sigma2
      acc.rates[i, "MFPT.sigma3"] <- results[[i]]$acc$mfpt.sigma3
      acc.rates[i, "MH.vus"] <- results[[i]]$acc$mh
    }
  }else{
    
    ci.estimates <- list(
      U = list(),
      MLE = list(),
      # MLE.Z = list(),
      LEH = list(),
      EP.Percent = list(),
      EP.BCa = list(),
      K1.Percent = list(),
      K1.BCa = list(),
      MLE.Percent = list(),
      MLE.BCa = list(),
      MLE.Student = list(),
      LEH.Percent = list(),
      LEH.BCa = list(),
      LEH.Student = list(),
      MFPT.Equi = list(),
      MFPT.HPD = list(),
      MH.Equi = list(),
      MH.HPD = list()
    )
    
    for(i in 1:n.sim) {
      vus.estimates[i, "EP"] <- results[[i]]$vus$ep
      vus.estimates[i, "K1"] <- results[[i]]$vus$k1
      vus.estimates[i, "MLE"] <- results[[i]]$vus$mle
      # vus.estimates[i, "MLE.Z"] <- results[[i]]$vus$mle.z
      vus.estimates[i, "LEH"] <- results[[i]]$vus$leh
      ci.estimates$U[[i]] <- results[[i]]$asymptotic$u
      ci.estimates$MLE[[i]] <- results[[i]]$asymptotic$mle
      # ci.estimates$MLE.Z[[i]] <- results[[i]]$asymptotic$mle.z
      ci.estimates$LEH[[i]] <- results[[i]]$asymptotic$leh
      vus.estimates[i, "MFPT"] <- results[[i]]$vus$mfpt
      vus.estimates[i, "MH"] <- results[[i]]$vus$mh
      ci.estimates$EP.Percent[[i]] <- results[[i]]$percent$ep
      ci.estimates$EP.BCa[[i]] <- results[[i]]$bca$ep
      ci.estimates$K1.Percent[[i]] <- results[[i]]$percent$k1
      ci.estimates$K1.BCa[[i]] <- results[[i]]$bca$k1
      ci.estimates$MLE.Percent[[i]] <- results[[i]]$percent$mle
      ci.estimates$MLE.BCa[[i]] <- results[[i]]$bca$mle
      ci.estimates$MLE.Student[[i]] <- results[[i]]$student$mle
      ci.estimates$LEH.Percent[[i]] <- results[[i]]$percent$leh
      ci.estimates$LEH.BCa[[i]] <- results[[i]]$bca$leh
      ci.estimates$LEH.Student[[i]] <- results[[i]]$student$leh
      ci.estimates$MFPT.Equi[[i]] <- results[[i]]$bayes$mfpt.equi
      ci.estimates$MFPT.HPD[[i]] <- results[[i]]$bayes$mfpt.hpd
      ci.estimates$MH.Equi[[i]] <- results[[i]]$bayes$mh.equi
      ci.estimates$MH.HPD[[i]] <- results[[i]]$bayes$mh.hpd
      acc.rates[i, "MFPT.mu1"] <- results[[i]]$acc$mfpt.mu1
      acc.rates[i, "MFPT.mu2"] <- results[[i]]$acc$mfpt.mu2
      acc.rates[i, "MFPT.mu3"] <- results[[i]]$acc$mfpt.mu3
      acc.rates[i, "MFPT.sigma1"] <- results[[i]]$acc$mfpt.sigma1
      acc.rates[i, "MFPT.sigma2"] <- results[[i]]$acc$mfpt.sigma2
      acc.rates[i, "MFPT.sigma3"] <- results[[i]]$acc$mfpt.sigma3
      acc.rates[i, "MH.vus"] <- results[[i]]$acc$mh
    }
  }
  
  bias.result <- matrix(NA, 2, length(method.names))
  rmse.result <- matrix(NA, 2, length(method.names))
  colnames(bias.result) <- method.names
  colnames(rmse.result) <- method.names
  rownames(bias.result) <- c("Bias", "SE")
  rownames(rmse.result) <- c("RMSE", "SE")
  
  acc.rate.result <- numeric(length(acc.rate.names))
  names(acc.rate.result) <- acc.rate.names
  
  for (i in 1:length(method.names)) {
    vus.est <- vus.estimates[, i]
    
    bias.result[1, i] <- calc.bias(vus.est, vus.true, n.sim)$bias * 1000
    bias.result[2, i] <- calc.bias(vus.est, vus.true, n.sim)$se * 1000
    rmse.result[1, i] <- calc.rmse(vus.est, vus.true, n.sim)$rmse * 1000
    rmse.result[2, i] <- calc.rmse(vus.est, vus.true, n.sim)$se * 1000
  }
  
  for (i in 1:length(acc.rate.names)){
    acc.rate.result[i] <- mean(acc.rates[, i])
  }
  
  coverage.result <- matrix(NA, 1, length(ci.method.names))
  avglength.result <- matrix(NA, 1, length(ci.method.names))
  condlength.result <- matrix(NA, 1, length(ci.method.names))
  
  colnames(coverage.result) <- ci.method.names
  colnames(avglength.result) <- ci.method.names
  colnames(condlength.result) <- ci.method.names
  
  for (i in 1:length(ci.method.names)) {
    ci.list <- ci.estimates[[ci.method.names[i]]]
    coverage.result[1, i] <- calc.coverage(ci.list, vus.true)
    avglength.result[1, i] <- calc.avg.length(ci.list)
    condlength.result[1, i] <- calc.cdtional.length(ci.list, vus.true)
  }
  
  return(list(
    vus.true = vus.true,
    bias = bias.result,
    rmse = rmse.result,
    coverage = coverage.result,
    avg.length = avglength.result,
    cond.length = condlength.result,
    acc.rate = acc.rate.result
  ))
}

results.case7.1 <- simulation.mixture(
  n = c(20, 20, 20),
  theta.X = c(0, 1, 3, 1),
  theta.Y = c(1, 1, 4, 1.5),
  theta.Z = c(2, 1, 5, 2),
  distbn = rnorm,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case7.2 <- simulation.mixture(
  n = c(20, 20, 30),
  theta.X = c(0, 1, 3, 1),
  theta.Y = c(1, 1, 4, 1.5),
  theta.Z = c(2, 1, 5, 2),
  distbn = rnorm,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case7.3 <- simulation.mixture(
  n = c(20, 30, 50),
  theta.X = c(0, 1, 3, 1),
  theta.Y = c(1, 1, 4, 1.5),
  theta.Z = c(2, 1, 5, 2),
  distbn = rnorm,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case7.4 <- simulation.mixture(
  n = c(30, 50, 50),
  theta.X = c(0, 1, 3, 1),
  theta.Y = c(1, 1, 4, 1.5),
  theta.Z = c(2, 1, 5, 2),
  distbn = rnorm,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case7.5 <- simulation.mixture(
  n = c(50, 50, 50),
  theta.X = c(0, 1, 3, 1),
  theta.Y = c(1, 1, 4, 1.5),
  theta.Z = c(2, 1, 5, 2),
  distbn = rnorm,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case8.1 <- simulation.mixture(
  n = c(20, 20, 20),
  theta.X = c(1, 1, 4, 1),
  theta.Y = c(2, 1, 5, 1.5),
  theta.Z = c(3, 1, 6, 2),
  distbn = rgamma,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = TRUE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case8.2 <- simulation.mixture(
  n = c(20, 20, 30),
  theta.X = c(1, 1, 4, 1),
  theta.Y = c(2, 1, 5, 1.5),
  theta.Z = c(3, 1, 6, 2),
  distbn = rgamma,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = TRUE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case8.3 <- simulation.mixture(
  n = c(20, 30, 50),
  theta.X = c(1, 1, 4, 1),
  theta.Y = c(2, 1, 5, 1.5),
  theta.Z = c(3, 1, 6, 2),
  distbn = rgamma,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = FALSE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case8.4 <- simulation.mixture(
  n = c(30, 50, 50),
  theta.X = c(1, 1, 4, 1),
  theta.Y = c(2, 1, 5, 1.5),
  theta.Z = c(3, 1, 6, 2),
  distbn = rgamma,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = FALSE,
  bc.inc = TRUE,
  n.sim = 50
)

results.case8.5 <- simulation.mixture(
  n = c(50, 50, 50),
  theta.X = c(1, 1, 4, 1),
  theta.Y = c(2, 1, 5, 1.5),
  theta.Z = c(3, 1, 6, 2),
  distbn = rgamma,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = FALSE,
  bc.inc = TRUE,
  n.sim = 50
)

################################################################################
# SCENARIO 9
################################################################################

calc.normal.gamma.vus <- function(theta.X, theta.Y, theta.Z){
  # pdfs and cdfs
  pdf.Y <- function(y) 0.5 * dnorm(y, theta.Y[1], theta.Y[2]) + 
    0.5 * dgamma(y, theta.Y[3], theta.Y[4])
  cdf.X <- function(x) 0.5 * pnorm(x, theta.X[1], theta.X[2]) + 
    0.5 * pgamma(x, theta.X[3], theta.X[4])
  cdf.Z <- function(z) 0.5 * pnorm(z, theta.Z[1], theta.Z[2]) + 
    0.5 * pgamma(z, theta.Z[3], theta.Z[4])
  
  integral <- function(y){
    cdf.X(y)*(1-cdf.Z(y))*pdf.Y(y)
  }
  vus <- integrate(integral, lower = -Inf, upper = Inf)$value
  return(vus)
}

simulation.normal.gamma <- function(n, theta.X, theta.Y, theta.Z, distbn.1, distbn.2, tune.X, tune.Y, tune.Z, u.var.inc, bc.inc, n.sim = 50) {
  
  vus.true <- calc.normal.gamma.vus(theta.X, theta.Y, theta.Z)
  
  if (bc.inc == TRUE){
    method.names <- c("EP", "K1", "MLE", "BC", "LEH", "MFPT", "MH")
    ci.method.names <- c("U", "MLE", "BC", "LEH",
                         "EP.Percent", "EP.BCa",
                         "K1.Percent", "K1.BCa",
                         "MLE.Percent", "MLE.BCa", "MLE.Student",
                         "BC.Percent", "BC.BCa", "BC.Student",
                         "LEH.Percent", "LEH.BCa", "LEH.Student",
                         "MFPT.Equi", "MFPT.HPD", "MH.Equi", "MH.HPD")
    acc.rate.names <- c("MFPT.mu1", "MFPT.mu2", "MFPT.mu3",
                        "MFPT.sigma1", "MFPT.sigma2", "MFPT.sigma3",
                        "MH.vus")
  }else{
    method.names <- c("EP", "K1", "MLE", "LEH", "MFPT", "MH")
    # method.names <- c("EP", "K1", "MLE", "MLE.Z", "LEH")
    ci.method.names <- c("U", "MLE", "LEH",
                         "EP.Percent", "EP.BCa",
                         "K1.Percent", "K1.BCa",
                         "MLE.Percent", "MLE.BCa", "MLE.Student",
                         "LEH.Percent", "LEH.BCa", "LEH.Student",
                         "MFPT.Equi", "MFPT.HPD", "MH.Equi", "MH.HPD")
    # ci.method.names <- c("U", "MLE", "MLE.Z", "LEH",
    #                      "EP.Percent", "EP.BCa",
    #                      "K1.Percent", "K1.BCa",
    #                      "MLE.Percent", "MLE.BCa", "MLE.Student",
    #                      "LEH.Percent", "LEH.BCa", "LEH.Student")
    acc.rate.names <- c("MFPT.mu1", "MFPT.mu2", "MFPT.mu3",
                        "MFPT.sigma1", "MFPT.sigma2", "MFPT.sigma3",
                        "MH.vus")
  }
  
  results <- list()
  for (i in 1:n.sim){
    X <- generate.mixture(n[1], theta.X, distbn.1, distbn.2)
    Y <- generate.mixture(n[2], theta.Y, distbn.1, distbn.2)
    Z <- generate.mixture(n[3], theta.Z, distbn.1, distbn.2)
    
    results[[i]] <- repetition(X, Y, Z, tune.X, tune.Y, tune.Z, u.var.inc, bc.inc)
    print(i)
  }
  
  vus.estimates <- matrix(NA, n.sim, length(method.names))
  colnames(vus.estimates) <- method.names
  
  acc.rates <- matrix(NA, n.sim, length(acc.rate.names))
  colnames(acc.rates) <- acc.rate.names
  
  if (bc.inc == TRUE){
    
    ci.estimates <- list(
      U = list(),
      MLE = list(),
      BC = list(),
      LEH = list(),
      EP.Percent = list(),
      EP.BCa = list(),
      K1.Percent = list(),
      K1.BCa = list(),
      MLE.Percent = list(),
      MLE.BCa = list(),
      MLE.Student = list(),
      BC.Percent = list(),
      BC.BCa = list(),
      BC.Student = list(),
      LEH.Percent = list(),
      LEH.BCa = list(),
      LEH.Student = list(),
      MFPT.Equi = list(),
      MFPT.HPD = list(),
      MH.Equi = list(),
      MH.HPD = list()
    )
    
    for(i in 1:n.sim) {
      vus.estimates[i, "EP"] <- results[[i]]$vus$ep
      vus.estimates[i, "K1"] <- results[[i]]$vus$k1
      vus.estimates[i, "MLE"] <- results[[i]]$vus$mle
      vus.estimates[i, "BC"] <- results[[i]]$vus$bc
      vus.estimates[i, "LEH"] <- results[[i]]$vus$leh
      vus.estimates[i, "MFPT"] <- results[[i]]$vus$mfpt
      vus.estimates[i, "MH"] <- results[[i]]$vus$mh
      ci.estimates$U[[i]] <- results[[i]]$asymptotic$u
      ci.estimates$MLE[[i]] <- results[[i]]$asymptotic$mle
      ci.estimates$BC[[i]] <- results[[i]]$asymptotic$bc
      ci.estimates$LEH[[i]] <- results[[i]]$asymptotic$leh
      ci.estimates$EP.Percent[[i]] <- results[[i]]$percent$ep
      ci.estimates$EP.BCa[[i]] <- results[[i]]$bca$ep
      ci.estimates$K1.Percent[[i]] <- results[[i]]$percent$k1
      ci.estimates$K1.BCa[[i]] <- results[[i]]$bca$k1
      ci.estimates$MLE.Percent[[i]] <- results[[i]]$percent$mle
      ci.estimates$MLE.BCa[[i]] <- results[[i]]$bca$mle
      ci.estimates$MLE.Student[[i]] <- results[[i]]$student$mle
      ci.estimates$BC.Percent[[i]] <- results[[i]]$percent$bc
      ci.estimates$BC.BCa[[i]] <- results[[i]]$bca$bc
      ci.estimates$BC.Student[[i]] <- results[[i]]$student$bc
      ci.estimates$LEH.Percent[[i]] <- results[[i]]$percent$leh
      ci.estimates$LEH.BCa[[i]] <- results[[i]]$bca$leh
      ci.estimates$LEH.Student[[i]] <- results[[i]]$student$leh
      ci.estimates$MFPT.Equi[[i]] <- results[[i]]$bayes$mfpt.equi
      ci.estimates$MFPT.HPD[[i]] <- results[[i]]$bayes$mfpt.hpd
      ci.estimates$MH.Equi[[i]] <- results[[i]]$bayes$mh.equi
      ci.estimates$MH.HPD[[i]] <- results[[i]]$bayes$mh.hpd
      acc.rates[i, "MFPT.mu1"] <- results[[i]]$acc$mfpt.mu1
      acc.rates[i, "MFPT.mu2"] <- results[[i]]$acc$mfpt.mu2
      acc.rates[i, "MFPT.mu3"] <- results[[i]]$acc$mfpt.mu3
      acc.rates[i, "MFPT.sigma1"] <- results[[i]]$acc$mfpt.sigma1
      acc.rates[i, "MFPT.sigma2"] <- results[[i]]$acc$mfpt.sigma2
      acc.rates[i, "MFPT.sigma3"] <- results[[i]]$acc$mfpt.sigma3
      acc.rates[i, "MH.vus"] <- results[[i]]$acc$mh
    }
  }else{
    
    ci.estimates <- list(
      U = list(),
      MLE = list(),
      # MLE.Z = list(),
      LEH = list(),
      EP.Percent = list(),
      EP.BCa = list(),
      K1.Percent = list(),
      K1.BCa = list(),
      MLE.Percent = list(),
      MLE.BCa = list(),
      MLE.Student = list(),
      LEH.Percent = list(),
      LEH.BCa = list(),
      LEH.Student = list(),
      MFPT.Equi = list(),
      MFPT.HPD = list(),
      MH.Equi = list(),
      MH.HPD = list()
    )
    
    for(i in 1:n.sim) {
      vus.estimates[i, "EP"] <- results[[i]]$vus$ep
      vus.estimates[i, "K1"] <- results[[i]]$vus$k1
      vus.estimates[i, "MLE"] <- results[[i]]$vus$mle
      # vus.estimates[i, "MLE.Z"] <- results[[i]]$vus$mle.z
      vus.estimates[i, "LEH"] <- results[[i]]$vus$leh
      ci.estimates$U[[i]] <- results[[i]]$asymptotic$u
      ci.estimates$MLE[[i]] <- results[[i]]$asymptotic$mle
      # ci.estimates$MLE.Z[[i]] <- results[[i]]$asymptotic$mle.z
      ci.estimates$LEH[[i]] <- results[[i]]$asymptotic$leh
      vus.estimates[i, "MFPT"] <- results[[i]]$vus$mfpt
      vus.estimates[i, "MH"] <- results[[i]]$vus$mh
      ci.estimates$EP.Percent[[i]] <- results[[i]]$percent$ep
      ci.estimates$EP.BCa[[i]] <- results[[i]]$bca$ep
      ci.estimates$K1.Percent[[i]] <- results[[i]]$percent$k1
      ci.estimates$K1.BCa[[i]] <- results[[i]]$bca$k1
      ci.estimates$MLE.Percent[[i]] <- results[[i]]$percent$mle
      ci.estimates$MLE.BCa[[i]] <- results[[i]]$bca$mle
      ci.estimates$MLE.Student[[i]] <- results[[i]]$student$mle
      ci.estimates$LEH.Percent[[i]] <- results[[i]]$percent$leh
      ci.estimates$LEH.BCa[[i]] <- results[[i]]$bca$leh
      ci.estimates$LEH.Student[[i]] <- results[[i]]$student$leh
      ci.estimates$MFPT.Equi[[i]] <- results[[i]]$bayes$mfpt.equi
      ci.estimates$MFPT.HPD[[i]] <- results[[i]]$bayes$mfpt.hpd
      ci.estimates$MH.Equi[[i]] <- results[[i]]$bayes$mh.equi
      ci.estimates$MH.HPD[[i]] <- results[[i]]$bayes$mh.hpd
      acc.rates[i, "MFPT.mu1"] <- results[[i]]$acc$mfpt.mu1
      acc.rates[i, "MFPT.mu2"] <- results[[i]]$acc$mfpt.mu2
      acc.rates[i, "MFPT.mu3"] <- results[[i]]$acc$mfpt.mu3
      acc.rates[i, "MFPT.sigma1"] <- results[[i]]$acc$mfpt.sigma1
      acc.rates[i, "MFPT.sigma2"] <- results[[i]]$acc$mfpt.sigma2
      acc.rates[i, "MFPT.sigma3"] <- results[[i]]$acc$mfpt.sigma3
      acc.rates[i, "MH.vus"] <- results[[i]]$acc$mh
    }
  }
  
  bias.result <- matrix(NA, 2, length(method.names))
  rmse.result <- matrix(NA, 2, length(method.names))
  colnames(bias.result) <- method.names
  colnames(rmse.result) <- method.names
  rownames(bias.result) <- c("Bias", "SE")
  rownames(rmse.result) <- c("RMSE", "SE")
  
  acc.rate.result <- numeric(length(acc.rate.names))
  names(acc.rate.result) <- acc.rate.names
  
  for (i in 1:length(method.names)) {
    vus.est <- vus.estimates[, i]
    
    bias.result[1, i] <- calc.bias(vus.est, vus.true, n.sim)$bias * 1000
    bias.result[2, i] <- calc.bias(vus.est, vus.true, n.sim)$se * 1000
    rmse.result[1, i] <- calc.rmse(vus.est, vus.true, n.sim)$rmse * 1000
    rmse.result[2, i] <- calc.rmse(vus.est, vus.true, n.sim)$se * 1000
  }
  
  for (i in 1:length(acc.rate.names)){
    acc.rate.result[i] <- mean(acc.rates[, i])
  }
  
  coverage.result <- matrix(NA, 1, length(ci.method.names))
  avglength.result <- matrix(NA, 1, length(ci.method.names))
  condlength.result <- matrix(NA, 1, length(ci.method.names))
  
  colnames(coverage.result) <- ci.method.names
  colnames(avglength.result) <- ci.method.names
  colnames(condlength.result) <- ci.method.names
  
  for (i in 1:length(ci.method.names)) {
    ci.list <- ci.estimates[[ci.method.names[i]]]
    coverage.result[1, i] <- calc.coverage(ci.list, vus.true)
    avglength.result[1, i] <- calc.avg.length(ci.list)
    condlength.result[1, i] <- calc.cdtional.length(ci.list, vus.true)
  }
  
  return(list(
    vus.true = vus.true,
    bias = bias.result,
    rmse = rmse.result,
    coverage = coverage.result,
    avg.length = avglength.result,
    cond.length = condlength.result,
    acc.rate = acc.rate.result
  ))
}

results.case9.1 <- simulation.normal.gamma(
  n = c(20, 20, 20),
  theta.X = c(0, 1, 1, 1),
  theta.Y = c(1, 1, 2, 1),
  theta.Z = c(2, 1, 3, 1),
  distbn.1 = rnorm,
  distbn.2 = rgamma,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case9.2 <- simulation.normal.gamma(
  n = c(20, 20, 30),
  theta.X = c(0, 1, 1, 1),
  theta.Y = c(1, 1, 2, 1),
  theta.Z = c(2, 1, 3, 1),
  distbn.1 = rnorm,
  distbn.2 = rgamma,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = TRUE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case9.3 <- simulation.normal.gamma(
  n = c(20, 30, 50),
  theta.X = c(0, 1, 1, 1),
  theta.Y = c(1, 1, 2, 1),
  theta.Z = c(2, 1, 3, 1),
  distbn.1 = rnorm,
  distbn.2 = rgamma,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case9.4 <- simulation.normal.gamma(
  n = c(30, 50, 50),
  theta.X = c(0, 1, 1, 1),
  theta.Y = c(1, 1, 2, 1),
  theta.Z = c(2, 1, 3, 1),
  distbn.1 = rnorm,
  distbn.2 = rgamma,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)

results.case9.5 <- simulation.normal.gamma(
  n = c(50, 50, 50),
  theta.X = c(0, 1, 1, 1),
  theta.Y = c(1, 1, 2, 1),
  theta.Z = c(2, 1, 3, 1),
  distbn.1 = rnorm,
  distbn.2 = rgamma,
  tune.X = c(0.2,2),
  tune.Y = c(0.1,3),
  tune.Z = c(10,0.5),
  u.var.inc = FALSE,
  bc.inc = FALSE,
  n.sim = 50
)
