################################################################################
# EMPIRICAL / MANN-WHITNEY
################################################################################

ep.boot <- function(data, indices){
  resampled.data <- data[indices,]
  X <- resampled.data[resampled.data$class == 1, "marker"]
  Y <- resampled.data[resampled.data$class == 2, "marker"]
  Z <- resampled.data[resampled.data$class == 3, "marker"]
  return(emp.vus(X, Y, Z, dat = NULL, old.version = TRUE))
}

u.var <- function(X, Y, Z, vus) {
  n1 <- length(X)
  n2 <- length(Y)
  n3 <- length(Z)
  
  indicator <- function(X, Y, Z) {
    X.less.Y <- ifelse(X < Y, 1, ifelse(X == Y, 0.5, 0))
    Y.less.Z <- ifelse(Y < Z, 1, ifelse(Y == Z, 0.5, 0))
    return(X.less.Y * Y.less.Z)
  }
  
  q12 <- 0
  for(i in 1:n1) {
    for(I in 1:n1) {
      if(I != i) {
        for(j in 1:n2) {
          for(J in 1:n2) {
            if(J != j) {
              for(k in 1:n3) {
                q12 <- q12 + (indicator(X[i], Y[j], Z[k]) * indicator(X[I], Y[J], Z[k]))
              }
            }
          }
        }
      }
    }
  }
  q12 <- q12 / (n1 * (n1 - 1) * n2 * (n2 - 1) * n3)
  
  
  q23 <- 0
  for(i in 1:n1) {
    for(j in 1:n2) {
      for(J in 1:n2) {
        if(J != j) {
          for(k in 1:n3) {
            for(K in 1:n3) {
              if(K != k) {
                q23 <- q23 + (indicator(X[i], Y[j], Z[k]) * indicator(X[i], Y[J], Z[K]))
              }
            }
          }
        }
      }
    }
  }
  q23 <- q23 / (n1 * n2 * (n2 - 1) * n3 * (n3 - 1))
  
  q13 <- 0
  for(i in 1:n1) {
    for(I in 1:n1) {
      if(I != i) {
        for(j in 1:n2) {
          for(k in 1:n3) {
            for(K in 1:n3) {
              if(K != k) {
                q13 <- q13 + (indicator(X[i], Y[j], Z[k]) * indicator(X[I], Y[j], Z[K]))
              }
            }
          }
        }
      }
    }
  }
  q13 <- q13 / (n1 * (n1 - 1) * n2 * n3 * (n3 - 1))
  
  q3 <- 0
  for(i in 1:n1) {
    for(j in 1:n2) {
      for(k in 1:n3) {
        for(K in 1:n3) {
          if(K != k) {
            q3 <- q3 + (indicator(X[i], Y[j], Z[k]) * indicator(X[i], Y[j], Z[K]))
          }
        }
      }
    }
  }
  q3 <- q3 / (n1 * n2 * n3 * (n3 - 1))
  
  q2 <- 0
  for(i in 1:n1) {
    for(j in 1:n2) {
      for(J in 1:n2) {
        if(J != j) {
          for(k in 1:n3) {
            q2 <- q2 + (indicator(X[i], Y[j], Z[k]) * indicator(X[i], Y[J], Z[k]))
          }
        }
      }
    }
  }
  q2 <- q2 / (n1 * n2 * (n2 - 1) * n3)
  
  q1 <- 0
  for(i in 1:n1) {
    for(I in 1:n1) {
      if(I != i) {
        for(j in 1:n2) {
          for(k in 1:n3) {
            q1 <- q1 + (indicator(X[i], Y[j], Z[k]) * indicator(X[I], Y[j], Z[k]))
          }
        }
      }
    }
  }
  q1 <- q1 / (n1 * (n1 - 1) * n2 * n3)
  
  var <- (vus * (1 - vus) +
            (n3 - 1) * (q3 - vus^2) +
            (n2 - 1) * (q2 - vus^2) +
            (n1 - 1) * (q1 - vus^2) +
            (n1 - 1) * (n2 - 1) * (q12 - vus^2) +
            (n2 - 1) * (n3 - 1) * (q23 - vus^2) +
            (n1 - 1) * (n3 - 1) * (q13 - vus^2)) / (n1 * n2 * n3) 
  
  s <- sqrt(var)
  ci <- c(vus - 1.96 * s, vus + 1.96 * s)
  
  return(list(sd = s, ci = ci))
}

################################################################################
# KERNEL SMOOTHING
################################################################################

bandwidth <- function(data) { 
  n <- length(data)
  s <- sd(data)
  iqr <- IQR(data)
  h <- (4/(3*n))^(1/5) * min(s, iqr/1.349)
  return(h)
}

f.hat <-  function(data, h, t) {
  mean(dnorm((t-data)/h)/h)
}

F.hat <- function(data, h, t) { 
  mean(pnorm((t-data)/h))
}

k1.method <- function(X, Y, Z){
  
  h1 <- bandwidth(X)
  h2 <- bandwidth(Y)
  h3 <- bandwidth(Z)
  
  integrand <- function(y) {
    F.hat(X, h1, y)*(1-F.hat(Z, h3, y))*f.hat(Y, h2, y)
  }
  vus <- integrate(Vectorize(integrand), -Inf, Inf)$value
  
  return(vus)
}

F.hat.inv <- function(data, h, p) {
  uniroot(function(x) F.hat(data, h, x) - p,
          lower = min(data)-3*bandwidth(data),
          upper = max(data)+3*bandwidth(data))$root
}

k1.roc.point <- function(X, Y, Z, p1, p3) {
  
  h1 <- bandwidth(X)
  h2 <- bandwidth(Y)
  h3 <- bandwidth(Z)
  
  F.1.inv <- F.hat.inv(X,h1,p1)
  F.3.inv <- F.hat.inv(Z,h3,1-p3)
  
  if (F.1.inv > F.3.inv){
    roc.point <- 0
  } else{
    roc.point <- F.hat(Y, h2, F.3.inv) - F.hat(Y,h2,F.1.inv)
  }
  return(roc.point)
}

k1.roc.points <- function(X, Y, Z){
  p1=seq(0.001,0.99,len=20)
  p3=seq(0.001,0.99,len=20)
  
  roc.points <- matrix(0, 20, 20)
  
  for (i in 1:20){
    for (j in 1:20) {
      roc.points[i,j] <- k1.roc.point(X,Y,Z,p1[i],p3[j])
    }
  }
  return(roc.points)
}

k1.boot <- function(data, indices){
  resampled.data <- data[indices,]
  X <- resampled.data[resampled.data$class == 1, "marker"]
  Y <- resampled.data[resampled.data$class == 2, "marker"]
  Z <- resampled.data[resampled.data$class == 3, "marker"]
  return(k1.method(X, Y, Z))
}

################################################################################
# MLE
################################################################################

mle.boot <- function(data, indices){
  resampled.data <- data[indices,]
  X <- resampled.data[resampled.data$class == 1, "marker"]
  Y <- resampled.data[resampled.data$class == 2, "marker"]
  Z <- resampled.data[resampled.data$class == 3, "marker"]
  out <- trinVUS.test(X, Y, Z)
  return(c(out$estimate, sqrt(out$Sigma[1,1])^2))
}

z.transformation <- function(vus) log((1+vus)/(1-vus))
variance <- function(vus, var.vus) (4*var.vus)/((1-vus^2)^2)
z.inverse <- function(theta) (exp(theta)-1)/(exp(theta)+1)

mle.z.boot <- function(data, indices){
  resampled.data <- data[indices,]
  X <- resampled.data[resampled.data$class == 1, "marker"]
  Y <- resampled.data[resampled.data$class == 2, "marker"]
  Z <- resampled.data[resampled.data$class == 3, "marker"]
  out <- trinVUS.test(X, Y, Z)
  theta <- z.transformation(out$estimate)
  var.theta <- variance(out$estimate, out$Sigma[1,1])
  return(c(theta, sqrt(var.theta)^2))
}

################################################################################
# BOX-COX
################################################################################

bc.boot <- function(data, indices){
  resampled.data <- data[indices,]
  X <- resampled.data[resampled.data$class == 1, "marker"]
  Y <- resampled.data[resampled.data$class == 2, "marker"]
  Z <- resampled.data[resampled.data$class == 3, "marker"]
  bc.XYZ <- boxcoxROC(X, Y, Z, verbose = FALSE)
  out <- trinVUS.test(bc.XYZ$xbc, bc.XYZ$ybc, bc.XYZ$zbc)
  return(c(out$estimate, sqrt(out$Sigma[1,1])^2))
}

bayes.cred <- function(vus.chain, alpha = 0.05) {
  ci.equi <- quantile(vus.chain, probs = c(alpha/2, 1 - alpha/2))
  
  sorted.vus <- sort(vus.chain)
  n <- length(sorted.vus)
  n.in.interval <- ceiling((1 - alpha) * n)
  interval.widths <- sorted.vus[n.in.interval:n] - sorted.vus[1:(n - n.in.interval + 1)]
  min.width <- which.min(interval.widths)
  ci.hpd <- c(sorted.vus[min.width], 
              sorted.vus[min.width + n.in.interval - 1])
  
  return(list(
    mean = mean(vus.chain),
    sd = sd(vus.chain),
    equi = ci.equi,
    hpd = ci.hpd
  ))
}

################################################################################
# LI AND ZHOU'S METHOD 1
################################################################################

library(nleqslv)
library(numDeriv)

bandwidth <- function(data) { 
  n <- length(data)
  s <- sd(data)
  iqr <- IQR(data)
  h <- (4/(3*n))^(1/5) * min(s, iqr/1.349)
  return(h)
}

f.hat <-  function(data, h, t) {
  mean(dnorm((t-data)/h)/h)
}

F.hat <- function(data, h, t) { 
  mean(pnorm((t-data)/h))
}

m1.method <- function(X, Y, Z){
  n1 <- length(X)
  n2 <- length(Y)
  n3 <- length(Z)
  
  mu1 <- mean(X)
  sigma1 <- sd(X)
  mu2 <- mean(Y)
  sigma2 <- sd(Y)
  mu3 <- mean(Z)
  sigma3 <- sd(Z)
  
  alpha1 = (mu3 - mu2) / sigma2
  alpha2 = sigma3 / sigma2
  alpha3 = (mu1 - mu2) / sigma2
  alpha4 = sigma1 / sigma2
  alpha.init <- c(alpha1, alpha2, alpha3, alpha4)
  
  N1 <- 20
  N3 <- 20
  p1 <- seq(0.001, 0.99, length.out = N1)
  p3 <- seq(0.001, 0.99, length.out = N3)
  
  F1 <- ecdf(X)
  F2 <- ecdf(Y)
  F3 <- ecdf(Z)
  quantile1 <- quantile(X, p1)
  quantile3 <- quantile(Z, 1-p3)
  
  roc.emp.grid <- matrix(0, N1, N3)
  for (i in 1:N1) {
    for (j in 1:N3) {
      if (quantile1[i] <= quantile3[j]) {
        roc.emp.grid[i, j] <- F2(quantile3[j]) - F2(quantile1[i])
      }
    }
  }
  
  method.1 <- function(alpha) {
    linear.eq <- matrix(0, nrow = 4, ncol = 1)
    r1 <- dnorm(alpha[1]+alpha[2]*qnorm(1-p3))
    r2 <- dnorm(alpha[1]+alpha[2]*qnorm(1-p3))*qnorm(1-p3)
    r3 <- -dnorm(alpha[3]+alpha[4]*qnorm(p1))
    r4 <- -dnorm(alpha[3]+alpha[4]*qnorm(p1))*qnorm(p1)
    for (i in 1:N1) {
      for (j in 1:N3) {
        if ((alpha[1] + alpha[2]*qnorm(1-p3[j])) > (alpha[3] + alpha[4]*qnorm(p1[i]))){
          roc.para <- pnorm(alpha[1] + alpha[2]*qnorm(1-p3[j])) - pnorm(alpha[3] + alpha[4]*qnorm(p1[i]))
        } else{
          roc.para <- 0
        }
        roc.emp <- roc.emp.grid[i, j]          
        roc.deriv <- matrix(c(r1[j], r2[j], r3[i], r4[i]), ncol=1)
        linear.eq <- linear.eq + roc.deriv * (roc.emp - roc.para)
      }
    }
    return(linear.eq)
  }
  
  solution <- nleqslv(alpha.init, method.1, method = "Newton", control = list(maxit = 100, trace = 0))
  
  alpha.hat <- solution$x
  roc <- function(prob1, prob3) { 
    ifelse((alpha.hat[1] + alpha.hat[2]*qnorm(1-prob3)) >= (alpha.hat[3] + alpha.hat[4]*qnorm(prob1)),
           pnorm(alpha.hat[1] + alpha.hat[2]*qnorm(1-prob3)) - pnorm(alpha.hat[3] + alpha.hat[4]*qnorm(prob1)),
           0)
  }
  integrand <- function(prob3) {
    integrate(function(prob1) roc(prob1, prob3), 0.001, 0.99)$value
  }
  vus <- integrate(Vectorize(integrand), 0.001, 0.99)$value
  
  r1.hat <- dnorm(alpha.hat[1]+alpha.hat[2]*qnorm(1-p3))
  r2.hat <- dnorm(alpha.hat[1]+alpha.hat[2]*qnorm(1-p3))*qnorm(1-p3)
  r3.hat <- -dnorm(alpha.hat[3]+alpha.hat[4]*qnorm(p1))
  r4.hat <- -dnorm(alpha.hat[3]+alpha.hat[4]*qnorm(p1))*qnorm(p1)
  
  roc.para.grid <- matrix(0, N1, N3)
  roc.deriv.grid <- matrix(0, nrow=N1*N3, ncol=4)
  for (i in 1:N1) {
    for (j in 1:N3) {
      roc.para.grid[i,j] <- roc(p1[i], p3[j])
      roc.deriv.grid[(i-1)*N3+j,] <- c(r1.hat[j], r2.hat[j], r3.hat[i], r4.hat[i])
    }
  }
  
  var.alpha <- function(FX, FY, FZ, fX, fY, fZ, FX.inv, FZ.inv){
    
    p3.inv <- FZ.inv(1-p3)
    p1.inv <- FX.inv(p1)
    fY.p3 <- fY(p3.inv)
    fZ.p3 <- fZ(p3.inv)
    fY.p1 <- fY(p1.inv)
    fX.p1 <- fX(p1.inv)
    
    cov1 <- matrix(0, N1 * N3, N1 * N3)
    cov2 <- matrix(0, N1 * N3, N1 * N3)
    cov3 <- matrix(0, N1 * N3, N1 * N3)
    
    for (i in 1:N1){
      for (j in 1:N3){
        idx1 <- (i-1)*N3 + j
        
        for (I in 1:N1){
          for (J in 1:N3){
            idx2 <- (I-1)*N3 + J
            
            Q1 <- roc.para.grid[i,j]
            Q2 <- roc.para.grid[I,J]
            cov1[idx1, idx2] <- min(Q1, Q2) - Q1 * Q2
            
            if (j == J) {
              cov2[idx1, idx2] <- (fY.p3[j] * fY.p3[J] / (fZ.p3[j] * fZ.p3[J])) * 
                (1-max(p3[j], p3[J]) - (1-p3[j])*(1-p3[J]))
            }
            
            if (i == I) {
              cov3[idx1, idx2] <- (fY.p1[i] * fY.p1[I] / (fX.p1[i] * fX.p1[I])) * 
                (min(p1[i], p1[I]) - p1[i]*p1[I])
            }
          }
        }
      }
    }
    
    lambda1 <- (n1 * n3) / (n1 * n2 + n1 * n3 + n2 * n3)
    lambda2 <- (n1 * n2) / (n1 * n2 + n1 * n3 + n2 * n3)
    lambda3 <- (n2 * n3) / (n1 * n2 + n1 * n3 + n2 * n3)
    cov <- lambda1*cov1 + lambda2*cov2 + lambda3*cov3
    
    # sandwich variance
    D <- t(roc.deriv.grid) %*% roc.deriv.grid
    
    eigenvalues <- eigen(D, only.values = TRUE)$values
    
    if (any(abs(eigenvalues) < 1e-10)) {
      lambda <- 0.001
      D <- D + lambda * diag(nrow(D))
    }
    
    D.inv <- solve(D)
    a.var <- D.inv %*% t(roc.deriv.grid) %*% cov %*% roc.deriv.grid %*% D.inv
    
    a <- (n1*n2*n3)/(n1*n2+n2*n3+n1*n3)
    var <- a.var/a
    
    return(var)
  }
  
  F1 <- function(t) pnorm(t, mu1, sigma1)
  F2 <- function(t) pnorm(t, mu2, sigma2)
  F3 <- function(t) pnorm(t, mu3, sigma3)
  f1 <- function(t) dnorm(t, mu1, sigma1)
  f2 <- function(t) dnorm(t, mu2, sigma2)
  f3 <- function(t) dnorm(t, mu3, sigma3)
  F1.inv <- function(p) qnorm(p, mu1, sigma1)
  F3.inv <- function(p) qnorm(p, mu3, sigma3)
  
  var.alpha.model <- var.alpha(F1, F2, F3, f1, f2, f3, F1.inv, F3.inv)
  
  h1 <- bandwidth(X)
  h2 <- bandwidth(Y)
  h3 <- bandwidth(Z)
  
  F1 <- function(x) sapply(x, function(t) F.hat(X, h1, t))
  F2 <- function(x) sapply(x, function(t) F.hat(Y, h2, t))
  F3 <- function(x) sapply(x, function(t) F.hat(Z, h3, t))
  f1 <- function(x) sapply(x, function(t) f.hat(X, h1, t))
  f2 <- function(x) sapply(x, function(t) f.hat(Y, h2, t))
  f3 <- function(x) sapply(x, function(t) f.hat(Z, h3, t))
  F1.inv <- function(p) quantile(X, p)
  F3.inv <- function(p) quantile(Z, p)
  
  var.alpha.emp <- var.alpha(F1, F2, F3, f1, f2, f3, F1.inv, F3.inv)
  
  vus.grad <- grad(function(alpha) {
    roc.parametric <- function(prob1, prob3) {
      ifelse((alpha[1] + alpha[2]*qnorm(1-prob3)) > (alpha[3] + alpha[4]*qnorm(prob1)),
             pnorm(alpha[1] + alpha[2]*qnorm(1-prob3)) - pnorm(alpha[3] + alpha[4]*qnorm(prob1)),
             0)
    }
    integrand <- function(prob3) {
      integrate(function(prob1) roc.parametric(prob1, prob3), 0.001, 0.99)$value
    }
    integrate(Vectorize(integrand), 0.001, 0.99)$value
  }, alpha.hat)
  
  var.model <- as.numeric(t(vus.grad) %*% var.alpha.model %*% vus.grad)
  s.model = sqrt(var.model)
  ci.model = c(vus - 1.96*s.model, vus + 1.96*s.model)
  
  var.emp <- as.numeric(t(vus.grad) %*% var.alpha.emp %*% vus.grad)
  s.emp = sqrt(var.emp)
  ci.emp = c(vus - 1.96*s.emp, vus + 1.96*s.emp)
  
  return(list(roc = roc,
              alpha = alpha.hat,
              vus = vus,
              sd.model = s.model,
              ci.model = ci.model,
              sd.emp = s.emp,
              ci.emp = ci.emp))
}

################################################################################
# LEHMANN
################################################################################

library(survival)
leh.method <- function(X, Y, Z) {
  n1 <- length(X)
  n2 <- length(Y)
  n3 <- length(Z)
  
  leh.data <- data.frame(
    D = c(rep(1,n1), rep(2, n2), rep(3,n3)),
    marker = c(X, Y, Z)
  )
  leh.data$D1 <- ifelse(leh.data$D >= 2, 1, 0)
  leh.data$D2 <- ifelse(leh.data$D == 3, 1, 0)
  leh.data$event <- 1
  
  cox.fit <- coxph(Surv(marker, event) ~ D1 + D2, leh.data)
  beta.hat <- coef(cox.fit)
  theta.hat1 <- exp(beta.hat["D1"])
  theta.hat2 <- exp(beta.hat["D2"])
  
  roc <- function(q1, q2) {
    ifelse(0 < q2 & q2 < (1-q1)^(theta.hat1*theta.hat2),
           (1 - q1)^theta.hat1 - q2^(1/theta.hat2),
           0)
  }
  
  vus <- 1 / ((1 + theta.hat2) * (1 + theta.hat1 * (1 + theta.hat2)))
  
  vcov.beta <- vcov(cox.fit)
  var.theta1 <- exp(2 * beta.hat["D1"]) * vcov.beta["D1", "D1"]
  var.theta2 <- exp(2 * beta.hat["D2"]) * vcov.beta["D2", "D2"]
  cov.theta <- exp(beta.hat["D1"] + beta.hat["D2"]) * vcov.beta["D1", "D2"]
  
  partial.theta1 <- -1 / ((1 + theta.hat1 * (1 + theta.hat2))^2)
  partial.theta2 <- -(1 + 2*theta.hat1*(1 + theta.hat2)) / ((1 + theta.hat2)^2 * (1 + theta.hat1*(1 + theta.hat2))^2)
  var <- (partial.theta1^2 * var.theta1 + partial.theta2^2 * var.theta2 + 2 * partial.theta1 * partial.theta2 * cov.theta)
  
  s <- sqrt(var)
  ci <- c(vus - 1.96 * s, vus + 1.96 * s)
  
  return(list(roc = roc, vus = vus, sd = s, var = var, ci = ci))
}

leh.boot <- function(data, indices){
  result <- tryCatch({
    resampled.data <- data[indices,]
    X <- resampled.data[resampled.data$class == 1, "marker"]
    Y <- resampled.data[resampled.data$class == 2, "marker"]
    Z <- resampled.data[resampled.data$class == 3, "marker"]
    out <- leh.method(X,Y,Z)
    c(out$vus, out$sd^2)
  }, error = function(e) {
    c(NA, NA)
  })
  return(result)
}

################################################################################
# MFPT
################################################################################

#Set at level J countaing y 
label<-function(y,mu,sigma){
  label<-0
  n<-length(y)  
  for(i in 1:n){
    if(floor(pnorm(y[i],mu,sigma)*16+1)>16.9999) label[i]<-16
    else label[i]<-floor(pnorm(y[i],mu,sigma)*16+1)
  }
  return(label)
}

#Count data in intervals at level j=1,2,3,4
counter1<-function(y,mu,sigma){
  q1<-qnorm(0.5,mu,sigma)
  count<-numeric(2)
  for(i in 1:length(y)){
    if(y[i]<=q1) count[1]<-count[1]+1
    else count[2]<-count[2]+1
  }  
  return(count)
}

counter2<-function(y,mu,sigma){
  q1<-qnorm(0.25,mu,sigma)	
  q2<-qnorm(0.5,mu,sigma)
  q3<-qnorm(0.75,mu,sigma)
  count<-numeric(4)
  n<-length(y)
  for(i in 1:n){
    if(y[i]<=q1) count[1]<-count[1]+1
    else
      if(y[i]>q1 & y[i]<=q2) count[2]<-count[2]+1
      else
        if(y[i]>q2 & y[i]<=q3) count[3]<-count[3]+1
        else count[4]<-count[4]+1	
  }
  return(count)
}

counter3<-function(y,mu,sigma){
  count<-numeric(8)
  q1<-qnorm(0.125,mu,sigma)
  q2<-qnorm(0.25,mu,sigma)	
  q3<-qnorm(0.375,mu,sigma)
  q4<-qnorm(0.5,mu,sigma)
  q5<-qnorm(0.625,mu,sigma)
  q6<-qnorm(0.75,mu,sigma)
  q7<-qnorm(0.875,mu,sigma)
  n<-length(y)
  for(i in 1:n){
    if(y[i]<=q1) count[1]<-count[1]+1	
    else
      if(y[i]>q1 & y[i]<=q2) count[2]<-count[2]+1
      else
        if(y[i]>q2 & y[i]<=q3) count[3]<-count[3]+1
        else
          if(y[i]>q3 & y[i]<=q4) count[4]<-count[4]+1
          else
            if(y[i]>q4 & y[i]<=q5) count[5]<-count[5]+1
            else
              if(y[i]>q5 & y[i]<=q6) count[6]<-count[6]+1
              else
                if(y[i]>q6 & y[i]<=q7) count[7]<-count[7]+1
                else count[8]<-count[8]+1
  }
  return(count)
}

counter4<-function(y,mu,sigma){
  count<-numeric(16)
  n<-length(y)
  q1<-qnorm(0.0625,mu,sigma)
  q2<-qnorm(0.125,mu,sigma)
  q3<-qnorm(0.1875,mu,sigma)	
  q4<-qnorm(0.25,mu,sigma)
  q5<-qnorm(0.3125,mu,sigma)
  q6<-qnorm(0.375,mu,sigma)
  q7<-qnorm(0.4375,mu,sigma)
  q8<-qnorm(0.5,mu,sigma)
  q9<-qnorm(0.5625,mu,sigma)
  q10<-qnorm(0.625,mu,sigma)
  q11<-qnorm(0.6875,mu,sigma)
  q12<-qnorm(0.75,mu,sigma)
  q13<-qnorm(0.8125,mu,sigma)
  q14<-qnorm(0.8750,mu,sigma)
  q15<-qnorm(0.9375,mu,sigma)
  for(i in 1:n){
    if(y[i]<=q1)	count[1]<-count[1]+1
    else
      if(y[i]>q1 & y[i]<=q2) count[2]<-count[2]+1
      else
        if(y[i]>q2 & y[i]<=q3) count[3]<-count[3]+1
        else
          if(y[i]>q3 & y[i]<=q4) count[4]<-count[4]+1
          else
            if(y[i]>q4 & y[i]<=q5) count[5]<-count[5]+1
            else
              if(y[i]>q5 & y[i]<=q6) count[6]<-count[6]+1
              else
                if(y[i]>q6 & y[i]<=q7) count[7]<-count[7]+1
                else
                  if(y[i]>q7 & y[i]<=q8) count[8]<-count[8]+1
                  else
                    if(y[i]>q8 & y[i]<=q9) count[9]<-count[9]+1
                    else
                      if(y[i]>q9 & y[i]<=q10) count[10]<-count[10]+1
                      else
                        if(y[i]>q10 & y[i]<=q11) count[11]<-count[11]+1
                        else
                          if(y[i]>q11 & y[i]<=q12) count[12]<-count[12]+1
                          else
                            if(y[i]>q12 & y[i]<=q13) count[13]<-count[13]+1
                            else
                              if(y[i]>q13 & y[i]<=q14) count[14]<-count[14]+1
                              else
                                if(y[i]>q14 & y[i]<=q15) count[15]<-count[15]+1
                                else count[16]<-count[16]+1
  }
  return(count)
}

#Probabilities down branches
probs.down<-function(x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416){
  p1<-x11*x21*x31*x41
  p2<-x11*x21*x31*x42
  p3<-x11*x21*x32*x43
  p4<-x11*x21*x32*x44
  p5<-x11*x22*x33*x45
  p6<-x11*x22*x33*x46
  p7<-x11*x22*x34*x47
  p8<-x11*x22*x34*x48
  p9<-x12*x23*x35*x49
  p10<-x12*x23*x35*x410
  p11<-x12*x23*x36*x411
  p12<-x12*x23*x36*x412
  p13<-x12*x24*x37*x413
  p14<-x12*x24*x37*x414
  p15<-x12*x24*x38*x415	
  p16<-x12*x24*x38*x416
  return(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16))
}

#Product of branch probabilities leading to interval containing y at level J
probs.label<-function(y,mu,sigma,x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416){
  p<-0
  p<-probs.down(x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416)[label(y,mu,sigma)]
  return(p)
}

#probs.label applied to all y in data
probs.data<-function(y,mu,sigma,x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416){
  probs.num<-0
  n<-length(y)
  for(i in 1:n){
    probs.num[i]<-as.numeric(probs.label(y[i],mu,sigma,x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416))
  }	
  return(probs.num)
}

#MCMC scheme
run.mfpt.single <- function(y, iterations, burn.in, thin, tune1, tune2, c=1){
  
  mu<-numeric(iterations)
  sigma<-numeric(iterations)
  x11<-x12<-numeric(iterations)
  x21<-x22<-x23<-x24<-numeric(iterations)
  x31<-x32<-x33<-x34<-x35<-x36<-x37<-x38<-numeric(iterations)
  x41<-x42<-x43<-x44<-x45<-x46<-x47<-x48<-x49<-x410<-x411<-x412<-x413<-x414<-x415<-x416<-numeric(iterations)
  
  current.mu<-mean(y)
  current.sigma<-sd(y)
  accept.mu<-0
  accept.sigma<-0
  n<-length(y)
  
  for(t in 1:iterations){
    
    counts1=counter1(y,current.mu,current.sigma)  
    counts2=counter2(y,current.mu,current.sigma)  
    counts3=counter3(y,current.mu,current.sigma)  
    counts4=counter4(y,current.mu,current.sigma)  
    
    x11[t]<-rbeta(1,c+counts1[1],c+counts1[2])  
    x12[t]<-1-x11[t]
    
    x21[t]<-rbeta(1,4*c+counts2[1],4*c+counts2[2])
    x22[t]<-1-x21[t]
    x23[t]<-rbeta(1,4*c+counts2[3],4*c+counts2[4])
    x24[t]<-1-x23[t]
    
    x31[t]<-rbeta(1,9*c+counts3[1],9*c+counts3[2])
    x32[t]<-1-x31[t]
    x33[t]<-rbeta(1,9*c+counts3[3],9*c+counts3[4])
    x34[t]<-1-x33[t]
    x35[t]<-rbeta(1,9*c+counts3[5],9*c+counts3[6])
    x36[t]<-1-x35[t]
    x37[t]<-rbeta(1,9*c+counts3[7],9*c+counts3[8])
    x38[t]<-1-x37[t]
    
    x41[t]<-rbeta(1,16*c+counts4[1],16*c+counts4[2])
    x42[t]<-1-x41[t]
    x43[t]<-rbeta(1,16*c+counts4[3],16*c+counts4[4])
    x44[t]<-1-x43[t]
    x45[t]<-rbeta(1,16*c+counts4[5],16*c+counts4[6])
    x46[t]<-1-x45[t]
    x47[t]<-rbeta(1,16*c+counts4[7],16*c+counts4[8])
    x48[t]<-1-x47[t]
    x49[t]<-rbeta(1,16*c+counts4[9],16*c+counts4[10])
    x410[t]<-1-x49[t]
    x411[t]<-rbeta(1,16*c+counts4[11],16*c+counts4[12])
    x412[t]<-1-x411[t]
    x413[t]<-rbeta(1,16*c+counts4[13],16*c+counts4[14])
    x414[t]<-1-x413[t]
    x415[t]<-rbeta(1,16*c+counts4[15],16*c+counts4[16])
    x416[t]<-1-x415[t]
    
    prev.mu<-rnorm(1,current.mu,tune1)
    logpriorimuc<--0.5*(100**(-2))*((current.mu-0)**2)
    logpriorimup<--0.5*(100**(-2))*((prev.mu-0)**2)
    loglikemuc<--0.5*(current.sigma**(-2))*sum(((y-current.mu)**2))+log(prod(probs.data(y,current.mu,current.sigma,x11[t],x12[t],x21[t],x22[t],x23[t],x24[t],x31[t],x32[t],x33[t],x34[t],x35[t],x36[t],x37[t],x38[t],x41[t],x42[t],x43[t],x44[t],x45[t],x46[t],x47[t],x48[t],x49[t],x410[t],x411[t],x412[t],x413[t],x414[t],x415[t],x416[t])))
    loglikemup<--0.5*(current.sigma**(-2))*sum(((y-prev.mu)**2))+log(prod(probs.data(y,prev.mu,current.sigma,x11[t],x12[t],x21[t],x22[t],x23[t],x24[t],x31[t],x32[t],x33[t],x34[t],x35[t],x36[t],x37[t],x38[t],x41[t],x42[t],x43[t],x44[t],x45[t],x46[t],x47[t],x48[t],x49[t],x410[t],x411[t],x412[t],x413[t],x414[t],x415[t],x416[t])))
    logalpha1<-(logpriorimup+loglikemup)-(logpriorimuc+loglikemuc)
    alpha1<-exp(logalpha1)
    u1<-runif(1)
    if(u1<min(1,alpha1)){
      current.mu<-prev.mu	
      accept.mu<-accept.mu+1
    }
    mu[t]<-current.mu
    
    prev.sigma<-rgamma(1,shape=current.sigma*tune2,rate=tune2)
    logpriorisc<-log(dgamma(current.sigma,shape=1,rate=0.1))
    logpriorisp<-log(dgamma(prev.sigma,shape=1,rate=0.1))
    loglikesc<--n*log(current.sigma)-0.5*(current.sigma**(-2))*sum(((y-current.mu)**2))+log(prod(probs.data(y,current.mu,current.sigma,x11[t],x12[t],x21[t],x22[t],x23[t],x24[t],x31[t],x32[t],x33[t],x34[t],x35[t],x36[t],x37[t],x38[t],x41[t],x42[t],x43[t],x44[t],x45[t],x46[t],x47[t],x48[t],x49[t],x410[t],x411[t],x412[t],x413[t],x414[t],x415[t],x416[t])))
    loglikesp<--n*log(prev.sigma)-0.5*(prev.sigma**(-2))*sum(((y-current.mu)**2))+log(prod(probs.data(y,current.mu,prev.sigma,x11[t],x12[t],x21[t],x22[t],x23[t],x24[t],x31[t],x32[t],x33[t],x34[t],x35[t],x36[t],x37[t],x38[t],x41[t],x42[t],x43[t],x44[t],x45[t],x46[t],x47[t],x48[t],x49[t],x410[t],x411[t],x412[t],x413[t],x414[t],x415[t],x416[t])))
    logajustasc<-log(dgamma(prev.sigma,shape=current.sigma*tune2,rate=tune2))
    logajustasp<-log(dgamma(current.sigma,shape=prev.sigma*tune2,rate=tune2))
    logalpha2<-(logpriorisp+loglikesp+logajustasp)-(logpriorisc+loglikesc+logajustasc)
    alpha2<-exp(logalpha2)
    u2<-runif(1)
    if(u2<min(1,alpha2)){
      current.sigma<-prev.sigma	
      accept.sigma<-accept.sigma+1
    }
    sigma[t]<-current.sigma
  }
  
  #Add burn-in and thin
  return(list(
    mu = mu[seq(burn.in+1,iterations,by = thin)],
    sigma = sigma[seq(burn.in+1,iterations,by = thin)],
    x11 = x11[seq(burn.in+1,iterations,by = thin)],
    x12 = x12[seq(burn.in+1,iterations,by = thin)],
    x21 = x21[seq(burn.in+1,iterations,by = thin)],
    x22 = x22[seq(burn.in+1,iterations,by = thin)],
    x23 = x23[seq(burn.in+1,iterations,by = thin)],
    x24 = x24[seq(burn.in+1,iterations,by = thin)],
    x31 = x31[seq(burn.in+1,iterations,by = thin)],
    x32 = x32[seq(burn.in+1,iterations,by = thin)],
    x33 = x33[seq(burn.in+1,iterations,by = thin)],
    x34 = x34[seq(burn.in+1,iterations,by = thin)],
    x35 = x35[seq(burn.in+1,iterations,by = thin)],
    x36 = x36[seq(burn.in+1,iterations,by = thin)],
    x37 = x37[seq(burn.in+1,iterations,by = thin)],
    x38 = x38[seq(burn.in+1,iterations,by = thin)],
    x41 = x41[seq(burn.in+1,iterations,by = thin)],
    x42 = x42[seq(burn.in+1,iterations,by = thin)],
    x43 = x43[seq(burn.in+1,iterations,by = thin)],
    x44 = x44[seq(burn.in+1,iterations,by = thin)],
    x45 = x45[seq(burn.in+1,iterations,by = thin)],
    x46 = x46[seq(burn.in+1,iterations,by = thin)],
    x47 = x47[seq(burn.in+1,iterations,by = thin)],
    x48 = x48[seq(burn.in+1,iterations,by = thin)],
    x49 = x49[seq(burn.in+1,iterations,by = thin)],
    x410 = x410[seq(burn.in+1,iterations,by = thin)],
    x411 = x411[seq(burn.in+1,iterations,by = thin)],
    x412 = x412[seq(burn.in+1,iterations,by = thin)],
    x413 = x413[seq(burn.in+1,iterations,by = thin)],
    x414 = x414[seq(burn.in+1,iterations,by = thin)],
    x415 = x415[seq(burn.in+1,iterations,by = thin)],
    x416 = x416[seq(burn.in+1,iterations,by = thin)],
    accept.mu = accept.mu / iterations,
    accept.sigma = accept.sigma / iterations
  ))
}

mfpt.method <- function(X, Y, Z, iterations=5000, burn.in=500, thin=10, tuneX=c(1,2), tuneY=c(5,3), tuneZ=c(100,0.2), c=1){
  
  
  out1 <- run.mfpt.single(X, iterations, burn.in, thin, tune1=tuneX[1], tune2=tuneX[2], c=1)
  mu1<-out1$mu
  sigma1<-out1$sigma
  x11_1<-out1$x11
  x12_1<-out1$x12
  x21_1<-out1$x21
  x22_1<-out1$x22
  x23_1<-out1$x23
  x24_1<-out1$x24
  x31_1<-out1$x31
  x32_1<-out1$x32
  x33_1<-out1$x33
  x34_1<-out1$x34
  x35_1<-out1$x35
  x36_1<-out1$x36
  x37_1<-out1$x37
  x38_1<-out1$x38
  x41_1<-out1$x41
  x42_1<-out1$x42
  x43_1<-out1$x43
  x44_1<-out1$x44
  x45_1<-out1$x45
  x46_1<-out1$x46
  x47_1<-out1$x47
  x48_1<-out1$x48
  x49_1<-out1$x49
  x410_1<-out1$x410
  x411_1<-out1$x411
  x412_1<-out1$x412
  x413_1<-out1$x413
  x414_1<-out1$x414
  x415_1<-out1$x415
  x416_1<-out1$x416
  accept.mu1 <- out1$accept.mu
  accept.sigma1 <- out1$accept.sigma
  
  out2 <- run.mfpt.single(Y, iterations, burn.in, thin, tune1=tuneY[1], tune2=tuneY[2], c=1)
  mu2<-out2$mu
  sigma2<-out2$sigma
  x11_2<-out2$x11
  x12_2<-out2$x12
  x21_2<-out2$x21
  x22_2<-out2$x22
  x23_2<-out2$x23
  x24_2<-out2$x24
  x31_2<-out2$x31
  x32_2<-out2$x32
  x33_2<-out2$x33
  x34_2<-out2$x34
  x35_2<-out2$x35
  x36_2<-out2$x36
  x37_2<-out2$x37
  x38_2<-out2$x38
  x41_2<-out2$x41
  x42_2<-out2$x42
  x43_2<-out2$x43
  x44_2<-out2$x44
  x45_2<-out2$x45
  x46_2<-out2$x46
  x47_2<-out2$x47
  x48_2<-out2$x48
  x49_2<-out2$x49
  x410_2<-out2$x410
  x411_2<-out2$x411
  x412_2<-out2$x412
  x413_2<-out2$x413
  x414_2<-out2$x414
  x415_2<-out2$x415
  x416_2<-out2$x416
  accept.mu2 <- out2$accept.mu
  accept.sigma2 <- out2$accept.sigma
  
  out3 <- run.mfpt.single(Z, iterations, burn.in, thin, tune1=tuneZ[1], tune2=tuneZ[2], c=1)
  mu3<-out3$mu
  sigma3<-out3$sigma
  x11_3<-out3$x11
  x12_3<-out3$x12
  x21_3<-out3$x21
  x22_3<-out3$x22
  x23_3<-out3$x23
  x24_3<-out3$x24
  x31_3<-out3$x31
  x32_3<-out3$x32
  x33_3<-out3$x33
  x34_3<-out3$x34
  x35_3<-out3$x35
  x36_3<-out3$x36
  x37_3<-out3$x37
  x38_3<-out3$x38
  x41_3<-out3$x41
  x42_3<-out3$x42
  x43_3<-out3$x43
  x44_3<-out3$x44
  x45_3<-out3$x45
  x46_3<-out3$x46
  x47_3<-out3$x47
  x48_3<-out3$x48
  x49_3<-out3$x49
  x410_3<-out3$x410
  x411_3<-out3$x411
  x412_3<-out3$x412
  x413_3<-out3$x413
  x414_3<-out3$x414
  x415_3<-out3$x415
  x416_3<-out3$x416
  accept.mu3 <- out3$accept.mu
  accept.sigma3 <- out3$accept.sigma
  
  #convert probs down branches into numeric
  probs.down.numeric<-function(x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416){
    p<-0  
    for(l in 1:16){
      p[l]<-as.numeric(probs.down(x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416)[l])  
    }  
    return(p)
  }
  
  cdf<-function(y,mu,sigma,x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416){
    n<-length(y)
    sum.down<-0
    sum.down[1]<-0
    sum<-0
    cum.distbn<-0
    for(i in 1:n){
      for(l in 2:16){  
        sum.down[l]<-probs.down.numeric(x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416)[l-1]+sum.down[l-1]
        if(label(y[i],mu,sigma)==1) sum[i]<-sum.down[1]
        else 
          if(label(y[i],mu,sigma)==l) sum[i]<-sum.down[l]
      }  
      cum.distbn[i]<-sum[i]+probs.data(y[i],mu,sigma,x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416)*(16*pnorm(y[i],mu,sigma)-label(y[i],mu,sigma)+1)
    }  
    return(cum.distbn)	
  }
  
  cdf.inv<-function(q,mu,sigma,x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416){
    N<-0
    n<-length(q)
    argument<-0
    quantile<-0
    sum.down<-0
    sum.down[1]<-probs.down.numeric(x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416)[1]
    for(i in 1:n){
      for(l in 2:16){
        sum.down[l]<-probs.down.numeric(x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416)[l]+sum.down[l-1]
        if(q[i]<=sum.down[1]) N[i]<-1
        else
          if(q[i]>sum.down[l-1] & q[i]<=sum.down[l]) N[i]<-l
      }
      argument[i]<-((q[i]-sum.down[N[i]])+(N[i]*probs.down.numeric(x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416)[N[i]]))/(16*probs.down.numeric(x11,x12,x21,x22,x23,x24,x31,x32,x33,x34,x35,x36,x37,x38,x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416)[N[i]])
      quantile[i]<-qnorm(argument[i],mu,sigma)
    }
    return(quantile)	
  }
  
  # plot the surface
  p1=seq(0.001,0.99,len=20)
  p3=seq(0.001,0.99,len=20)
  nu=length(p1)
  m=length(mu1)
  
  roc=array(0,c(nu,nu,m))
  
  for(k in 1:m){
    for(i in 1:nu){
      for(j in 1:nu){
        q1=cdf.inv(p1[i],mu1[k],sigma1[k],x11_1[k],x12_1[k],x21_1[k],x22_1[k],x23_1[k],x24_1[k],x31_1[k],x32_1[k],x33_1[k],x34_1[k],x35_1[k],x36_1[k],x37_1[k],x38_1[k],x41_1[k],x42_1[k],x43_1[k],x44_1[k],x45_1[k],x46_1[k],x47_1[k],x48_1[k],x49_1[k],x410_1[k],x411_1[k],x412_1[k],x413_1[k],x414_1[k],x415_1[k],x416_1[k])
        q3=cdf.inv(1-p3[j],mu3[k],sigma3[k],x11_3[k],x12_3[k],x21_3[k],x22_3[k],x23_3[k],x24_3[k],x31_3[k],x32_3[k],x33_3[k],x34_3[k],x35_3[k],x36_3[k],x37_3[k],x38_3[k],x41_3[k],x42_3[k],x43_3[k],x44_3[k],x45_3[k],x46_3[k],x47_3[k],x48_3[k],x49_3[k],x410_3[k],x411_3[k],x412_3[k],x413_3[k],x414_3[k],x415_3[k],x416_3[k])
        d21=cdf(q3,mu2[k],sigma2[k],x11_2[k],x12_2[k],x21_2[k],x22_2[k],x23_2[k],x24_2[k],x31_2[k],x32_2[k],x33_2[k],x34_2[k],x35_2[k],x36_2[k],x37_2[k],x38_2[k],x41_2[k],x42_2[k],x43_2[k],x44_2[k],x45_2[k],x46_2[k],x47_2[k],x48_2[k],x49_2[k],x410_2[k],x411_2[k],x412_2[k],x413_2[k],x414_2[k],x415_2[k],x416_2[k])  
        d22=cdf(q1,mu2[k],sigma2[k],x11_2[k],x12_2[k],x21_2[k],x22_2[k],x23_2[k],x24_2[k],x31_2[k],x32_2[k],x33_2[k],x34_2[k],x35_2[k],x36_2[k],x37_2[k],x38_2[k],x41_2[k],x42_2[k],x43_2[k],x44_2[k],x45_2[k],x46_2[k],x47_2[k],x48_2[k],x49_2[k],x410_2[k],x411_2[k],x412_2[k],x413_2[k],x414_2[k],x415_2[k],x416_2[k])  
        roc[i,j,k]=d21-d22
        if(roc[i,j,k]<0) roc[i,j,k]=0
      }  
    }
  }
  
  roc.m=matrix(0,nrow=nu,ncol=nu)
  for(i in 1:nu){
    for(j in 1:nu){
      roc.m[i,j]=mean(roc[i,j,])  
    }  
  }
  
  vus=numeric(m)
  for(k in 1:m){
    vus[k]=sum(sum(roc[,,k]))/(nu*nu)
  }
  
  return(list(vus=mean(vus), roc=roc.m, 
              vus.chain = vus,
              mu1 = mu1,
              mu2 = mu2, 
              mu3 = mu3,
              sigma1 = sigma1, 
              sigma2 = sigma2, 
              sigma3 = sigma3,
              accept1 = c(accept.mu1, accept.sigma1),
              accept2 = c(accept.mu2, accept.sigma2),
              accept3 = c(accept.mu3, accept.sigma3)))
}

################################################################################
# MH
################################################################################


library(MASS)
#theta=log(mu1,mu2,mu3,sigma1,sigma2,sigma3)
logprior <- function(param){
  mu <- param[1:3]
  sigma <- exp(param[4:6])
  logpriormu <- sum(dnorm(mu, mean = 0, sd = 100, log = TRUE))
  logpriorsigma <- sum(dgamma(sigma, shape = 1, rate = 0.1, log = TRUE))
  logjacobian <- sum(param[4:6])
  return(logpriormu + logpriorsigma + logjacobian)
}

loglike <- function(theta,y1,y2,y3){
  mu1 <- theta[1]
  mu2 <- theta[2]
  mu3 <- theta[3]
  sigma1 <- exp(theta[4])
  sigma2 <- exp(theta[5])
  sigma3 <- exp(theta[6])
  
  loglike1 <- sum(dnorm(y1, mean = mu1, sd = sigma1, log = TRUE))
  loglike2 <- sum(dnorm(y2, mean = mu2, sd = sigma2, log = TRUE))
  loglike3 <- sum(dnorm(y3, mean = mu3, sd = sigma3, log = TRUE))
  return(loglike1 + loglike2 + loglike3)
}

mh <- function(N, y1, y2, y3, Sigma){
  mu1 <- mean(y1)
  mu2 <- mean(y2)
  mu3 <- mean(y3)
  sigma1 <- sd(y1)
  sigma2 <- sd(y2)
  sigma3 <- sd(y3)
  init <- c(mu1, mu2, mu3, log(sigma1), log(sigma2), log(sigma3))
  
  curr <- init
  theta.mat <- matrix(0, nrow = N, ncol = 6)
  colnames(theta.mat) <- c("mu1", "mu2", "mu3", "log sigma1", "log sigma2", "log sigma3")
  theta.mat[1, ] <- curr
  count <- 0
  
  for (j in 2:N){
    can <- mvrnorm(1,curr,Sigma)
    laprob <- logprior(can) - logprior(curr) + loglike(can, y1, y2, y3) - loglike(curr, y1, y2, y3)
    if (log(runif(1)) < laprob) {
      curr <- can
      count <- count + 1
    }
    theta.mat[j, ] <- curr
  }
  
  return(list(theta = theta.mat,
              acc.prob = count/(N-1)))
}

mh.method <- function(y1, y2, y3, iterations=10000, burn.in=1000, thin=10, Sigma=diag(c(0.01, 0.01, 0.01, 0.005, 0.005, 0.005))){
  out <- mh(1000,y1,y2,y3,Sigma)
  SigmaOpt <- (2.38^2)/6*var(out$theta)
  out2 <- mh(iterations,y1,y2,y3,SigmaOpt)
  
  mu1 <- out2$theta[,1]
  mu2 <- out2$theta[,2]
  mu3 <- out2$theta[,3]
  sigma1 <- exp(out2$theta[,4])
  sigma2 <- exp(out2$theta[,5])
  sigma3 <- exp(out2$theta[,6])
  
  #Apply burn-in and thin
  mu1 = mu1[seq(burn.in+1,iterations,by = thin)]
  mu2 = mu2[seq(burn.in+1,iterations,by = thin)]
  mu3 = mu3[seq(burn.in+1,iterations,by = thin)]
  sigma1 = sigma1[seq(burn.in+1,iterations,by = thin)]
  sigma2 = sigma2[seq(burn.in+1,iterations,by = thin)]
  sigma3 = sigma3[seq(burn.in+1,iterations,by = thin)]
  
  beta1 <- sigma2/sigma1
  beta2 <- (mu1-mu2)/sigma1
  beta3 <- sigma2/sigma3
  beta4 <- (mu3-mu2)/sigma3
  
  p1=seq(0.001,0.99,len=20)
  p3=seq(0.001,0.99,len=20)
  nu=length(p1)
  m=length(mu1)
  
  roc=array(0,c(nu,nu,m))
  
  for(k in 1:m){
    for(i in 1:nu){
      for(j in 1:nu){
        q1 = (qnorm(1-p3[j]) + beta4[k]) / beta3[k]
        q2 = (qnorm(p1[i]) + beta2[k]) / beta1[k]
        roc[i,j,k] = pnorm(q1) - pnorm(q2)
      }  
    }
  }
  roc.m=matrix(0,nrow=nu,ncol=nu)
  for(i in 1:nu){
    for(j in 1:nu){
      roc.m[i,j]=mean(roc[i,j,])  
    }  
  }
  
  vus=numeric(m)
  for(k in 1:m){
    integrand <- function(s) {
      pnorm(beta1[k]*s - beta2[k]) * pnorm(-beta3[k]*s + beta4[k]) * dnorm(s)
    }
    vus[k] <- integrate(integrand, lower=-Inf, upper=Inf)$value
  }
  
  return(list(vus=mean(vus), roc=roc.m, 
              vus.chain = vus,
              out = out2$theta,
              acc.prob = out2$acc.prob
  ))
}
