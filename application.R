################################################################################
# DATASET
################################################################################

library(readxl)
pancreatic_data <- read_excel("Documents/Project/Data Set.xlsx")

# group 1 contains 40 pancreatic cancer patients
# group 7 contains 26 pancreatitis patients
# group 2 contains 40 healthy controls

ca199_data <- pancreatic_data[,c("Group", "CA199")]

X <- ca199_data[ca199_data$Group== 2,]
X <- X$CA199
X <- as.numeric(X)

Y <- ca199_data[ca199_data$Group== 7,]
Y <- Y$CA199
Y <- as.numeric(Y)

Z <- ca199_data[ca199_data$Group== 1,]
Z <- Z$CA199
Z <- as.numeric(Z)

n_1 <- length(X)
n_2 <- length(Y)
n_3 <- length(Z)
n_total <- n_1 + n_2 + n_3

################################################################################
# STATISTICS
################################################################################

function_group <-function(group) {
  group_data <- ca199_data[ca199_data$Group==group,]
  group_data$CA199 <- as.numeric(group_data$CA199)
  group_mean <- mean(group_data$CA199)
  group_sd <- sd(group_data$CA199)
  group_median <- median(group_data$CA199)
  group_min <- min(group_data$CA199)
  group_max <- max(group_data$CA199)
  group_n <- length(group_data$CA199)
  return(c(group_mean, group_sd, group_median, group_min, group_max, group_n))
}
ca199_statistics <- lapply(c(1,7,2), function_group)
ca199_statistics

################################################################################
# EMPIRICAL / MANN-WHITNEY
################################################################################

library(trinROC)
ep.vus <- emp.vus(X, Y, Z, dat = NULL, old.version = TRUE)
ep.roc <- rocsurf.emp(X, Y, Z, plot = TRUE, saveVUS = FALSE)

################################################################################
# KERNEL SMOOTHING
################################################################################

k1.vus <- k1.method(X,Y,Z)
k1.roc <- k1.roc.points(X, Y, Z)

# install.packages("rgl")
library(rgl)
p1 <- seq(0.001, 0.99, len=20)
p3 <- seq(0.001, 0.99, len=20)
trinROC.colours <- cm.colors(50) #colours
k1.colours <- trinROC.colours[cut(k1.roc[, ncol(k1.roc):1], 50, labels = FALSE)]
surface3d(p1, p3[length(p3):1], k1.roc[, ncol(k1.roc):1], 
          color = k1.colours, 
          shade = 0.75, 
          smooth = FALSE, 
          shininess = 100)
grid3d(c("x", "y", "z"), n = 10)
axes3d(c("x+", "y+", "z+"), labels = TRUE, color = "darkgray")
rgl::title3d(zlab = "TCF2", color = "darkgray", level = 2.5)
rgl::mtext3d("TCF1", "x+", color = "darkgray", level = 2.5, line = 2.5)
rgl::mtext3d("TCF3", "y+", color = "darkgray", line = 2.5, level = 2.5)

################################################################################
# MFPT
################################################################################

mfpt.out <- mfpt.method(X, Y, Z, 
                        iterations=5000, 
                        burn.in=500, 
                        thin=10, 
                        tuneX=c(1, 2), 
                        tuneY=c(5, 3), 
                        tuneZ=c(100, 0.2),
                        c=1)
mfpt.vus <- mfpt.out$vus

mfpt.colours <- trinROC.colours[cut(mfpt.out$roc[, ncol(mfpt.out$roc):1], 50, labels = FALSE)]
surface3d(p1, p3[length(p3):1], mfpt.out$roc[, ncol(mfpt.out$roc):1], 
          color = mfpt.colours, 
          shade = 0.75, 
          smooth = FALSE, 
          shininess = 100)
grid3d(c("x", "y", "z"), n = 10)
axes3d(c("x+", "y+", "z+"), labels = TRUE, color = "darkgray")
rgl::title3d(zlab = "TCF2", color = "darkgray", level = 2.5)
rgl::mtext3d("TCF1", "x+", color = "darkgray", level = 2.5, line = 2.5)
rgl::mtext3d("TCF3", "y+", color = "darkgray", line = 2.5, level = 2.5)

# diagnostic plots
n <- length(mfpt.out$mu1)
mfpt.chain <- matrix(0, nrow = n, ncol = 6)
colnames(mfpt.chain) <- c("mu1", "mu2", "mu3", "sigma1", "sigma2", "sigma3")
for (i in 1:n){
  mfpt.chain[i, ] <- c(mfpt.out$mu1[i], mfpt.out$mu2[i], mfpt.out$mu3[i], 
                       mfpt.out$sigma1[i], mfpt.out$sigma2[i], mfpt.out$sigma3[i])
}
plot(ts(mfpt.chain),main="Trace plots", col="darkorchid4")
plot(ts(mfpt.out$vus.chain),main="Trace plot", ylab="VUS", ylim=c(0.25,0.85), col="darkorchid4")
acf(mfpt.out$vus.chain, lag.max=50, main="VUS ACF")
library(coda)
effectiveSize(mfpt.out$vus.chain)

################################################################################
# MLE AND BC
################################################################################

mle.out <- trinVUS.test(X, Y, Z)
mle.vus <- mle.out$estimate
rocsurf.trin(X, Y, Z)

bc.XYZ <- boxcoxROC(X, Y, Z, verbose = FALSE)
lambda.opt <- bc.XYZ$lambda
lambda.opt
bc.out <- trinVUS.test(bc.XYZ$xbc, bc.XYZ$ybc, bc.XYZ$zbc)
bc.vus <- bc.out$estimate
rocsurf.trin(bc.XYZ$xbc, bc.XYZ$ybc, bc.XYZ$zbc)

################################################################################
# MH
################################################################################

Sigma <- diag(c(0.1, 0.5, 20, 0.005, 0.005, 0.005))
mh.out <- mh.method(X, Y, Z, iterations=10000, burn.in=1000, thin=10, Sigma)
mh.vus <- mh.out$vus

mh.colours <- trinROC.colours[cut(mh.out$roc[, ncol(mh.out$roc):1], 50, labels = FALSE)]
surface3d(p1, p3[length(p3):1], mh.out$roc[, ncol(mh.out$roc):1], 
          color = mh.colours, 
          shade = 0.75, 
          smooth = FALSE, 
          shininess = 100)
grid3d(c("x", "y", "z"), n = 10)
axes3d(c("x+", "y+", "z+"), labels = TRUE, color = "darkgray")
rgl::title3d(zlab = "TCF2", color = "darkgray", level = 2.5)
rgl::mtext3d("TCF1", "x+", color = "darkgray", level = 2.5, line = 2.5)
rgl::mtext3d("TCF3", "y+", color = "darkgray", line = 2.5, level = 2.5)

# diagnostic plots
plot(ts(mh.out$out),main="Trace plots", col="darkorchid4")
acf(mh.out$vus.chain, lag.max=50, main="VUS ACF")
library(coda)
effectiveSize(mh.out$vus.chain)

################################################################################
# LI AND ZHOU'S METHOD 1
################################################################################

m1.out <- m1.method(X, Y, Z)
m1.vus <- m1.out$vus

m1.roc=matrix(0,20,20)
for (i in 1:20){
  for (j in 1:20){
    m1.roc[i,j] <- m1.out$roc(p1[i],p3[j])
  }
}
m1.colours <- trinROC.colours[cut(m1.roc[, ncol(m1.roc):1], 50, labels = FALSE)]
surface3d(p1, p3[length(p3):1], m1.roc[, ncol(m1.roc):1], 
          color = m1.colours, 
          shade = 0.75, 
          smooth = FALSE, 
          shininess = 100)
grid3d(c("x", "y", "z"), n = 10)
axes3d(c("x+", "y+", "z+"), labels = TRUE, color = "darkgray")
rgl::title3d(zlab = "TCF2", color = "darkgray", level = 2.5)
rgl::mtext3d("TCF1", "x+", color = "darkgray", level = 2.5, line = 2.5)
rgl::mtext3d("TCF3", "y+", color = "darkgray", line = 2.5, level = 2.5)

################################################################################
# LEHMANN
################################################################################

leh.out <-leh.method(X, Y, Z)
leh.vus <- leh.out$vus

leh.roc=matrix(0,20,20)
for (i in 1:20){
  for (j in 1:20){
    leh.roc[i,j] <- leh.out$roc(p1[i],p3[j])
  }
}
leh.colours <- trinROC.colours[cut(leh.roc[, ncol(leh.roc):1], 50, labels = FALSE)]
surface3d(p1, p3[length(p3):1], leh.roc[, ncol(leh.roc):1], 
          color = leh.colours, 
          shade = 0.75, 
          smooth = FALSE, 
          shininess = 100)
grid3d(c("x", "y", "z"), n = 10)
axes3d(c("x+", "y+", "z+"), labels = TRUE, color = "darkgray")
rgl::title3d(zlab = "TCF2", color = "darkgray", level = 2.5)
rgl::mtext3d("TCF1", "x+", color = "darkgray", level = 2.5, line = 2.5)
rgl::mtext3d("TCF3", "y+", color = "darkgray", line = 2.5, level = 2.5)

################################################################################
# ASYMPTOTIC CONFIDENCE INTERVALS
################################################################################

u.out <- u.var(X,Y,Z,ep.vus)
u.ci <- u.out$ci
u.var <- u.out$sd^2
mle.ci <- mle.out$conf.int
mle.var <- mle.out$Sigma[1,1]
bc.ci <- bc.out$conf.int
bc.var <- bc.out$Sigma[1,1]
m1.model.ci <- m1.out$ci.model
m1.model.var <- m1.out$sd.model^2
m1.emp.ci <- m1.out$ci.emp
m1.emp.var <- m1.out$sd.emp^2
leh.ci <- leh.out$ci
leh.var <- leh.out$var

################################################################################
# BAYESIAN CONFIDENCE INTERVALS
################################################################################

mfpt.ci <- bayes.cred(mfpt.out$vus.chain)
mh.ci <- bayes.cred(mh.out$vus.chain)

################################################################################
# BOOTSTRAP CONFIDENCE INTERVALS
################################################################################

data.boot <- data.frame(
  marker = c(X, Y, Z),
  class = c(rep(1, length(X)), rep(2, length(Y)), rep(3, length(Z)))
)

ep.out.boot <- boot(data.boot, ep.boot, R = 1000)
ep.ci.boot <- boot.ci(ep.out.boot, type = c("perc", "bca"), conf = 0.95)
k1.out.boot <- boot(data.boot, k1.boot, R = 1000)
k1.ci.boot <- boot.ci(k1.out.boot, type = c("perc", "bca"), conf = 0.95)
mle.out.boot <- boot(data.boot, mle.boot, R = 1000)
mle.ci.boot <- boot.ci(mle.out.boot, type = c("perc", "bca", "stud"), conf = 0.95)
bc.out.boot <- boot(data.boot, bc.boot, R = 1000)
bc.ci.boot <- boot.ci(bc.out.boot, type = c("perc", "bca", "stud"), conf = 0.95)
leh.out.boot <- boot(data.boot, leh.boot, R = 1000)
leh.ci.boot <- boot.ci(leh.out.boot, type = c("perc", "bca", "stud"), conf = 0.95)

################################################################################
# RESULTS
################################################################################

vus.est <- list(ep = ep.vus, k1 = k1.vus, mfpt = mfpt.vus, mle = mle.vus, bc = bc.vus, mh = mh.vus, m1 = m1.vus, leh = leh.vus)
bayes.est <- list(mfpt.equi = mfpt.ci$equi, mfpt.hpd = mfpt.ci$hpd, mh.equi = mh.ci$equi, mh.hpd = mh.ci$hpd)
asymptotic.est <- list(u = u.ci, mle = mle.ci, bc = bc.ci, m1.model = m1.model.ci, m1.emp = m1.emp.ci, leh = leh.ci)
percent.est <- list(ep = ep.ci.boot$percent[4:5], k1 = k1.ci.boot$percent[4:5], mle = mle.ci.boot$percent[4:5], bc = bc.ci.boot$percent[4:5], leh = leh.ci.boot$percent[4:5])
bca.est <- list(ep = ep.ci.boot$bca[4:5], k1 = k1.ci.boot$bca[4:5], mle = mle.ci.boot$bca[4:5], bc = bc.ci.boot$bca[4:5], leh = leh.ci.boot$bca[4:5])
student.est <- list(mle = mle.ci.boot$student[4:5], bc = bc.ci.boot$student[4:5], leh = leh.ci.boot$student[4:5])
var.est <- list(u = u.var, mle = mle.var, bc = bc.var, m1.model = m1.model.var, m1.emp = m1.emp.var, leh = leh.var)

################################################################################
# DISTRIBUTIONAL TESTS
################################################################################

par(mfrow = c(1, 3))
hist(X, main = "Healthy", xlab = "Intensity", ylab = "Density", 
     col = "lightblue", breaks = 15)
hist(Y, main = "Pancreatitis", xlab = "Intensity", ylab = "Density", 
     col = "lightblue", breaks = 15)
hist(Z, main = "Cancer", xlab = "Intensity", ylab = "Density", 
     col = "lightblue", breaks = 15)
hist(bc.XYZ$xbc, main = "Healthy", xlab = "Intensity", ylab = "Density", 
     col = "lightblue", breaks = 15)
hist(bc.XYZ$ybc, main = "Pancreatitis", xlab = "Intensity", ylab = "Density", 
     col = "lightblue", breaks = 15)
hist(bc.XYZ$zbc, main = "Cancer", xlab = "Intensity", ylab = "Density", 
     col = "lightblue", breaks = 15)


hist(X, main = "Healthy", xlab = "Intensity", ylab = "Density", 
     col = "lightblue", breaks = 15)
hist(Y, main = "Pancreatitis", xlab = "Intensity", ylab = "Density", 
     col = "lightblue", breaks = 15)
hist(Z, main = "Cancer", xlab = "Intensity", ylab = "Density", 
     col = "lightblue", breaks = 15)

par(mfrow=c(1,2))
boxplot(X, Y, Z,
        names = c("Control", "Pancreatitis", "Cancer"),
        ylab = "CA19-9 Measurements",
        main = "Boxplot of CA19-9 Marker Measurements\nBefore Transformation",
        col = "grey80",
        border = "black",
        outpch = 16,
        outcex = 0.8,
        boxwex = 0.5,
        whisklty = 1,
        staplelty = 0,
        frame = FALSE)

boxplot(bc.XYZ$xbc, bc.XYZ$ybc, bc.XYZ$zbc,
        names = c("Control", "Pancreatitis", "Cancer"),
        ylab = "CA19-9 Measurements",
        main = "Boxplot of CA19-9 Marker Measurements\nAfter Transformation",
        col = "grey80",
        border = "black",
        outpch = 16,
        outcex = 0.8,
        boxwex = 0.5,
        whisklty = 1,
        staplelty = 0,
        frame = FALSE)

shapiro.test(X)$p.value > 0.05
shapiro.test(Y)$p.value > 0.05
shapiro.test(Z)$p.value > 0.05

plot(density(Z), ylim=c(0,0.12))
lines(density(Y))
lines(density(X))

plot(density(bc.XYZ$zbc), ylim=c(0,0.6))
lines(density(bc.XYZ$ybc))
lines(density(bc.XYZ$xbc))

