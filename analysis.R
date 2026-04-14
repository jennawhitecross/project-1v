################################################################################
# SCENARIO 1
################################################################################
# Bayesian cover
################################################################################

NP.BA.1 <- c(0.94,0.88)
P.BA.1 <- c(0.30,0.12)

NP.BA.2 <- c(0.94,0.88)
P.BA.2 <- c(0.20,0.10)

NP.BA.3 <- c(0.76,0.68)
P.BA.3 <- c(0.04,0.00)

NP.BA.4 <- c(0.50,0.30)
P.BA.4 <- c(0.00,0.00)

NP.BA.5 <- c(0.32,0.26)
P.BA.5 <- c(0.00,0.00)

NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))
P.BA <- c(mean(P.BA.1), mean(P.BA.2), mean(P.BA.3), mean(P.BA.4), mean(P.BA.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP.BA, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Bayesian Mean Coverage Probability", cex.main = 1)

lines(x.axis, P.BA, type = "b", col = "red", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("topright", 
       legend = c("MFPT", "MH", "Nominal Level (0.95)"),
       col = c("blue", "red", "black"),
       lty = c(1, 1, 2), 
       pch = c(20, 20, NA), 
       lwd = 2,
       cex = 0.9)

################################################################################
# Scenario 1 - Bayesian alen
################################################################################

NP.BA.1 <- c(0.078,0.064)
P.BA.1 <- c(0.016,0.012)

NP.BA.2 <- c(0.063,0.051)
P.BA.2 <- c(0.015,0.011)

NP.BA.3 <- c(0.055,0.043)
P.BA.3 <- c(0.010,0.007)

NP.BA.4 <- c(0.026,0.020)
P.BA.4 <- c(0.005,0.004)

NP.BA.5 <- c(0.020,0.016)
P.BA.5 <- c(0.004,0.003)

NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))
P.BA <- c(mean(P.BA.1), mean(P.BA.2), mean(P.BA.3), mean(P.BA.4), mean(P.BA.5))

x.axis <- 1:5

plot(x.axis, NP.BA, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.1),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Bayesian Mean Average Length", cex.main = 1)

lines(x.axis, P.BA, type = "b", col = "red", lwd = 2, pch = 20)

legend("topright", 
       legend = c("MFPT", "MH"),
       col = c("blue", "red"),
       lty = c(1, 1), 
       pch = c(20, 20), 
       lwd = 2,
       cex = 1)

################################################################################
# Scenario 1 - Bayesian clen
################################################################################

NP.BA.1 <- c(0.083,0.071)
P.BA.1 <- c(0.030,0.029)

NP.BA.2 <- c(0.067,0.056)
P.BA.2 <- c(0.031, 0.031)

NP.BA.3 <- c(0.067,0.056)

NP.BA.4 <- c(0.039,0.038)

NP.BA.5 <- c(0.037,0.031)

NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))
P.BA <- c(mean(P.BA.1), mean(P.BA.2))

x.axis <- 1:5

plot(x.axis, NP.BA, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.1),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Bayesian Mean Conditional Average Length", cex.main = 1)

lines(1:2, P.BA, type = "b", col = "red", lwd = 2, pch = 20)

legend("topright", 
       legend = c("MFPT", "MH"),
       col = c("blue", "red"),
       lty = c(1, 1), 
       pch = c(20, 20), 
       lwd = 2,
       cex = 1)

################################################################################
# SCENARIO 2
################################################################################
# cover
################################################################################

NP.1 <- c(0.92,0.94,0.94,0.88,0.86,0.90,0.86)
P.1 <- c(0.92,0.92,0.92,0.94,0.96,0.92)
SP.1 <- c(0.94,0.78,0.90,0.88,0.76)

NP.2 <- c(0.90,0.94,0.94,0.94,0.88,0.90,0.86)
P.2 <- c(0.92,0.94,0.94,0.94,0.94,0.90)
SP.2 <- c(0.84,0.82,0.91,0.84,0.84)

NP.3 <- c(0.96,0.96,0.88,0.84,1.00,1.00)
P.3 <- c(0.92,0.90,0.92,0.94,0.92,0.90)
SP.3 <- c(0.86,0.72,0.94,0.84,0.82)

NP.4 <- c(0.92,0.92,0.80,0.74,0.92,0.92)
P.4 <- c(0.92,0.94,0.94,0.94,0.96,0.96)
SP.4 <- c(0.80,0.58,0.80,0.68,0.56)

NP.5 <- c(0.90,0.88,0.80,0.82,0.82,0.80)
P.5 <- c(0.88,0.88,0.90,0.90,0.94,0.92)
SP.5 <- c(0.82,0.46,0.78,0.64,0.56)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.45,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Mean Coverage Probability\nNonparametric vs Parametric vs Semi-Parametric", cex.main=1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "Semi-Parametric", "Nominal Level (0.95)"),
       col = c("blue", "red", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.9)

################################################################################
#  Scenario 2 -alen
################################################################################

NP.1 <- c(0.286,0.324,0.324,0.283,0.283,0.247,0.240)
P.1 <- c(0.288,0.293,0.293,0.306,0.287,0.283)
SP.1 <- c(0.262,0.252,0.298,0.289,0.295)

NP.2 <- c(0.264,0.297,0.297,0.261,0.260,0.243,0.236)
P.2 <- c(0.268,0.269,0.269,0.279,0.272,0.268)
SP.2 <- c(0.244,0.230,0.268,0.263,0.265)

NP.3 <- c(0.266,0.266,0.237,0.235,0.211,0.204)
P.3 <- c(0.241,0.245,0.247,0.253,0.238,0.235)
SP.3 <- c(0.206,0.198,0.247,0.248,0.246)

NP.4 <- c(0.218,0.217,0.199,0.198,0.171,0.165)
P.4 <- c(0.203,0.205,0.207,0.210,0.203,0.200)
SP.4 <- c(0.163,0.173,0.206,0.205,0.205)

NP.5 <- c(0.197,0.197,0.180,0.179,0.157,0.153)
P.5 <- c(0.183,0.187,0.187,0.190,0.182,0.179)
SP.5 <- c(0.158,0.156,0.183,0.181,0.180)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.1,0.3),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Mean Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "Semi-Parametric"),
       col = c("blue", "red", "green"),
       lty = c(1, 1), 
       pch = c(20, 20), 
       lwd = 2,
       cex = 1)

################################################################################
# Scenario 2 - clen
################################################################################

NP.1 <- c(0.288,0.325,0.326,0.286,0.286,0.249,0.243)
P.1 <- c(0.289,0.296,0.296,0.308,0.288,0.285)
SP.1 <- c(0.262,0.252,0.298,0.289,0.295)

NP.2 <- c(0.266,0.299,0.299,0.263,0.262,0.247,0.240)
P.2 <- c(0.270,0.271,0.270,0.280,0.271,0.269)
SP.2 <- c(0.256,0.230,0.268,0.263,0.265)

NP.3 <- c(0.267,0.267,0.240,0.239,0.211,0.204)
P.3 <- c(0.242,0.246,0.247,0.254,0.239,0.236)
SP.3 <- c(0.208,0.198,0.247,0.248,0.246)

NP.4 <- c(0.219,0.219,0.200,0.199,0.172,0.166)
P.4 <- c(0.204,0.206,0.207,0.210,0.204,0.201)
SP.4 <- c(0.165,0.173,0.206,0.205,0.205)

NP.5 <- c(0.198,0.198,0.182,0.181,0.158,0.153)
P.5 <- c(0.183,0.188,0.188,0.191,0.182,0.180)
SP.5 <- c(0.159,0.156,0.183,0.181,0.180)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.1,0.3),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Mean Conditional Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "Semi-Parametric"),
       col = c("blue", "red", "green"),
       lty = c(1, 1), 
       pch = c(20, 20), 
       lwd = 2,
       cex = 1)

################################################################################
# Scenario 2 - Parametric cover
################################################################################

P.A.1 <- c(0.92)
P.B.1 <- c(0.92,0.92,0.94)
P.BA.1 <- c(0.96,0.92)

P.A.2 <- c(0.92)
P.B.2 <- c(0.94,0.94,0.94)
P.BA.2 <- c(0.94,0.90)

P.A.3 <- c(0.92)
P.B.3 <- c(0.90,0.92,0.94)
P.BA.3 <- c(0.92,0.90)

P.A.4 <- c(0.92)
P.B.4 <- c(0.94,0.94,0.94)
P.BA.4 <- c(0.96,0.96)

P.A.5 <- c(0.88)
P.B.5 <- c(0.88,0.90,0.90)
P.BA.5 <- c(0.94,0.92)

P.A <- c(mean(P.A.1), mean(P.A.2), mean(P.A.3), mean(P.A.4), mean(P.A.5))
P.B <- c(mean(P.B.1), mean(P.B.2), mean(P.B.3), mean(P.B.4), mean(P.B.5))
P.BA <- c(mean(P.BA.1), mean(P.BA.2), mean(P.BA.3), mean(P.BA.4), mean(P.BA.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, P.A, type = "b", col = "red", lwd = 2, pch = 20,
     ylim = c(0.83,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Parametric Mean Coverage Probability", cex.main = 1)

lines(x.axis, P.B, type = "b", col = "blue", lwd = 2, pch = 20)
lines(x.axis, P.BA, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian", "Nominal Level (0.95)"),
       col = c("red", "blue", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.75)

################################################################################
# Scenario 2 - Parametric alen
################################################################################

P.A.1 <- c(0.288)
P.B.1 <- c(0.293,0.293,0.306)
P.BA.1 <- c(0.287,0.283)

P.A.2 <- c(0.268)
P.B.2 <- c(0.269,0.269,0.279)
P.BA.2 <- c(0.272,0.268)

P.A.3 <- c(0.241)
P.B.3 <- c(0.245,0.247,0.253)
P.BA.3 <- c(0.238,0.235)

P.A.4 <- c(0.203)
P.B.4 <- c(0.205,0.207,0.210)
P.BA.4 <- c(0.203,0.200)

P.A.5 <- c(0.183)
P.B.5 <- c(0.187,0.187,0.190)
P.BA.5 <- c(0.182,0.179)

P.A <- c(mean(P.A.1), mean(P.A.2), mean(P.A.3), mean(P.A.4), mean(P.A.5))
P.B <- c(mean(P.B.1), mean(P.B.2), mean(P.B.3), mean(P.B.4), mean(P.B.5))
P.BA <- c(mean(P.BA.1), mean(P.BA.2), mean(P.BA.3), mean(P.BA.4), mean(P.BA.5))

x.axis <- 1:5

plot(x.axis, P.A, type = "b", col = "red", lwd = 2, pch = 20,
     ylim = c(0.1,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Parametric Mean Average Length", cex.main = 1)

lines(x.axis, P.B, type = "b", col = "blue", lwd = 2, pch = 20)
lines(x.axis, P.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.75)

################################################################################
# Scenario 2 - Parametric clen
################################################################################

P.A.1 <- c(0.289)
P.B.1 <- c(0.296,0.296,0.308)
P.BA.1 <- c(0.288,0.285)

P.A.2 <- c(0.270)
P.B.2 <- c(0.271,0.270,0.280)
P.BA.2 <- c(0.273,0.269)

P.A.3 <- c(0.242)
P.B.3 <- c(0.246,0.247,0.254)
P.BA.3 <- c(0.239,0.236)

P.A.4 <- c(0.204)
P.B.4 <- c(0.206,0.207,0.210)
P.BA.4 <- c(0.204,0.201)

P.A.5 <- c(0.183)
P.B.5 <- c(0.188,0.188,0.191)
P.BA.5 <- c(0.182,0.180)

P.A <- c(mean(P.A.1), mean(P.A.2), mean(P.A.3), mean(P.A.4), mean(P.A.5))
P.B <- c(mean(P.B.1), mean(P.B.2), mean(P.B.3), mean(P.B.4), mean(P.B.5))
P.BA <- c(mean(P.BA.1), mean(P.BA.2), mean(P.BA.3), mean(P.BA.4), mean(P.BA.5))

x.axis <- 1:5

plot(x.axis, P.A, type = "b", col = "red", lwd = 2, pch = 20,
     ylim = c(0.1,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Parametric Mean Conditional Average Length", cex.main = 1)

lines(x.axis, P.B, type = "b", col = "blue", lwd = 2, pch = 20)
lines(x.axis, P.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.75)

################################################################################
# SCENARIO 3
################################################################################
# cover
################################################################################

NP.1 <- c(0.92,0.92,0.98,0.96,0.98,0.96,0.96)
P.1 <- c(0.92,0.92,0.94,0.98,0.96,0.88)
SP.1 <- c(0.56,0.92,0.94,0.94,0.94)

NP.2 <- c(0.96,0.98,0.96,0.96,0.96,0.96,0.96)
P.2 <- c(0.94,0.96,0.96,0.96,1.00,1.00)
SP.2 <- c(0.52,0.92,0.94,0.94,0.98)

NP.3 <- c(0.92,0.96,0.92,0.94,0.90,0.92)
P.3 <- c(0.90,0.90,0.92,0.94,0.88,0.88)
SP.3 <- c(0.50,0.88,0.90,0.96,0.98)

NP.4 <- c(0.98,1.00,0.98,1.00,0.94,0.90)
P.4 <- c(0.94,0.96,0.98,1.00,0.96,0.96)
SP.4 <- c(0.64,0.98,0.98,0.98,0.98)

NP.5 <- c(0.96,0.94,0.96,0.98,0.90,0.88)
P.5 <- c(0.94,0.98,0.98,0.98,0.98,0.96)
SP.5 <- c(0.62,0.94,0.94,0.94,0.96)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
M1 <- c(mean(M1.1), mean(M1.2), mean(M1.3), mean(M1.4), mean(M1.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.75,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Mean Coverage Probability\nNonparametric vs Parametric vs Semi-Parametric", cex.main=1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "Semi-Parametric", "Nominal Level (0.95)"),
       col = c("blue", "red", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.8)

################################################################################
# Scenario 3 - alen
################################################################################

NP.1 <- c(0.180,0.206,0.215,0.177,0.181,0.161,0.154)
P.1 <- c(0.183,0.189,0.196,0.213,0.181,0.176)
SP.1 <- c(0.120,0.159,0.169,0.173,0.180)

NP.2 <- c(0.177,0.200,0.206,0.174,0.178,0.151,0.145)
P.2 <- c(0.174,0.181,0.185,0.197,0.173,0.167)
SP.2 <- c(0.111,0.147,0.159,0.161,0.163)

NP.3 <- c(0.174,0.179,0.154,0.158,0.137,0.132)
P.3 <- c(0.157,0.163,0.167,0.176,0.156,0.151)
SP.3 <- c(0.111,0.130,0.137,0.139,0.137)

NP.4 <- c(0.144,0.147,0.130,0.132,0.117,0.112)
P.4 <- c(0.133,0.136,0.138,0.144,0.133,0.131)
SP.4 <- c(0.082,0.112,0.117,0.118,0.118)

NP.5 <- c(0.128,0.130,0.116,0.118,0.105,0.102)
P.5 <- c(0.119,0.121,0.123,0.125,0.118,0.116)
SP.5 <- c(0.080,0.100,0.105,0.106,0.106)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.2),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Mean Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "Semi-Parametric"),
       col = c("blue", "red", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 1)

################################################################################
# Scenario 3 - clen
################################################################################

NP.1 <- c(0.180,0.206,0.215,0.177,0.181,0.161,0.154)
P.1 <- c(0.183,0.189,0.196,0.213,0.181,0.176)
SP.1 <- c(0.133,0.159,0.169,0.173,0.180)

NP.2 <- c(0.177,0.200,0.206,0.174,0.178,0.151,0.145)
P.2 <- c(0.174,0.181,0.185,0.197,0.173,0.167)
SP.2 <- c(0.130,0.147,0.159,0.161,0.163)

NP.3 <- c(0.174,0.179,0.154,0.158,0.137,0.132)
P.3 <- c(0.157,0.163,0.167,0.176,0.156,0.151)
SP.3 <- c(0.115,0.130,0.137,0.139,0.137)

NP.4 <- c(0.144,0.147,0.130,0.132,0.117,0.112)
P.4 <- c(0.133,0.136,0.138,0.144,0.133,0.131)
SP.4 <- c(0.085,0.112,0.117,0.118,0.118)

NP.5 <- c(0.128,0.130,0.116,0.118,0.105,0.102)
P.5 <- c(0.119,0.121,0.123,0.125,0.118,0.116)
SP.5 <- c(0.082,0.100,0.105,0.106,0.106)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.2),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Mean Conditional Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "Semi-Parametric"),
       col = c("blue", "red", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 1)

################################################################################
# Scenario 3 - Nonparametric cover
################################################################################

NP.A.1 <- c(0.92)
NP.B.1 <- c(0.92,0.98,0.96,0.98)
NP.BA.1 <- c(0.96,0.96)

NP.A.2 <- c(0.96)
NP.B.2 <- c(0.98,0.96,0.96,0.96)
NP.BA.2 <- c(0.96,0.96)

NP.B.3 <- c(0.92,0.96, 0.92,0.94)
NP.BA.3 <- c(0.90,0.92)

NP.B.4 <- c(0.98,1.00,0.98,1.00)
NP.BA.4 <- c(0.94,0.90)

NP.B.5 <- c(0.96,0.94,0.96,0.98)
NP.BA.5 <- c(0.90,0.88)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))
LEH <- c(mean(LEH.1), mean(LEH.2), mean(LEH.3), mean(LEH.4), mean(LEH.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.83,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Nonparametric Mean Coverage Probability", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian", "Nominal Level (0.95)"),
       col = c("red", "blue", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 3 - Nonparametric alen
################################################################################

NP.A.1 <- c(0.180)
NP.B.1 <- c(0.206,0.215,0.177,0.181)
NP.BA.1 <- c(0.161,0.154)

NP.A.2 <- c(0.177)
NP.B.2 <- c(0.200,0.206,0.174,0.178)
NP.BA.2 <- c(0.151,0.145)

NP.B.3 <- c(0.174,0.179,0.154,0.158)
NP.BA.3 <- c(0.137,0.132)

NP.B.4 <- c(0.144,0.147,0.130,0.132)
NP.BA.4 <- c(0.117,0.112)

NP.B.5 <- c(0.128,0.130,0.116,0.118)
NP.BA.5 <- c(0.105,0.102)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.05,0.2),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Nonparametric Mean Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 3 - Nonparametric clen
################################################################################

NP.A.1 <- c(0.182)
NP.B.1 <- c(0.209,0.214,0.179,0.180)
NP.BA.1 <- c(0.159,0.152)

NP.A.2 <- c(0.176)
NP.B.2 <- c(0.198,0.203,0.171,0.175)
NP.BA.2 <- c(0.150,0.144)

NP.B.3 <- c(0.175,0.180,0.155,0.157)
NP.BA.3 <- c(0.134,0.129)

NP.B.4 <- c(0.145,0.147,0.130,0.132)
NP.BA.4 <- c(0.117,0.112)

NP.B.5 <- c(0.128,0.128,0.116,0.117)
NP.BA.5 <- c(0.104,0.100)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.05,0.2),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Nonparametric Mean Conditional Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# SCENARIO 4
################################################################################
# cover
################################################################################

NP.1 <- c(0.88,0.94,0.94,0.82,0.68,0.70,0.68)
P.1 <- c(0.78,0.90,0.84,0.88,0.72,0.70)
BC.1 <- c(0.96,0.96,0.98,1.00)
SP.1 <- c(0.58,0.82,0.90,0.86,0.84)

NP.2 <- c(0.86,0.92,0.92,0.90,0.86,0.72,0.70)
P.2 <- c(0.74,0.84,0.72,0.74,0.76,0.74)
BC.2 <- c(0.92,0.92,0.92,0.96)
SP.2 <- c(0.58,0.74,0.94,0.92,0.82)

NP.3 <- c(0.98,0.98,0.82,0.72,0.68,0.70)
P.3 <- c(0.70,0.78,0.54,0.58,0.64,0.60)
BC.3 <- c(0.94,0.92,0.92,0.96)
SP.3 <- c(0.64,0.72,0.92,0.92,0.80)

NP.4 <- c(0.98,0.98,0.72,0.72,0.58,0.60)
P.4 <- c(0.56,0.52,0.40,0.48,0.32,0.32)
BC.4 <- c(0.96,0.96,0.96,0.96)
SP.4 <- c(0.50,0.64,0.84,0.78,0.68)

NP.5 <- c(0.90,0.96,0.74,0.72,0.64,0.62)
P.5 <- c(0.52,0.50,0.28,0.30,0.40,0.40)
BC.5 <- c(0.90,0.92,0.94,0.92)
SP.5 <- c(0.54,0.66,0.78,0.66,0.62)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Mean Coverage Probability\nNonparametric vs Parametric vs Semi-Parametric", cex.main=1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, BC, type = "b", col = "darkorange", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "With Box-Cox Transformation", "Semi-Parametric", "Nominal Level (0.95)"),
       col = c("blue", "red", "darkorange", "green", "black"),
       lty = c(1, 1, 1, 1, 2), 
       pch = c(20, 20, 20, 20, NA), 
       lwd = 2,
       cex = 0.80)

################################################################################
# Scenario 4 - alen
################################################################################

NP.1 <- c(0.285,0.321,0.322,0.284,0.287,0.243,0.236)
P.1 <- c(0.292,0.297,0.301,0.307,0.281,0.277)
BC.1 <- c(0.290,0.303,0.304,0.318)
SP.1 <- c(0.124,0.252,0.300,0.292,0.301)

NP.2 <- c(0.268,0.299,0.299,0.266,0.266,0.227,0.218)
P.2 <- c(0.266,0.264,0.262,0.269,0.263,0.259)
BC.2 <- c(0.271,0.278,0.279,0.289)
SP.2 <- c(0.140,0.231,0.271,0.267,0.268)

NP.3 <- c(0.261,0.261,0.228,0.229,0.189,0.183)
P.3 <- c(0.222,0.211,0.211,0.213,0.215,0.212)
BC.3 <- c(0.238,0.246,0.247,0.254)
SP.3 <- c(0.101,0.199,0.245,0.244,0.244)

NP.4 <- c(0.218,0.218,0.188,0.187,0.162,0.157)
P.4 <- c(0.191,0.186,0.190,0.187,0.188,0.186)
BC.4 <- c(0.204,0.209,0.209,0.214)
SP.4 <- c(0.086,0.174,0.202,0.200,0.200)

NP.5 <- c(0.197,0.198,0.175,0.175,0.150,0.146)
P.5 <- c(0.186,0.184,0.186,0.186,0.182,0.180)
BC.5 <- c(0.183,0.185,0.184,0.188)
SP.5 <- c(0.081,0.157,0.182,0.180,0.180)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.3),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Mean Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, BC, type = "b", col = "darkorange", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "With Box-Cox Transformation", "Semi-Parametric"),
       col = c("blue", "red", "darkorange", "green"),
       lty = c(1, 1, 1, 1), 
       pch = c(20, 20, 20, 20), 
       lwd = 2,
       cex = 0.8)

################################################################################
# Scenario 4 - clen
################################################################################

NP.1 <- c(0.287,0.323,0.323,0.286,0.293,0.249,0.242)
P.1 <- c(0.301,0.298,0.299,0.300,0.296,0.291)
BC.1 <- c(0.291,0.304,0.305,0.318)
SP.1 <- c(0.144,0.257,0.306,0.299,0.309)

NP.2 <- c(0.271,0.301,0.302,0.269,0.269,0.230,0.221)
P.2 <- c(0.276,0.271,0.254,0.248,0.270,0.267)
BC.2 <- c(0.273,0.281,0.283,0.291)
SP.2 <- c(0.173,0.236,0.275,0.271,0.268)

NP.3 <- c(0.262,0.262,0.231,0.233, 0.191,0.184)
P.3 <- c(0.226,0.213,0.212,0.217,0.220,0.219)
BC.3 <- c(0.239,0.249,0.250,0.256)
SP.3 <- c(0.114,0.201,0.247,0.246,0.241)

NP.4 <- c(0.218,0.218,0.191,0.190, 0.163,0.160)
P.4 <- c(0.197,0.190,0.182,0.186,0.195,0.193)
BC.4 <- c(0.205,0.210,0.210,0.215)
SP.4 <- c(0.108,0.176,0.204,0.203,0.203)

NP.5 <- c(0.197,0.198,0.178,0.177,0.152,0.148)
P.5 <- c(0.192,0.187,0.198,0.190,0.192,0.190)
BC.5 <- c(0.183,0.185,0.184,0.188)
SP.5 <- c(0.095,0.159,0.184,0.181,0.182)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.3),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Mean Conditional Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, BC, type = "b", col = "darkorange", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "With Box-Cox Transformation", "Semi-Parametric"),
       col = c("blue", "red", "darkorange", "green"),
       lty = c(1, 1, 1, 1), 
       pch = c(20, 20, 20, 20), 
       lwd = 2,
       cex = 0.8)

################################################################################
# Scenario 4 - Nonparametric cover
################################################################################

NP.A.1 <- c(0.88)
NP.B.1 <- c(0.94,0.94,0.82,0.68)
NP.BA.1 <- c(0.70,0.68)

NP.A.2 <- c(0.86)
NP.B.2 <- c(0.92,0.92,0.90,0.86)
NP.BA.2 <- c(0.72,0.70)

NP.B.3 <- c(0.98,0.98,0.82,0.72)
NP.BA.3 <- c(0.68,0.70)

NP.B.4 <- c(0.98,0.98,0.72,0.72)
NP.BA.4 <- c(0.58,0.60)

NP.B.5 <- c(0.90,0.96,0.74,0.72)
NP.BA.5 <- c(0.64,0.62)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.45,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Nonparametric Mean Coverage Probability", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian", "Nominal Level (0.95)"),
       col = c("red", "blue", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 4 - Nonparametric alen
################################################################################

NP.A.1 <- c(0.285)
NP.B.1 <- c(0.321,0.322,0.284,0.287)
NP.BA.1 <- c(0.243,0.236)

NP.A.2 <- c(0.268)
NP.B.2 <- c(0.299,0.299,0.266,0.266)
NP.BA.2 <- c(0.227,0.218)

NP.B.3 <- c(0.261,0.261,0.228,0.229)
NP.BA.3 <- c(0.189,0.183)

NP.B.4 <- c(0.218,0.218,0.188,0.187)
NP.BA.4 <- c(0.162,0.157)

NP.B.5 <- c(0.197,0.198,0.175,0.175)
NP.BA.5 <- c(0.150,0.146)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Nonparametric Mean Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 4 - Nonparametric clen
################################################################################

NP.A.1 <- c(0.287)
NP.B.1 <- c(0.323,0.323,0.286,0.293)
NP.BA.1 <- c(0.249,0.242)

NP.A.2 <- c(0.271)
NP.B.2 <- c(0.301,0.302,0.269,0.269)
NP.BA.2 <- c(0.230,0.221)

NP.B.3 <- c(0.262,0.262,0.231,0.233)
NP.BA.3 <- c(0.191,0.184)

NP.B.4 <- c(0.218,0.218,0.191,0.190)
NP.BA.4 <- c(0.163,0.160)

NP.B.5 <- c(0.197,0.198,0.178,0.177)
NP.BA.5 <- c(0.152,0.148)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Nonparametric Mean Conditional Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# SCENARIO 5
################################################################################
# cover
################################################################################

NP.1 <- c(0.96,0.98,0.98,0.90,0.90,0.78,0.74)
P.1 <- c(0.86,0.96,0.90,0.96,0.78,0.78)
BC.1 <- c(0.90,0.90,0.90,0.96)
SP.1 <- c(0.70,0.68,0.84,0.78,0.70)

NP.2 <- c(0.90,0.92,0.92,0.88,0.84,0.86,0.82)
P.2 <- c(0.86,0.92,0.88,0.90,0.84,0.84)
BC.2 <- c(0.84,0.84,0.88,0.94)
SP.2 <- c(0.80,0.72,0.84,0.86,0.82)

NP.3 <- c(0.94,0.96,0.84,0.80,0.84,0.82)
P.3 <- c(0.72,0.88,0.86,0.88,0.82,0.80)
BC.3 <- c(0.92,0.92,0.92,0.92)
SP.3 <- c(0.81,0.58,0.84,0.76,0.72)

NP.4 <- c(0.96,0.96,0.86,0.84,0.84,0.82)
P.4 <- c(0.84,0.92,0.88,0.90,0.82,0.80)
BC.4 <- c(0.94,0.92,0.94,0.94)
SP.4 <- c(0.72,0.60,0.76,0.72,0.66)

NP.5 <- c(1.00,1.00,0.90,0.90,0.84,0.80)
P.5 <- c(0.90,0.94,0.90,0.96,0.84,0.86)
BC.5 <- c(0.96,0.96,0.98,0.98)
SP.5 <- c(0.84,0.74,0.66,0.58)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.5,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Mean Coverage Probability\nNonparametric vs Parametric vs Semi-Parametric", cex.main=1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, BC, type = "b", col = "darkorange", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "With Box-Cox Transformation", "Semi-Parametric", "Nominal Level (0.95)"),
       col = c("blue", "red", "darkorange", "green", "black"),
       lty = c(1, 1, 1, 1, 2), 
       pch = c(20, 20, 20, 20, NA), 
       lwd = 2,
       cex = 0.80)

################################################################################
# Scenario 5 - alen
################################################################################

NP.1 <- c(0.281,0.316,0.316,0.277,0.275,0.233,0.226)
P.1 <- c(0.279,0.296,0.294,0.309,0.269,0.266)
BC.1 <- c(0.286,0.296,0.296,0.307)
SP.1 <- c(0.207,0.218,0.222,0.243,0.276,0.271,0.270)

NP.2 <- c(0.262,0.292,0.292,0.265,0.264,0.226,0.219)
P.2 <- c(0.258,0.276,0.275,0.283,0.255,0.252)
BC.2 <- c(0.271,0.277,0.275,0.284)
SP.2 <- c(0.210,0.218,0.269,0.264,0.263)

NP.3 <- c(0.258,0.259,0.229,0.228,0.163,0.160)
P.3 <- c(0.218,0.243,0.241,0.249,0.191,0.188)
BC.3 <- c(0.238,0.245,0.244,0.251)
SP.3 <- c(0.182,0.187,0.237,0.232,0.233)

NP.4 <- c(0.218,0.218,0.197,0.196,0.163,0.160)
P.4 <- c(0.192,0.209,0.210,0.214,0.191,0.188)
BC.4 <- c(0.202,0.204,0.204,0.207)
SP.4 <- c(0.142,0.167,0.200,0.197,0.197)

NP.5 <- c(0.195,0.194,0.179,0.179,0.153,0.148)
P.5 <- c(0.180,0.189,0.189,0.192,0.177,0.175)
BC.5 <- c(0.181,0.186,0.185,0.188)
SP.5 <- c(0.133,0.152,0.173,0.170,0.169)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.3),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Mean Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, BC, type = "b", col = "darkorange", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "With Box-Cox Transformation", "Semi-Parametric"),
       col = c("blue", "red", "darkorange", "green"),
       lty = c(1, 1, 1, 1), 
       pch = c(20, 20, 20, 20), 
       lwd = 2,
       cex = 0.8)

################################################################################
# Scenario 5 - clen
################################################################################

NP.1 <- c(0.281,0.316,0.316,0.279,0.276,0.238,0.233)
P.1 <- c(0.284,0.296,0.293,0.308,0.276,0.272)
BC.1 <- c(0.286,0.299,0.298,0.308)
SP.1 <- c(0.223,0.254,0.280,0.275,0.275)

NP.2 <- c(0.265,0.294,0.294,0.268,0.269,0.230,0.223)
P.2 <- c(0.262,0.279,0.278,0.285,0.261,0.257)
BC.2 <- c(0.270,0.278,0.276,0.285)
SP.2 <- c(0.227,0.224,0.273,0.268,0.267)

NP.3 <- c(0.259,0.260,0.233,0.233,0.166,0.162)
P.3 <- c(0.222,0.247,0.245,0.251,0.194,0.191)
BC.3 <- c(0.240,0.249,0.248,0.256)
SP.3 <- c(0.196,0.194,0.242,0.237,0.240)

NP.4 <- c(0.218,0.218,0.199,0.198,0.166,0.162)
P.4 <- c(0.194,0.210,0.211,0.215,0.194,0.191)
BC.4 <- c(0.203,0.205,0.204,0.208)
SP.4 <- c(0.150,0.173,0.206,0.204,0.201)

NP.5 <- c(0.195,0.194,0.181,0.181,0.155,0.150)
P.5 <- c(0.181,0.189,0.190,0.192,0.179,0.177)
BC.5 <- c(0.181,0.186,0.185,0.188)
SP.5 <- c(0.147,0.157,0.176,0.173,0.172)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.3),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Mean Conditional Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, BC, type = "b", col = "darkorange", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "With Box-Cox Transformation", "Semi-Parametric"),
       col = c("blue", "red", "darkorange", "green"),
       lty = c(1, 1, 1, 1), 
       pch = c(20, 20, 20, 20), 
       lwd = 2,
       cex = 0.8)

################################################################################
# Scenario 5 - Nonparametric cover
################################################################################

NP.A.1 <- c(0.90)
NP.B.1 <- c(0.98,0.98,0.90,0.90)
NP.BA.1 <- c(0.78,0.74)

NP.A.2 <- c(0.94)
NP.B.2 <- c(0.92,0.92,0.88,0.84)
NP.BA.2 <- c(0.86,0.82)

NP.B.3 <- c(0.94,0.96,0.84,0.80)
NP.BA.3 <- c(0.84,0.82)

NP.B.4 <- c(0.96,0.96,0.86,0.84)
NP.BA.4 <- c(0.84,0.82)

NP.B.5 <- c(1.00,1.00,0.90,0.90)
NP.BA.5 <- c(0.84,0.80)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.45,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Nonparametric Mean Coverage Probability", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian", "Nominal Level (0.95)"),
       col = c("red", "blue", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 5 - Nonparametric alen
################################################################################

NP.A.1 <- c(0.281)
NP.B.1 <- c(0.316,0.316,0.277,0.275)
NP.BA.1 <- c(0.233,0.226)

NP.A.2 <- c(0.262)
NP.B.2 <- c(0.292,0.292,0.265,0.264)
NP.BA.2 <- c(0.226,0.219)

NP.B.3 <- c(0.258,0.259,0.229,0.228)
NP.BA.3 <- c(0.163,0.160)

NP.B.4 <- c(0.218,0.218,0.197,0.196)
NP.BA.4 <- c(0.163,0.160)

NP.B.5 <- c(0.195,0.194,0.179,0.179)
NP.BA.5 <- c(0.153,0.148)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Nonparametric Mean Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 5 - Nonparametric clen
################################################################################

NP.A.1 <- c(0.281)
NP.B.1 <- c(0.316,0.316,0.279,0.276)
NP.BA.1 <- c(0.238,0.233)

NP.A.2 <- c(0.265)
NP.B.2 <- c(0.294,0.294,0.268,0.269)
NP.BA.2 <- c(0.230,0.223)

NP.B.3 <- c(0.259,0.260,0.233,0.233)
NP.BA.3 <- c(0.166,0.162)

NP.B.4 <- c(0.218,0.218,0.199,0.198)
NP.BA.4 <- c(0.166,0.162)

NP.B.5 <- c(0.195,0.194,0.181,0.181)
NP.BA.5 <- c(0.155,0.150)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Nonparametric Mean Conditional Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# SCENARIO 6
################################################################################
# cover
################################################################################

NP.1 <- c(0.86,0.94,0.96,0.88,0.76,0.40,0.44)
P.1 <- c(0.82,0.86,0.50,0.44,0.72,0.70)
SP.1 <- c(0.66,0.82,0.90,0.84,0.82)

NP.2 <- c(0.98,1.00,1.00,0.74,0.62,0.32,0.34)
P.2 <- c(0.66,0.64,0.38,0.40,0.46,0.50)
SP.2 <- c(0.56,0.82,0.86,0.86,0.82)

NP.3 <- c(0.90,0.96,0.72,0.60,0.46,0.46)
P.3 <- c(0.56,0.50,0.30,0.30,0.44,0.48)
SP.3 <- c(0.34,0.54,0.30,0.80,0.86,0.82,0.78)

NP.4 <- c(0.96,0.98,0.60,0.50,0.46,0.48)
P.4 <- c(0.36,0.36,0.12,0.18,0.32,0.32)
SP.4 <- c(0.48,0.64,0.78,0.68,0.60)

NP.5 <- c(0.96,0.98,0.54,0.50,0.32,0.30)
P.5 <- c(0.42,0.32,0.16,0.10,0.26,0.22)
SP.5 <- c(0.54,0.56,0.62,0.56,0.52)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.1,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Mean Coverage Probability\nNonparametric vs Parametric vs Semi-Parametric", cex.main=1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomleft", 
       legend = c("Nonparametric", "Parametric", "Semi-Parametric", "Nominal Level (0.95)"),
       col = c("blue", "red", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.70)

################################################################################
# Scenario 6 - alen
################################################################################

NP.1 <- c(0.267,0.298,0.306,0.266,0.269,0.258,0.249)
P.1 <- c(0.311,0.258,0.256,0.280,0.308,0.303)
SP.1 <- c(0.250,0.264,0.293,0.303,0.322)

NP.2 <- c(0.277,0.306,0.313,0.258,0.264,0.249,0.241)
P.2 <- c(0.286,0.249,0.260,0.261,0.287,0.282)
SP.2 <- c(0.395,0.257,0.289,0.294,0.312)

NP.3 <- c(0.286,0.292,0.238,0.244,0.227,0.220)
P.3 <- c(0.246,0.219,0.226,0.226,0.252,0.249)
SP.3 <- c(0.179,0.226,0.258,0.265,0.286)

NP.4 <- c(0.237,0.241,0.205,0.207,0.188,0.181)
P.4 <- c(0.217,0.193,0.198,0.198,0.213,0.210)
SP.4 <- c(0.172,0.190,0.214,0.217,0.226)

NP.5 <- c(0.187,0.190,0.166,0.168,0.161,0.155)
P.5 <- c(0.198,0.167,0.164,0.172,0.195,0.192)
SP.5 <- c(0.219,0.167,0.170,0.171,0.173)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.1,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Mean Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomleft", 
       legend = c("Nonparametric", "Parametric","Semi-Parametric"),
       col = c("blue", "red", "green"),
       lty = c(1, 1, 1, 1), 
       pch = c(20, 20, 20, 20), 
       lwd = 2,
       cex = 0.9)

################################################################################
# Scenario 6 - clen
################################################################################

NP.1 <- c(0.274,0.300,0.307,0.264,0.269,0.251,0.242)
P.1 <- c(0.309,0.253,0.244,0.264,0.306,0.301)
SP.1 <- c(0.274,0.258,0.291,0.301,0.319)

NP.2 <- c(0.278,0.306,0.313,0.255,0.258,0.237,0.232)
P.2 <- c(0.283,0.245,0.250,0.251,0.283,0.278)
SP.2 <- c(0.592,0.257,0.289,0.294,0.312)

NP.3 <- c(0.296,0.296,0.234,0.237,0.212,0.206)
P.3 <- c(0.240,0.209,0.214,0.222,0.247,0.246)
SP.3 <- c(0.262,0.223,0.254,0.258,0.278)

NP.4 <- c(0.237,0.240,0.205,0.205,0.184,0.177)
P.4 <- c(0.212,0.196,0.179,0.180,0.212,0.209)
SP.4 <- c(0.217,0.185,0.215,0.215,0.222)

NP.5 <- c(0.189,0.191,0.164,0.165,0.151,0.145)
P.5 <- c(0.195,0.165,0.146,0.153,0.188,0.184)
SP.5 <- c(0.319,0.163,0.165,0.165,0.167)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.1,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Mean Conditional Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomleft", 
       legend = c("Nonparametric", "Parametric","Semi-Parametric"),
       col = c("blue", "red", "green"),
       lty = c(1, 1, 1, 1), 
       pch = c(20, 20, 20, 20), 
       lwd = 2,
       cex = 0.9)

################################################################################
# Scenario 6 - Nonparametric cover
################################################################################

NP.A.1 <- c(0.86)
NP.B.1 <- c(0.94,0.96,0.88,0.76)
NP.BA.1 <- c(0.40,0.44)

NP.A.2 <- c(0.98)
NP.B.2 <- c(1.00,1.00,0.74,0.62)
NP.BA.2 <- c(0.32,0.34)

NP.B.3 <- c(0.90,0.96,0.72,0.60)
NP.BA.3 <- c(0.46,0.46)

NP.B.4 <- c(0.96,0.98,0.60,0.50)
NP.BA.4 <- c(0.46,0.48)

NP.B.5 <- c(0.96,0.98,0.54,0.50)
NP.BA.5 <- c(0.32,0.30)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Nonparametric Mean Coverage Probability", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian", "Nominal Level (0.95)"),
       col = c("red", "blue", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 6 - Nonparametric alen
################################################################################

NP.A.1 <- c(0.267)
NP.B.1 <- c(0.298,0.306,0.266,0.269)
NP.BA.1 <- c(0.258,0.249)

NP.A.2 <- c(0.277)
NP.B.2 <- c(0.306,0.313,0.258,0.264)
NP.BA.2 <- c(0.249,0.241)

NP.B.3 <- c(0.286,0.292,0.238,0.244)
NP.BA.3 <- c(0.227,0.220)

NP.B.4 <- c(0.237,0.241,0.205,0.207)
NP.BA.4 <- c(0.188,0.181)

NP.B.5 <- c(0.187,0.190,0.166,0.168)
NP.BA.5 <- c(0.161,0.155)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Nonparametric Mean Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 6 - Nonparametric clen
################################################################################

NP.A.1 <- c(0.274)
NP.B.1 <- c(0.300,0.307,0.264,0.269)
NP.BA.1 <- c(0.251,0.242)

NP.A.2 <- c(0.278)
NP.B.2 <- c(0.306,0.313,0.255,0.258)
NP.BA.2 <- c(0.237,0.232)

NP.B.3 <- c(0.296,0.296,0.234,0.237)
NP.BA.3 <- c(0.212,0.206)

NP.B.4 <- c(0.237,0.240,0.205,0.205)
NP.BA.4 <- c(0.184,0.177)

NP.B.5 <- c(0.189,0.191,0.164,0.165)
NP.BA.5 <- c(0.151,0.145)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Nonparametric Mean Conditional Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# SCENARIO 7
################################################################################
# cover
################################################################################

NP.1 <- c(0.86,0.96,0.96,0.86,0.94,0.96,0.90)
P.1 <- c(0.88,0.90,0.94,0.94,0.96,0.94)
SP.1 <- c(0.70,0.92,0.92,0.96,0.98)

NP.2 <- c(0.88,0.92,0.92,0.94,0.96,0.98,0.98)
P.2 <- c(0.88,0.88,0.88,0.90,1.00,1.00)
SP.2 <- c(0.76,0.90,0.90,1.00,1.00)

NP.3 <- c(0.96,0.96,0.94,0.94,0.94,0.92)
P.3 <- c(0.88,0.88,0.92,0.92,0.96,0.96)
SP.3 <- c(0.72,0.94,0.94,0.96,0.98)

NP.4 <- c(0.96,0.98,0.94,0.96,0.92,0.90)
P.4 <- c(0.98,0.92,0.94,0.92,0.94,0.94)
SP.4 <- c(0.72,0.96,0.94,0.94,1.00)

NP.5 <- c(0.96,0.94,0.94,0.94,0.90,0.90)
P.5 <- c(0.88,0.86,0.84,0.84,0.92,0.92)
SP.5 <- c(0.78,0.96,0.94,0.94,0.96)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.8,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Mean Coverage Probability\nNonparametric vs Parametric vs Semi-Parametric", cex.main=1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "Semi-Parametric", "Nominal Level (0.95)"),
       col = c("blue", "red", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.70)

################################################################################
# Scenario 7 - alen
################################################################################

NP.1 <- c(0.251,0.285,0.288,0.243,0.244,0.216,0.210)
P.1 <- c(0.260,0.265,0.267,0.283,0.250,0.246)
SP.1 <- c(0.197,0.225,0.229,0.230,0.232)

NP.2 <- c(0.237,0.266,0.268,0.227,0.227,0.205,0.198)
P.2 <- c(0.243,0.245,0.245,0.253,0.238,0.234)
SP.2 <- c(0.163,0.205,0.211,0.210,0.210)

NP.3 <- c(0.242,0.243,0.210,0.209,0.185,0.181)
P.3 <- c(0.216,0.220,0.221,0.228,0.210,0.207)
SP.3 <- c(0.158,0.180,0.182,0.181,0.182)

NP.4 <- c(0.199,0.200,0.175,0.176,0.159,0.155)
P.4 <- c(0.184,0.184,0.184,0.189,0.182,0.178)
SP.4 <- c(0.120,0.160,0.156,0.155,0.156)

NP.5 <- c(0.177,0.177,0.156,0.156,0.143,0.138)
P.5 <- c(0.168,0.165,0.165,0.168,0.164,0.162)
SP.5 <- c(0.115,0.145,0.142,0.141,0.141)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.30),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Mean Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric","Semi-Parametric"),
       col = c("blue", "red", "green"),
       lty = c(1, 1, 1, 1), 
       pch = c(20, 20, 20, 20), 
       lwd = 2,
       cex = 0.9)

################################################################################
# Scenario 7 - clen
################################################################################

NP.1 <- c(0.253,0.289,0.289,0.247,0.246,0.218,0.214)
P.1 <- c(0.261,0.264,0.265,0.282,0.253,0.250)
SP.1 <- c(0.202,0.228,0.230,0.231,0.234)

NP.2 <- c(0.237,0.265,0.265,0.229,0.228,0.205,0.198)
P.2 <- c(0.239,0.242,0.242,0.252,0.238,0.234)
SP.2 <- c(0.190,0.206,0.210,0.210,0.210)

NP.3 <- c(0.242,0.243,0.210,0.209,0.185,0.181)
P.3 <- c(0.216,0.220,0.221,0.228,0.210,0.207)
SP.3 <- c(0.165,0.181,0.183,0.187,0.183)

NP.4 <- c(0.200,0.200,0.177,0.177,0.160,0.156)
P.4 <- c(0.184,0.184,0.184,0.189,0.182,0.179)
SP.4 <- c(0.128,0.161,0.157,0.182,0.156)

NP.5 <- c(0.177,0.177,0.156,0.156,0.142,0.138)
P.5 <- c(0.167,0.164,0.164,0.167,0.164,0.161)
SP.5 <- c(0.122,0.145,0.141,0.141,0.141)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.3),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Mean Conditional Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric","Semi-Parametric"),
       col = c("blue", "red", "green"),
       lty = c(1, 1, 1, 1), 
       pch = c(20, 20, 20, 20), 
       lwd = 2,
       cex = 0.9)

################################################################################
# Scenario 7 - Non-parametric cover
################################################################################

NP.A.1 <- c(0.86)
NP.B.1 <- c(0.96,0.96,0.86,0.94)
NP.BA.1 <- c(0.96,0.90)

NP.A.2 <- c(0.88)
NP.B.2 <- c(0.92,0.91,0.94,0.96)
NP.BA.2 <- c(0.98,0.98)

NP.B.3 <- c(0.96,0.96,0.94,0.94)
NP.BA.3 <- c(0.94,0.92)

NP.B.4 <- c(0.96,0.98,0.94,0.96)
NP.BA.4 <- c(0.92,0.90)

NP.B.5 <- c(0.96,0.94,0.94,0.94)
NP.BA.5 <- c(0.90,0.90)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.4,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Non-Parametric Mean Coverage Probability", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian", "Nominal Level (0.95)"),
       col = c("red", "blue", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 7 - Semi-parametric alen
################################################################################

NP.A.1 <- c(0.252)
NP.B.1 <- c(0.285,0.288,0.243,0.244)
NP.BA.1 <- c(0.216,0.210)

NP.A.2 <- c(0.237)
NP.B.2 <- c(0.266,0.268,0.227,0.227)
NP.BA.2 <- c(0.205,0.198)

NP.B.3 <- c(0.242,0.243,0.209,0.208)
NP.BA.3 <- c(0.184,0.180)

NP.B.4 <- c(0.199,0.200,0.175,0.176)
NP.BA.4 <- c(0.159,0.155)

NP.B.5 <- c(0.177,0.177,0.156,0.156)
NP.BA.5 <- c(0.143,0.138)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.3),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Non-Parametric Mean Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# SCENARIO 8
################################################################################
# cover
################################################################################

NP.1 <- c(0.88,0.94,0.96,0.88,0.92,0.90,0.86)
P.1 <- c(0.72,0.90,0.88,0.92,0.86,0.86)
BC.1 <- c(0.92,0.96,0.96,0.94)
SP.1 <- c(0.78,0.54,0.64,0.68,0.64)

NP.2 <- c(0.94,0.98,0.98,0.92,0.94,0.94,0.94)
P.2 <- c(0.80,0.96,0.94,0.94,0.88,0.82)
BC.2 <- c(0.94,0.98,0.98,0.98)
SP.2 <- c(0.74,0.50,0.68,0.78,0.78)

NP.3 <- c(0.98,1.00,0.98,0.98,0.88,0.88)
P.3 <- c(0.92,0.96,0.96,0.94,0.92,0.88)
BC.3 <- c(0.94,0.94,0.96,0.96)
SP.3 <- c(0.82,0.62,0.76,0.74,0.74)

NP.4 <- c(0.90,0.94,0.86,0.86,0.90,0.86)
P.4 <- c(0.76,0.84,0.84,0.88,0.92,0.86)
BC.4 <- c(0.94,0.98,0.98,0.96)
SP.4 <- c(0.80,0.38,0.40,0.40,0.42)

NP.5 <- c(0.92,0.94,0.84,0.86,0.88,0.84)
P.5 <- c(0.78,0.78,0.80,0.80,0.86,0.82)
BC.5 <- c(0.88,0.94,0.92,0.94)
SP.5 <- c(0.82,0.18,0.28,0.26,0.24)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.1,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Mean Coverage Probability\nNonparametric vs Parametric vs Semi-Parametric", cex.main=1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)
lines(x.axis, BC, type = "b", col = "orange", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "With Box-Cox Transformation","Semi-Parametric", "Nominal Level (0.95)"),
       col = c("blue", "red",  "orange","green", "black"),
       lty = c(1, 1, 1, 1, 2), 
       pch = c(20, 20, 20, 20, NA), 
       lwd = 2,
       cex = 0.70)

################################################################################
# Scenario 8 - alen
################################################################################

NP.1 <- c(0.234,0.262,0.268,0.220,0.224,0.195,0.188)
P.1 <- c(0.218,0.232,0.233,0.248,0.209,0.204)
BC.1 <- c(0.236,0.268,0.271,0.289)
SP.1 <- c(0.161,0.165,0.188,0.190,0.190)

NP.2 <- c(0.223,0.249,0.253,0.206,0.206,0.183,0.177)
P.2 <- c(0.204,0.217,0.217,0.227,0.200,0.195)
BC.2 <- c(0.224,0.251,0.253,0.265)
SP.2 <- c(0.148,0.151,0.180,0.181,0.182)

NP.3 <- c(0.233,0.236,0.197,0.198,0.168,0.163)
P.3 <- c(0.194,0.205,0.203,0.212,0.187,0.183)
BC.3 <- c(0.211,0.238,0.240,0.251)
SP.3 <- c(0.137,0.139,0.182,0.179,0.179)

NP.4 <- c(0.190,0.192,0.160,0.161,0.148,0.143)
P.4 <- c(0.160,0.166,0.166,0.172,0.161,0.159)
BC.4 <- c(0.176,0.193,0.194,0.200)
SP.4 <- c(0.107,0.120,0.137,0.135,0.137)

NP.5 <- c(0.165,0.167,0.141,0.141,0.129,0.124)
P.5 <- c(0.140,0.139,0.140,0.143,0.137,0.135)
BC.5 <- c(0.151,0.164,0.164,0.166)
SP.5 <- c(0.103,0.116,0.124,0.123,0.124)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.30),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Mean Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)
lines(x.axis, BC, type = "b", col = "orange", lwd = 2, pch = 20)

legend("bottomleft", 
       legend = c("Nonparametric", "Parametric","With Box-Cox Transformation","Semi-Parametric"),
       col = c("blue", "red","orange", "green"),
       lty = c(1, 1, 1, 1), 
       pch = c(20, 20, 20, 20), 
       lwd = 2,
       cex = 0.9)

################################################################################
# Scenario 8 - clen
################################################################################

NP.1 <- c(0.243,0.268,0.271,0.228,0.229,0.201,0.195)
P.1 <- c(0.235,0.239,0.241,0.252,0.216,0.212)
BC.1 <- c(0.242,0.272,0.275,0.293)
SP.1 <- c(0.178,0.182,0.212,0.209,0.210)

NP.2 <- c(0.226,0.251,0.255,0.211,0.209,0.186,0.180)
P.2 <- c(0.214,0.219,0.219,0.229,0.205,0.203)
BC.2 <- c(0.225,0.250,0.253,0.264)
SP.2 <- c(0.162,0.164,0.191,0.188,0.189)

NP.3 <- c(0.234,0.236,0.198,0.199,0.170,0.166)
P.3 <- c(0.197,0.208,0.205,0.215,0.190,0.187)
BC.3 <- c(0.210,0.239,0.240,0.251)
SP.3 <- c(0.139,0.148,0.193,0.189,0.190)

NP.4 <- c(0.192,0.192,0.163,0.163,0.149,0.145)
P.4 <- c(0.165,0.170,0.170,0.175,0.164,0.162)
BC.4 <- c(0.175,0.193,0.193,0.199)
SP.4 <- c(0.106,0.134,0.155,0.149,0.152)

NP.5 <- c(0.167,0.168,0.145,0.144,0.131,0.126)
P.5 <- c(0.145,0.145,0.145,0.148,0.140,0.138)
BC.5 <- c(0.152,0.164,0.163,0.166)
SP.5 <- c(0.106,0.116,0.124,0.123,0.124)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.30),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Mean Conditional Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)
lines(x.axis, BC, type = "b", col = "orange", lwd = 2, pch = 20)

legend("bottomleft", 
       legend = c("Nonparametric", "Parametric","With Box-Cox Transformation", "Semi-Parametric"),
       col = c("blue", "red", "orange", "green"),
       lty = c(1, 1, 1, 1), 
       pch = c(20, 20, 20, 20), 
       lwd = 2,
       cex = 0.9)

################################################################################
# Scenario 8 - Nonparametric cover
################################################################################

NP.A.1 <- c(0.88)
NP.B.1 <- c(0.94,0.96,0.88,0.92)
NP.BA.1 <- c(0.90,0.86)

NP.A.2 <- c(0.94)
NP.B.2 <- c(0.98,0.98,0.92,0.94)
NP.BA.2 <- c(0.94,0.94)

NP.B.3 <- c(0.98,1.00,0.98,0.98)
NP.BA.3 <- c(0.88,0.88)

NP.B.4 <- c(0.90,0.94,0.86,0.86)
NP.BA.4 <- c(0.90,0.86)

NP.B.5 <- c(0.92,0.94,0.84,0.86)
NP.BA.5 <- c(0.88,0.84)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.8,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Nonparametric Mean Coverage Probability", cex.main = 1)

lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)
lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomleft", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian", "Nominal Level (0.95)"),
       col = c("red", "blue", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 8 - Nonparametric alen
################################################################################

NP.A.1 <- c(0.234)
NP.B.1 <- c(0.262,0.268,0.220,0.224)
NP.BA.1 <- c(0.195,0.188)

NP.A.2 <- c(0.223)
NP.B.2 <- c(0.249,0.253,0.206,0.206)
NP.BA.2 <- c(0.183,0.177)

NP.B.3 <- c(0.233,0.236,0.197,0.198)
NP.BA.3 <- c(0.168,0.163)

NP.B.4 <- c(0.190,0.192,0.160,0.161)
NP.BA.4 <- c(0.148,0.143)

NP.B.5 <- c(0.165,0.167,0.141,0.141)
NP.BA.5 <- c(0.129,0.124)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Nonparametric Mean Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 8 - Nonparametric clen
################################################################################

NP.A.1 <- c(0.243)
NP.B.1 <- c(0.268,0.271,0.228,0.229)
NP.BA.1 <- c(0.201,0.195)

NP.A.2 <- c(0.226)
NP.B.2 <- c(0.251,0.255,0.211,0.209)
NP.BA.2 <- c(0.186,0.180)

NP.B.3 <- c(0.234,0.236,0.198,0.199)
NP.BA.3 <- c(0.170,0.166)

NP.B.4 <- c(0.192,0.192,0.163,0.163)
NP.BA.4 <- c(0.149,0.145)

NP.B.5 <- c(0.167,0.168,0.145,0.144)
NP.BA.5 <- c(0.131,0.126)

NP.A <- c(mean(NP.A.1), mean(NP.A.2))
NP.B <- c(mean(NP.B.1), mean(NP.B.2), mean(NP.B.3), mean(NP.B.4), mean(NP.B.5))
NP.BA <- c(mean(NP.BA.1), mean(NP.BA.2), mean(NP.BA.3), mean(NP.BA.4), mean(NP.BA.5))

x.axis <- 1:5

plot(x.axis, NP.B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.35),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Nonparametric Mean Conditional Average Length", cex.main = 1)

lines(1:2, NP.A, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, NP.BA, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# SCENARIO 9
################################################################################
# cover
################################################################################

NP.1 <- c(0.92,0.96,0.98,0.90,0.90,0.86,0.76)
P.1 <- c(0.88,0.96,0.92,0.94,0.88,0.86)
SP.1 <- c(0.70,0.86,0.82,0.70)

NP.2 <- c(0.98,0.98,0.98,0.96,0.96,0.84,0.82)
P.2 <- c(0.88,0.94,0.96,0.92,0.88,0.90)
SP.2 <- c(0.68,0.92,0.86,0.70)

NP.3 <- c(0.94,0.94,0.86,0.88,0.82,0.76)
P.3 <- c(0.88,0.92,0.90,0.90,0.84,0.82)
SP.3 <- c(0.54,0.88,0.90,0.74)

NP.4 <- c(1.00,1.00,0.94,0.94,0.82,0.82)
P.4 <- c(0.88,0.94,0.94,0.96,0.88,0.88)
SP.4 <- c(0.54,0.88,0.70,0.58)

NP.5 <- c(0.88,0.88,0.90,0.90,0.86,0.84)
P.5 <- c(0.86,0.88,0.88,0.90,0.86,0.86)
SP.5 <- c(0.52,0.68,0.54,0.52)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

par(mfrow=c(1,1))
x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.3,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Mean Coverage Probability\nNonparametric vs Parametric vs Semi-Parametric", cex.main=1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric","Semi-Parametric", "Nominal Level (0.95)"),
       col = c("blue", "red", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.70)

################################################################################
# Scenario 9 - alen
################################################################################

NP.1 <- c(0.283,0.319,0.320,0.278,0.277,0.236,0.229)
P.1 <- c(0.285,0.291,0.289,0.303,0.278,0.274)
SP.1 <- c(0.242,0.291,0.286,0.285)

NP.2 <- c(0.274,0.306,0.307,0.273,0.273,0.223,0.216)
P.2 <- c(0.260,0.281,0.280,0.286,0.257,0.254)
SP.2 <- c(0.218,0.275,0.268,0.265)

NP.3 <- c(0.264,0.263,0.236,0.235,0.194,0.188)
P.3 <- c(0.231,0.248,0.248,0.255,0.223,0.220)
SP.3 <- c(0.188,0.258,0.252,0.253)

NP.4 <- c(0.222,0.221,0.199,0.198,0.166,0.162)
P.4 <- c(0.198,0.219,0.218,0.224,0.194,0.191)
SP.4 <- c(0.166,0.215,0.211,0.211)

NP.5 <- c(0.196,0.195,0.181,0.181,0.153,0.149)
P.5 <- c(0.180,0.192,0.193,0.194,0.178,0.176)
SP.5 <- c(0.149,0.179,0.175,0.174)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.30),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Mean Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric","Semi-Parametric"),
       col = c("blue", "red", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.9)

################################################################################
# Scenario 9 - clen
################################################################################

NP.1 <- c(0.285,0.320,0.320,0.280,0.279,0.239,0.231)
P.1 <- c(0.289,0.292,0.290,0.304,0.283,0.280)
SP.1 <- c(0.250,0.298,0.291,0.292)

NP.2 <- c(0.275,0.307,0.308,0.274,0.274,0.230,0.222)
P.2 <- c(0.264,0.284,0.283,0.284,0.260,0.256)
SP.2 <- c(0.225,0.280,0.273,0.267)

NP.3 <- c(0.265,0.265,0.239,0.239,0.198,0.191)
P.3 <- c(0.233,0.252,0.252,0.259,0.227,0.224)
SP.3 <- c(0.194,0.260,0.255,0.251)

NP.4 <- c(0.222,0.221,0.200,0.199,0.168,0.164)
P.4 <- c(0.200,0.220,0.219,0.224,0.196,0.193)
SP.4 <- c(0.171,0.220,0.215,0.211)

NP.5 <- c(0.197,0.197,0.183,0.183,0.155,0.151)
P.5 <- c(0.182,0.194,0.195,0.196,0.180,0.177)
SP.5 <- c(0.154,0.184,0.182,0.179)

NP <- c(mean(NP.1), mean(NP.2), mean(NP.3), mean(NP.4), mean(NP.5))
P <- c(mean(P.1), mean(P.2), mean(P.3), mean(P.4), mean(P.5))
BC <- c(mean(BC.1), mean(BC.2), mean(BC.3), mean(BC.4), mean(BC.5))
SP <- c(mean(SP.1), mean(SP.2), mean(SP.3), mean(SP.4), mean(SP.5))

x.axis <- 1:5

plot(x.axis, NP, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.30),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Mean Conditional Average Length\nNonparametric vs Parametric vs Semi-Parametric", cex.main = 1)

lines(x.axis, P, type = "b", col = "red", lwd = 2, pch = 20)
lines(x.axis, SP, type = "b", col = "green", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Nonparametric", "Parametric", "Semi-Parametric"),
       col = c("blue", "red", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.9)

################################################################################
# Scenario 9 - Nonparametric and parametric cover
################################################################################

A.1 <- c(0.92,0.88)
B.1 <- c(0.96,0.98,0.90,0.90,0.96,0.92,0.94)
BA.1 <- c(0.86,0.76,0.88,0.86)

A.2 <- c(0.98,0.88)
B.2 <- c(0.98,0.98,0.96,0.96,0.94,0.96,0.92)
BA.2 <- c(0.84,0.82,0.88,0.90)

A.3 <- c(0.88)
B.3 <- c(0.94,0.94,0.86,0.88,0.92,0.90,0.90)
BA.3 <- c(0.82,0.76,0.84,0.82)

A.4 <- c(0.88)
B.4 <- c(1.00,1.00,0.94,0.94,0.94,0.94,0.96)
BA.4 <- c(0.82,0.82,0.88,0.88)

A.5 <- c(0.86)
B.5 <- c(0.88,0.88,0.90,0.90,0.88,0.88,0.90)
BA.5 <- c(0.86,0.84,0.86,0.86)

A <- c(mean(A.1), mean(A.2), mean(A.3), mean(A.4), mean(A.5))
B <- c(mean(B.1), mean(B.2), mean(B.3), mean(B.4), mean(B.5))
BA <- c(mean(BA.1), mean(BA.2), mean(BA.3), mean(BA.4), mean(BA.5))

par(mfrow=c(1,1))

plot(x.axis, B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0.75,1),
     xlab = "Sample Size Scenario", ylab = "Mean Coverage Probability", 
     main = "Nonparametric and Parametric Mean Coverage Probability", cex.main = 1)

lines(x.axis, BA, type = "b", col = "green", lwd = 2, pch = 20)
lines(x.axis, A, type = "b", col = "red", lwd = 2, pch = 20)

abline(h = 0.95, col = "black", lty = 2, lwd = 2)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian", "Nominal Level (0.95)"),
       col = c("red", "blue", "green", "black"),
       lty = c(1, 1, 1, 2), 
       pch = c(20, 20, 20, NA), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 9 - Nonparametric and parametric alen
################################################################################

A.1 <- c(0.283,0.285)
B.1 <- c(0.319,0.320,0.278,0.277,0.291,0.289,0.303)
BA.1 <- c(0.236,0.229,0.278,0.274)

A.2 <- c(0.274,0.260)
B.2 <- c(0.306,0.307,0.273,0.273,0.281,0.280,0.286)
BA.2 <- c(0.223,0.216,0.257,0.254)

A.3 <- c(0.231)
B.3 <- c(0.264,0.263,0.236,0.235,0.248,0.248,0.255)
BA.3 <- c(0.194,0.188,0.223,0.220)

A.4 <- c(0.198)
B.4 <- c(0.222,0.221,0.199,0.198,0.219,0.218,0.224)
BA.4 <- c(0.166,0.162,0.194,0.191)

A.5 <- c(0.180)
B.5 <- c(0.196,0.195,0.181,0.181,0.192,0.193,0.194)
BA.5 <- c(0.153,0.149,0.178,0.176)

A <- c(mean(A.1), mean(A.2), mean(A.3), mean(A.4), mean(A.5))
B <- c(mean(B.1), mean(B.2), mean(B.3), mean(B.4), mean(B.5))
BA <- c(mean(BA.1), mean(BA.2), mean(BA.3), mean(BA.4), mean(BA.5))

x.axis <- 1:5

plot(x.axis, B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.3),
     xlab = "Sample Size Scenario", ylab = "Mean Average Length", 
     main = "Nonparametric and Parametric Mean Average Length", cex.main = 1)

lines(x.axis, BA, type = "b", col = "green", lwd = 2, pch = 20)
lines(x.axis, A, type = "b", col = "red", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)

################################################################################
# Scenario 9 - Nonparametric and parametric clen
################################################################################

A.1 <- c(0.285,0.289)
B.1 <- c(0.320,0.320,0.280,0.279,0.292,0.290,0.304)
BA.1 <- c(0.239,0.231,0.283,0.280)

A.2 <- c(0.275,0.264)
B.2 <- c(0.307,0.308,0.274,0.274,0.284,0.283,0.284)
BA.2 <- c(0.230,0.222,0.260,0.256)

A.3 <- c(0.233)
B.3 <- c(0.265,0.265,0.239,0.239,0.252,0.252,0.259)
BA.3 <- c(0.198,0.191,0.227,0.224)

A.4 <- c(0.200)
B.4 <- c(0.222,0.221,0.200,0.199,0.220,0.219,0.224)
BA.4 <- c(0.168,0.164,0.196,0.193)

A.5 <- c(0.182)
B.5 <- c(0.197,0.197,0.183,0.183,0.194,0.195,0.196)
BA.5 <- c(0.155,0.151,0.180,0.177)

A <- c(mean(A.1), mean(A.2), mean(A.3), mean(A.4), mean(A.5))
B <- c(mean(B.1), mean(B.2), mean(B.3), mean(B.4), mean(B.5))
BA <- c(mean(BA.1), mean(BA.2), mean(BA.3), mean(BA.4), mean(BA.5))

x.axis <- 1:5

plot(x.axis, B, type = "b", col = "blue", lwd = 2, pch = 20,
     ylim = c(0,0.3),
     xlab = "Sample Size Scenario", ylab = "Mean Conditional Average Length", 
     main = "Nonparametric and Parametric Mean Conditional Average Length", cex.main = 1)

lines(x.axis, BA, type = "b", col = "green", lwd = 2, pch = 20)
lines(x.axis, A, type = "b", col = "red", lwd = 2, pch = 20)

legend("bottomright", 
       legend = c("Asymptotic", "Bootstrap", "Bayesian"),
       col = c("red", "blue", "green"),
       lty = c(1, 1, 1), 
       pch = c(20, 20, 20), 
       lwd = 2,
       cex = 0.65)