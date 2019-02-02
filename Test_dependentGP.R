rm(list = ls()) #Start with a clean working directory

###Load Library (set the path that you save the package. You end with the package in the file path)
path.package <- "/Users/sam/Documents/Postdoc/Software/dependentGP/"
devtools::load_all(path.package, export_all = TRUE)
devtools::document(path.package)
library(dependentGP)

# install.packages("dependentGP_1.0.tar.gz", lib = "/home/sib2/R", type = "source", repo = NULL)
# devtools::build_win(pkg = path.package)
# devtools::use_vignette("spCP-example", pkg = path.package)
# devtools::use_build_ignore("spCP-example.pdf", pkg = path.package)
# devtools::use_build_ignore(files = "cran-comments.md", pkg = path.package)
# devtools::use_build_ignore("NEWS.md", pkg = path.package)
# devtools::use_build_ignore("README.md", pkg = path.package)
# devtools::build_win(pkg = path.package)
# devtools::spell_check()
# devtools::release(path.package) # for submitting package for the first time
# devtools::submit_cran(path.package) # submit to CRAN without going through checks

###Simulate data for MCMC sampler
load(file = "/Users/sam/Box Sync/Postdoc/Projects/VAE/Data/latent.RData")
Latent <- Latent[Latent$data == 2, ]
ids <- NULL
for (i in 1:length(unique(Latent$id))) {
  temp <- Latent[Latent$id == unique(Latent$id)[i], ]
  if (dim(temp)[1] == 1) ids <- c(ids, unique(Latent$id)[i])
  if (length(temp$time) != length(unique(temp$time))) ids <- c(ids, unique(Latent$id)[i])
}
Latent <- Latent[!(Latent$id %in% ids), ]
Latent <-  Latent[, -3]
Latent[, 2] <- Latent[, 2] / 365
K <- dim(Latent[, -c(1, 2)])[2]
dat <- Latent[Latent$id == 320, ]
Time <- dat$time
Y <- matrix(c(dat[, 3], dat[, 4]), ncol = 1)
T <- length(Time)
Min <- apply(Latent[, -c(1, 2)], 2, min)
Max <- apply(Latent[, -c(1, 2)], 2, max)

###Play around with the Gaussian process prior (squared exponential)
library(mvtnorm)
library(pscl)
library(bayesm)
Seq <- seq(min(Time), max(Time), length.out = 250)
TimeDist <- outer(Seq, Seq, "-")^2 / 2
minDiff <- min(TimeDist[TimeDist > 0])
maxDiff <- max(TimeDist[TimeDist > 0])
A <- sqrt(-maxDiff / log(0.25)) #longest diff goes down to 25%
B <- sqrt(-minDiff / log(0.999)) #shortest diff goes up to 99.9%
Lower <- min(A, B)
Upper <- max(A, B)
Mean <- (Lower + Upper) / 2
CovUpper <- H(Upper, TimeDist, 1, 250)
SamplesUpper <- rmvnorm(10, rep(0, dim(CovUpper)[1]), CovUpper)
CovLower <- H(Lower, TimeDist, 1, 250)
SamplesLower <- rmvnorm(10, rep(0, dim(CovUpper)[1]), CovLower)
CovMean <- H(Mean, TimeDist, 1, 250)
SamplesMean <- rmvnorm(10, rep(0, dim(CovUpper)[1]), CovMean)
par(mfcol = c(1, 3))
plot(Seq, SamplesUpper[1, ], type = "n", lwd = 3, ylim = c(-5, 5), xlab = "Time", ylab = "Latent Space", main = "Prior Samples of the GP (Upper Bound)")
for (i in 1:10) lines(Seq, SamplesUpper[i, ], lwd = 3, col = "gray")
abline(h = 0, lty = 4, lwd = 4)
plot(Seq, SamplesLower[1, ], type = "n", lwd = 3, ylim = c(-5, 5), xlab = "Time", ylab = "Latent Space", main = "Prior Samples of the GP (Lower Bound)")
for (i in 1:10) lines(Seq, SamplesLower[i, ], lwd = 3, col = "gray")
abline(h = 0, lty = 4, lwd = 4)
plot(Seq, SamplesMean[1, ], type = "n", lwd = 3, ylim = c(-5, 5), xlab = "Time", ylab = "Latent Space", main = "Prior Samples of the GP (Mean)")
for (i in 1:10) lines(Seq, SamplesMean[i, ], lwd = 3, col = "gray")
abline(h = 0, lty = 4, lwd = 4)

###Play around with the Gaussian process prior (exponential)
# Seq <- seq(min(Time), max(Time), length.out = 250)
# TimeDist <- abs(outer(Seq, Seq, "-"))
# minDiff <- min(TimeDist[TimeDist > 0])
# maxDiff <- max(TimeDist[TimeDist > 0])
# A <- -log(0.25) / maxDiff #longest diff goes down to 25%
# B <- -log(0.999) / minDiff #shortest diff goes up to 99.9%
# Lower <- min(A, B)
# Upper <- max(A, B)
# Mean <- (Lower + Upper) / 2
# CovUpper <- H(Upper, TimeDist, 0)
# SamplesUpper <- rmvnorm(10, rep(0, dim(CovUpper)[1]), CovUpper)
# CovLower <- H(Lower, TimeDist, 0)
# SamplesLower <- rmvnorm(10, rep(0, dim(CovUpper)[1]), CovLower)
# CovMean <- H(Mean, TimeDist, 0)
# SamplesMean <- rmvnorm(10, rep(0, dim(CovUpper)[1]), CovMean)
# par(mfcol = c(1, 3))
# plot(Seq, SamplesUpper[1, ], type = "n", lwd = 3, ylim = c(-5, 5), xlab = "Time", ylab = "Latent Space", main = "Prior Samples of the GP (Upper Bound)")
# for (i in 1:10) lines(Seq, SamplesUpper[i, ], lwd = 3, col = "gray")
# abline(h = 0, lty = 4, lwd = 4)
# plot(Seq, SamplesLower[1, ], type = "n", lwd = 3, ylim = c(-5, 5), xlab = "Time", ylab = "Latent Space", main = "Prior Samples of the GP (Lower Bound)")
# for (i in 1:10) lines(Seq, SamplesLower[i, ], lwd = 3, col = "gray")
# abline(h = 0, lty = 4, lwd = 4)
# plot(Seq, SamplesMean[1, ], type = "n", lwd = 3, ylim = c(-5, 5), xlab = "Time", ylab = "Latent Space", main = "Prior Samples of the GP (Mean)")
# for (i in 1:10) lines(Seq, SamplesMean[i, ], lwd = 3, col = "gray")
# abline(h = 0, lty = 4, lwd = 4)

###Play around with the Gaussian process prior (ar1)
# Seq <- seq(min(Time), max(Time), length.out = 250)
# TimeDist <- abs(outer(Seq, Seq, "-"))
# minDiff <- min(TimeDist[TimeDist > 0])
# maxDiff <- max(TimeDist[TimeDist > 0])
# A <- (0.5)^(1 / maxDiff) #longest diff goes down to 25%
# B <- (0.5)^(1 / minDiff) #shortest diff goes up to 99.9%
# Lower <- min(A, B)
# Upper <- max(A, B)
# Mean <- (Lower + Upper) / 2
# CovUpper <- H(Upper, TimeDist, 0)
# SamplesUpper <- rmvnorm(10, rep(0, dim(CovUpper)[1]), CovUpper)
# CovLower <- H(Lower, TimeDist, 0)
# SamplesLower <- rmvnorm(10, rep(0, dim(CovUpper)[1]), CovLower)
# CovMean <- H(Mean, TimeDist, 0)
# SamplesMean <- rmvnorm(10, rep(0, dim(CovUpper)[1]), CovMean)
# par(mfcol = c(1, 3))
# plot(Seq, SamplesUpper[1, ], type = "n", lwd = 3, ylim = c(-5, 5), xlab = "Time", ylab = "Latent Space", main = "Prior Samples of the GP (Upper Bound)")
# for (i in 1:10) lines(Seq, SamplesUpper[i, ], lwd = 3, col = "gray")
# abline(h = 0, lty = 4, lwd = 4)
# plot(Seq, SamplesLower[1, ], type = "n", lwd = 3, ylim = c(-5, 5), xlab = "Time", ylab = "Latent Space", main = "Prior Samples of the GP (Lower Bound)")
# for (i in 1:10) lines(Seq, SamplesLower[i, ], lwd = 3, col = "gray")
# abline(h = 0, lty = 4, lwd = 4)
# plot(Seq, SamplesMean[1, ], type = "n", lwd = 3, ylim = c(-5, 5), xlab = "Time", ylab = "Latent Space", main = "Prior Samples of the GP (Mean)")
# for (i in 1:10) lines(Seq, SamplesMean[i, ], lwd = 3, col = "gray")
# abline(h = 0, lty = 4, lwd = 4)

###Initial values
# Starting <- list(T = diag(K),
#                  Sigma2 = rep(1, K),
#                  Gamma = 0)

###Hyperparameters
# Hypers <- list(T = list(Xi = K + 1, Psi = diag(K)),
#                Sigma2 = list(Alpha = 0.001, Beta = 0.001))

###Metropolis tuners
# Tuning <- list(Phi = rep(1, K), T = rep(1, ((K + 1) * K) / 2))

###MCMC inputs
# MCMC <- list(NBurn = 10000, NSims = 250000, NThin = 25, NPilot = 10)
# MCMC <- list(NBurn = 10000, NSims = 10000, NThin = 4, NPilot = 10)

###Fit sampler
reg <- lmcGP(Y, Time)





P <- 50
NewTimes <- seq(0, Time[T] + 1, length.out = P)

derv <- predict(reg, NewTimes = Time, type = "derivative")
MeanDerv <- UpperDerv <- LowerDerv <- matrix(nrow = T, ncol = K)
for  (t in 1:T) {
  MeanDerv[t, ] <- apply(derv[t, , ], 1, mean)
  UpperDerv[t, ] <- apply(derv[t, , ], 1, function(x) quantile(x, probs = 0.975))
  LowerDerv[t, ] <- apply(derv[t, , ], 1, function(x) quantile(x, probs = 0.025))
}
Norm <- Norm2 <- matrix(nrow = dim(derv)[3], ncol = T)
for (s in 1:dim(derv)[3]) {
  Norm[s, ] <- apply(derv[ , , s], 1, function(x) sqrt(x %*% x))
  Norm2[s, ] <- apply(derv[ , , s], 1, function(x) sqrt((x - derv[ , , s][1, ]) %*% (x - derv[ , , s][1, ])))
}

pdf("/Users/sam/Desktop/VAE/4_derv.pdf", height = 6, width = 8)
plot(Time, MeanDerv[, 1], type = "l", ylim = c(-0.01, 0.01), lwd = 3, ylab = "Derivative")
lines(Time, UpperDerv[, 1], lty = 3, lwd = 3)
lines(Time, LowerDerv[, 1], lty = 3, lwd = 3)
lines(Time, MeanDerv[, 2], col = 2, lwd = 3)
lines(Time, UpperDerv[, 2], lty = 3, col = 2, lwd = 3)
lines(Time, LowerDerv[, 2], lty = 3, col = 2, lwd = 3)
abline(h = 0, lwd = 3, lty = 4, col = "gray")
dev.off()

pdf("/Users/sam/Desktop/VAE/5_norm.pdf", height = 6, width = 8)
plot(Time, apply(Norm, 2, mean), type = "l", ylim = c(0, 0.025), lwd = 3, ylab = "Gradient Norm")
lines(Time, apply(Norm, 2, function(x) quantile(x, probs = 0.975)), lty = 3, lwd = 3)
lines(Time, apply(Norm, 2, function(x) quantile(x, probs = 0.025)), lty = 3, lwd = 3)
dev.off()

pdf("/Users/sam/Desktop/VAE/6_norm2.pdf", height = 6, width = 8)
plot(Time, apply(Norm2, 2, mean), type = "l", ylim = c(0, 0.025), lwd = 3, ylab = "Gradient Norm2")
lines(Time, apply(Norm2, 2, function(x) quantile(x, probs = 0.975)), lty = 3, lwd = 3)
lines(Time, apply(Norm2, 2, function(x) quantile(x, probs = 0.025)), lty = 3, lwd = 3)
dev.off()


library(coda)
par(mfcol = c(1, 2))
traceplot(as.mcmc(reg$phi))


pred1 <- predict(reg, NewTimes)
Mean <- Lower <- Upper <- matrix(nrow = P, ncol = K)
MeanT <- LowerT <- UpperT <- matrix(nrow = P, ncol = K)
for (t in 1:P) {
  foot <- pred1$y[t, , ]
  fooT <- pred1$theta[t, , ]
  Mean[t, ] <- apply(foot, 1, mean)
  Lower[t, ] <- apply(foot, 1, function(x) quantile(x, probs = 0.025))
  Upper[t, ] <- apply(foot, 1, function(x) quantile(x, probs = 0.975))
  MeanT[t, ] <- apply(fooT, 1, mean)
  LowerT[t, ] <- apply(fooT, 1, function(x) quantile(x, probs = 0.025))
  UpperT[t, ] <- apply(fooT, 1, function(x) quantile(x, probs = 0.975))
}
pred2 <- predict(reg, NewTimes = (Time[T] + 1))
mean <- apply(pred2$y[1, , ], 1, mean)
upper <- apply(pred2$y[1, , ], 1, function(x) quantile(x, probs = 0.975))
lower <- apply(pred2$y[1, , ], 1, function(x) quantile(x, probs = 0.025))
meanT <- apply(pred2$theta[1, , ], 1, mean)
upperT <- apply(pred2$theta[1, , ], 1, function(x) quantile(x, probs = 0.975))
lowerT <- apply(pred2$theta[1, , ], 1, function(x) quantile(x, probs = 0.025))

pdf("/Users/sam/Desktop/VAE/1_latentspace.pdf", height = 6, width = 8)
par(mfcol = c(1, 1))
plot(Time, Y[1:T], ylim = c(min(Min[1], Min[2]), max(Max[1], Max[2])), ylab = "Latent Space", xlim = c(min(NewTimes), max(NewTimes)), pch = 16)
points(Time, Y[(T + 1):(2 * T)], col = 2, pch = 16)
lines(NewTimes, Mean[, 1], lwd = 3)
lines(NewTimes, Mean[, 2], col = 2, lwd = 3)
lines(NewTimes, Lower[, 1], lty = 3, lwd = 3)
lines(NewTimes, Lower[, 2], lty = 3, col = 2, lwd = 3)
lines(NewTimes, Upper[, 1], lty = 3, lwd = 3)
lines(NewTimes, Upper[, 2], lty = 3, col = 2, lwd = 3)
dev.off()

pdf("/Users/sam/Desktop/VAE/0_mean.pdf", height = 6, width = 8)
par(mfcol = c(1, 1))
plot(Time, Y[1:T], ylim = c(min(Min[1], Min[2]), max(Max[1], Max[2])), ylab = "Latent Space", xlim = c(min(NewTimes), max(NewTimes)), pch = 16)
points(Time, Y[(T + 1):(2 * T)], col = 2, pch = 16)
lines(NewTimes, MeanT[, 1], lwd = 3)
lines(NewTimes, MeanT[, 2], col = 2, lwd = 3)
lines(NewTimes, LowerT[, 1], lty = 3, lwd = 3)
lines(NewTimes, LowerT[, 2], lty = 3, col = 2, lwd = 3)
lines(NewTimes, UpperT[, 1], lty = 3, lwd = 3)
lines(NewTimes, UpperT[, 2], lty = 3, col = 2, lwd = 3)
dev.off()

pdf("/Users/sam/Desktop/VAE/3_2dls.pdf", height = 6, width = 8)
plot(Y[1:T], Y[(T + 1):(2 * T)], xlab = "X", ylab = "Y", ylim = c(Min[2], Max[2]), xlim = c(Min[1], Max[1]), pch = 16)
lines(Mean[1:P, 1], Mean[1:P, 2], col = "red", pch = 16, lwd = 3)
points(mean[1], mean[2], col = "green", pch = 16)
segments(x0 = mean[1], y0 = lower[2], x1 = mean[1], y1 = upper[2], col = "green", lty = 4)
segments(y0 = mean[2], x0 = lower[1], y1 = mean[2], x1 = upper[1], col = "green", lty = 4)
pred3 <- predict(reg, NewTimes = Time[T])
mean3 <- apply(pred3$y[1, , ], 1, mean)
points(mean3[1], mean3[2], col = "blue", pch = 16)
dev.off()

pdf("/Users/sam/Desktop/VAE/2_2dmean.pdf", height = 6, width = 8)
plot(Y[1:T], Y[(T + 1):(2 * T)], xlab = "X", ylab = "Y", ylim = c(Min[2], Max[2]), xlim = c(Min[1], Max[1]), pch = 16)
lines(MeanT[1:P, 1], MeanT[1:P, 2], col = "red", pch = 16, lwd = 3)
points(meanT[1], meanT[2], col = "green", pch = 16)
segments(x0 = meanT[1], y0 = lowerT[2], x1 = meanT[1], y1 = upperT[2], col = "green", lty = 4)
segments(y0 = meanT[2], x0 = lowerT[1], y1 = meanT[2], x1 = upperT[1], col = "green", lty = 4)
pred3 <- predict(reg, NewTimes = Time[T])
mean3T <- apply(pred3$theta[1, , ], 1, mean)
points(mean3T[1], mean3T[2], col = "blue", pch = 16)
dev.off()