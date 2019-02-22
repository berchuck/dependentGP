###Function for reading in sampler inputs and creating a list object that contains all relavent data objects--------------------
CreateDatObj <- function(Y, Time, Kernel) {

  ###Data objects
  Y <- matrix(Y, ncol = 1)
  T <- length(Time)
  N <- dim(Y)[1]
  K <- N / T
  YWide <- matrix(Y, nrow = T, ncol = K)

  ###Matrix objects
  EyeT <- diag(T)
  EyeK <- diag(K)
  EyeTK <- diag(T * K)
  ZeroT <- matrix(0, nrow = T, ncol = 1)

  ###Gaussian process kernel
  if (Kernel == "exp") {
    KernInd <- 0
    TimeDist <- abs(outer(Time, Time, "-"))
  }
  if (Kernel == "sqexp") {
    KernInd <- 1
    TimeDist <- 0.5 * (outer(Time, Time, "-"))^2
  }
  if (Kernel == "ar1") {
    KernInd <- 2
    TimeDist <- abs(outer(Time, Time, "-"))
  }

  ###A Indices
  AInd <- matrix(nrow = K, ncol = K)
  AInd[upper.tri(AInd, diag = TRUE)] <- 1:(((K + 1) * K) / 2)
  AInd <- t(AInd)
  AInd <- which(!is.na(AInd), arr.ind = TRUE)
  AInd <- AInd[order(AInd[, 1]), ]

  ###Make parameters global
  DatObj <- list()
  DatObj$Y <- Y
  DatObj$YWide <- YWide
  DatObj$T <- T
  DatObj$N <- N
  DatObj$K <- K
  DatObj$EyeT <- EyeT
  DatObj$EyeK <- EyeK
  DatObj$EyeTK <- EyeTK
  DatObj$ZeroT <- ZeroT
  DatObj$TimeDist <- TimeDist
  DatObj$Time <- Time
  DatObj$AInd <- matrix(AInd - 1, nrow = (K * (K + 1)) / 2, ncol = 2)
  DatObj$KernInd <- KernInd
  return(DatObj)

}



###Function for creating an object containing relevant Hyperparameter information---------------------------------------------------
CreateHyPara <- function(Hypers, DatObj) {

  ###Set data objects
  K <- DatObj$K
  TimeDist <- DatObj$TimeDist
  KernInd <-  DatObj$KernInd

  ###Which parameters are user defined?
  UserHypers <- names(Hypers)

  ###Set hyperparameters for Phi
  if (KernInd == 0) { # exponential
    if ("Phi" %in% UserHypers) {
      APhi <- matrix(Hypers$Phi$APhi, nrow = K, ncol = 1)
      BPhi <- matrix(Hypers$Phi$BPhi, nrow = K, ncol = 1)
    }
    if (!("Phi" %in% UserHypers)) {
      minDiff <- min(TimeDist[TimeDist > 0])
      maxDiff <- max(TimeDist[TimeDist > 0])
      Lower <- -log(0.50) / maxDiff #longest diff goes up to 50%
      Upper <- -log(0.25) / minDiff #shortest diff goes down to 25%
      APhi <- matrix(Lower, nrow = K, ncol = 1)
      BPhi <- matrix(Upper, nrow = K, ncol = 1)
    }
  }
  if (KernInd == 1) { # squared exponential
    if ("Phi" %in% UserHypers) {
      APhi <- matrix(Hypers$Phi$APhi, nrow = K, ncol = 1)
      BPhi <- matrix(Hypers$Phi$BPhi, nrow = K, ncol = 1)
    }
    if (!("Phi" %in% UserHypers)) {
      minDiff <- min(TimeDist[TimeDist > 0])
      maxDiff <- max(TimeDist[TimeDist > 0])
      A <- sqrt(-maxDiff / log(0.95)) #longest diff goes down to 25%
      B <- sqrt(-minDiff / log(0.999)) #shortest diff goes up to 99.9%
      Lower <- min(A, B)
      Upper <- max(A, B)
      APhi <- matrix(Lower, nrow = K, ncol = 1)
      BPhi <- matrix(Upper, nrow = K, ncol = 1)
    }
  }
  if (KernInd == 2) { # ar1
    if ("Phi" %in% UserHypers) {
      APhi <- matrix(Hypers$Phi$APhi, nrow = K, ncol = 1)
      BPhi <- matrix(Hypers$Phi$BPhi, nrow = K, ncol = 1)
    }
    if (!("Phi" %in% UserHypers)) {
      minDiff <- min(TimeDist[TimeDist > 0])
      maxDiff <- max(TimeDist[TimeDist > 0])
      Lower <- 0.5 ^ (-1 / maxDiff)  #longest diff goes up to 50%
      Upper <- 0.25 ^ (-1 / minDiff) #shortest diff goes down to 1%
      APhi <- matrix(Lower, nrow = K, ncol = 1)
      BPhi <- matrix(Upper, nrow = K, ncol = 1)
    }
  }

  ###Set hyperparameters for Sigma2
  if ("Sigma2" %in% UserHypers) {
    Alpha <- Hypers$Sigma2$Alpha
    Beta <- Hypers$Sigma2$Beta
  }
  if (!("Sigma2" %in% UserHypers)) {
    Alpha <- 0.001
    Beta <- 0.001
  }

  ###Set hyperparameters for T
  if ("T" %in% UserHypers) {
    Xi <- Hypers$T$Xi
    Psi <- Hypers$T$Psi
  }
  if (!("T" %in% UserHypers)) {
    Xi <- K + 1
    Psi <- diag(K)
  }

  ###Create object for hyperparameters
  HyPara <- list()
  HyPara$Alpha <- Alpha
  HyPara$Beta <- Beta
  HyPara$APhi <- APhi
  HyPara$BPhi <- BPhi
  HyPara$Xi <- Xi
  HyPara$Psi <- Psi
  return(HyPara)

}



###Function for creating an object containing relevant Metropolis information---------------------------------------------------
CreateMetrObj <- function(Tuning, DatObj) {

  ###Set data objects
  K <- DatObj$K

  ###Which parameters are user defined?
  UserTuners <- names(Tuning)

  ###Set tuning parameters for Phi
  if ("Phi" %in% UserTuners) MetropPhi <- matrix(Tuning$Phi, nrow = K, ncol = 1)
  if (!("Phi" %in% UserTuners)) MetropPhi <- matrix(1, nrow = K, ncol = 1)

  ###Set tuning parameters for T
  if ("T" %in% UserTuners) MetropT <- matrix(Tuning$T, nrow = ((K + 1) * K) / 2, ncol = 1)
  if (!("T" %in% UserTuners)) MetropT <- matrix(1, nrow = ((K + 1) * K) / 2, ncol = 1)

  ###Set acceptance rate counters
  AcceptanceT <- matrix(0, nrow = ((K + 1) * K) / 2, ncol = 1)
  AcceptancePhi <- matrix(0, nrow = K, ncol = 1)

  ###Return metropolis object
  MetrObj <- list()
  MetrObj$MetropT <- MetropT
  MetrObj$AcceptanceT <- AcceptanceT
  MetrObj$MetropPhi <- MetropPhi
  MetrObj$AcceptancePhi <- AcceptancePhi
  MetrObj$OriginalTuners <- c(MetropPhi, MetropT)
  return(MetrObj)

}



###Function for creating inital parameter object-------------------------------------------------------------------------------
CreatePara <- function(Starting, DatObj, HyPara) {

  ###Set data objects
  K <- DatObj$K
  T <- DatObj$T
  EyeT <- DatObj$EyeT
  EyeK <- DatObj$EyeK
  TimeDist <- DatObj$TimeDist
  KernInd <- DatObj$KernInd

  ###Set hyperparameter objects
  APhi <- HyPara$APhi
  BPhi <- HyPara$BPhi

  ###Which parameters are user defined?
  UserStarters <- names(Starting)

  ###Set initial values of Gamma
  if ("Gamma" %in% UserStarters) GammaWide <- matrix(Starting$Gamma, nrow = K, ncol = T)
  if ((!"Gamma" %in% UserStarters)) GammaWide <- matrix(0, nrow = K, ncol = T)

  ###Set initial values of Sigma2
  if ("Sigma2" %in% UserStarters) Sigma2 <- matrix(Starting$Sigma2, nrow = K, ncol = 1)
  if ((!"Sigma2" %in% UserStarters)) Sigma2 <- matrix(1, nrow = K, ncol = 1)

  ###Set initial value of Phi
  if ("Phi" %in% UserStarters) {
    Phi <- matrix(Starting$Phi, nrow = K, ncol = 1)
    if (any(Phi <= APhi) | any(Phi >= BPhi)) stop('Starting: "Phi" must be in the interval (APhi, BPhi)')
  }
  if (!("Phi" %in% UserStarters)) Phi <-  matrix((APhi + BPhi) / 2, nrow = K, ncol = 1)

  ###Set initial value of T
  if ("T" %in% UserStarters) TMat <- matrix(Starting$T, nrow = K, ncol = K)
  if (!("T" %in% UserStarters)) TMat <- diag(K)

  ###Set Gaussian process parameters
  Gamma <- matrix(t(GammaWide), ncol = 1)
  A <- t(chol(TMat))
  HPhi <- HPhiInv <- CholHPhi <- array(dim = c(T, T, K))
  for (k in 1:K) {
    HPhi[ , , k] <- H(Phi[k], TimeDist, KernInd, T)
    CholHPhi[ , , k] <- chol(HPhi[ , , k])
    HPhiInv[ , , k] <- chol2inv(CholHPhi[ , , k])
  }

  ###Covariance
  Sigma <- diag(as.numeric(Sigma2))
  SigmaInv <- diag(as.numeric(1 / Sigma2))
  CholSigma <- diag(as.numeric(sqrt(Sigma2)))

  ###T
  logDet <- determinant(TMat, logarithm = TRUE)
  logDetT <- as.numeric(logDet$modulus * logDet$sign)
  AInv <- forwardsolve(A, EyeK)
  TInv <- t(AInv) %*% AInv

  ###Save parameter objects
  Para <- list()
  Para$Gamma <- Gamma
  Para$GammaWide <- GammaWide
  Para$Theta <- kronecker(A, EyeT) %*% Gamma
  Para$Sigma2 <- Sigma2
  Para$Sigma <- Sigma
  Para$SigmaInv <- SigmaInv
  Para$CholSigma <- CholSigma
  Para$Phi <- Phi
  Para$HPhi <- HPhi
  Para$HPhiInv <- HPhiInv
  Para$CholHPhi  <- CholHPhi
  Para$TMat <- TMat
  Para$logDetT <- logDetT
  Para$TInv <- TInv
  Para$A <- A
  return(Para)

}



###Function that creates inputs for MCMC sampler--------------------------------------------------------------------------------
CreateMcmc <- function(MCMC, DatObj) {

  ###Which parameters are user defined?
  UserMCMC <- names(MCMC)

  ###Set MCMC objects
  if ("NBurn" %in% UserMCMC) NBurn <- MCMC$NBurn
  if (!("NBurn" %in% UserMCMC)) NBurn <- 5000
  if ("NSims" %in% UserMCMC) NSims <- MCMC$NSims
  if (!("NSims" %in% UserMCMC)) NSims <- 5000
  if ("NThin" %in% UserMCMC) NThin <- MCMC$NThin
  if (!("NThin" %in% UserMCMC)) NThin <- 1
  if ("NPilot" %in% UserMCMC) NPilot <- MCMC$NPilot
  if (!("NPilot" %in% UserMCMC)) NPilot <- 10

  ###One last check of MCMC user inputs
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!(is.wholenumber(NSims / NThin))) stop('MCMC: "NThin" must be a factor of "NSims"')
  if (!(is.wholenumber(NBurn / NPilot))) stop('MCMC: "NPilot" must be a factor of "NBurn"')

  ###Create MCMC objects
  NTotal <- NBurn + NSims
  WhichKeep <- NBurn + (1:(NSims / NThin)) * NThin
  NKeep <- length(WhichKeep)

  ###Pilot adaptation objects
  WhichPilotAdapt <- (1:NPilot) * NBurn / NPilot
  PilotAdaptDenominator <- WhichPilotAdapt[1]

  ###Burn-in progres bar
  BarLength <- 50 #Burn-in bar length (arbitrary)
  BurnInProgress <- seq(1 / BarLength, 1, 1 / BarLength)
  WhichBurnInProgress <- sapply(BurnInProgress, function(x) tail(which(1 : NBurn <= x * NBurn), 1))

  ###Progress output objects
  SamplerProgress <- seq(0.1, 1.0, 0.1) #Intervals of progress update (arbitrary)
  WhichSamplerProgress <- sapply(SamplerProgress, function(x) tail(which(1:NSims <= x * NSims), 1)) + NBurn
  WhichBurnInProgressInt <- sapply(SamplerProgress, function(x) tail(which(1:NBurn <= x * NBurn), 1))

  ###Save objects
  MCMC <- list()
  MCMC$NBurn <- NBurn
  MCMC$NSims <- NSims
  MCMC$NThin <- NThin
  MCMC$NPilot <- NPilot
  MCMC$NTotal <- NTotal
  MCMC$WhichKeep <- WhichKeep
  MCMC$NKeep <- NKeep
  MCMC$WhichPilotAdapt <- WhichPilotAdapt
  MCMC$PilotAdaptDenominator <- PilotAdaptDenominator
  MCMC$BurnInProgress <- BurnInProgress
  MCMC$WhichBurnInProgress <- WhichBurnInProgress
  MCMC$WhichBurnInProgressInt <- WhichBurnInProgressInt
  MCMC$BarLength <- BarLength
  MCMC$WhichSamplerProgress <- WhichSamplerProgress
  return(MCMC)

}



###Function that creates a storage object for raw samples-----------------------------------------------------------------------
CreateStorage <- function(DatObj, McmcObj) {

  ###Set data objects
  K <- DatObj$K
  N <- DatObj$N

  ###Set MCMC objects
  NKeep <- McmcObj$NKeep

  ###Create storage object
  # Gamma: N
  # Theta: N
  # Sigma2: K
  # Phi: K
  # A: ((K + 1) * K) / 2
  # T: ((K + 1) * K) / 2
  Out <- matrix(nrow = (N + N + K + K + ((K + 1) * K) / 2 + ((K + 1) * K) / 2), ncol = NKeep)
  return(Out)

}
