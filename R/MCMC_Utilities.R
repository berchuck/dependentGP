###Function for summarizing the raw MCMC samples-------------------------------------------------------------------
FormatSamples <- function(DatObj, RawSamples) {

  ###Set data objects
  K <- DatObj$K
  N <- DatObj$N
  T <- DatObj$T
  AInd <- DatObj$AInd

  ###Format raw samples
  RawSamples <- t(RawSamples)
  Gamma <- RawSamples[, 1:N]
  Theta <- RawSamples[, (N + 1):(N + N)]
  Sigma2 <- RawSamples[, (N + N + 1):(N + N + K)]
  Phi <- RawSamples[, (N + N + K + 1):(N + N + K + K)]
  A <- RawSamples[, (N + N + K + K + 1):(N + N + K + K + ((K + 1) * K) / 2)]
  TMat <- RawSamples[, (N + N + K + K + ((K + 1) * K) / 2 + 1):(N + N + K + K + ((K + 1) * K) / 2 + ((K + 1) * K) / 2)]
  colnames(Gamma) <- paste0("Gamma", rep(1:K, each = T), "_", rep(1:T, K))
  colnames(Theta) <- paste0("Theta", rep(1:K, each = T), "_", rep(1:T, K))
  colnames(Sigma2) <- paste0("Sigma2_", 1:K)
  colnames(Phi) <- paste0("Phi", 1:K)
  colnames(A) <- paste0("A",  AInd[, 1], "_", AInd[, 2])
  colnames(TMat) <- paste0("T",  AInd[, 1], "_", AInd[, 2])
  Out <- list(Gamma = Gamma, Theta = Theta, Sigma2 = Sigma2, Phi = Phi, A = A, T = TMat)
  return(Out)
}



###Function for creating a data object that contains objects needed for ModelFit-----------------------------------
OutputDatObj <- function(DatObj) {

  ###Collect needed objects
  # DatObjOut <- list(M = DatObj$M,
  #                   Nu = DatObj$Nu,
  #                   AdjacentEdgesBoolean = DatObj$AdjacentEdgesBoolean,
  #                   W = DatObj$W,
  #                   EyeM = DatObj$EyeM,
  #                   EyeNu = DatObj$EyeNu,
  #                   OneM = DatObj$OneM,
  #                   OneN = DatObj$OneN,
  #                   OneNu = DatObj$OneNu,
  #                   YStarWide = DatObj$YStarWide,
  #                   Rho = DatObj$Rho,
  #                   FamilyInd = DatObj$FamilyInd,
  #                   ScaleY = DatObj$ScaleY,
  #                   YObserved = DatObj$YObserved,
  #                   ScaleDM = DatObj$ScaleDM,
  #                   Time = DatObj$Time,
  #                   TimeVec = DatObj$TimeVec,
  #                   YObserved = DatObj$YObserved,
  #                   tNu = DatObj$tNu,
  #                   t1 = DatObj$t1,
  #                   XThetaInd = DatObj$XThetaInd,
  #                   N = DatObj$N,
  #                   EyeN = DatObj$EyeN)
  DatObjOut <- DatObj
  return(DatObjOut)

}



###Function for summarizing Metropolis objects post sampler--------------------------------------------------------
SummarizeMetropolis <- function(DatObj, MetrObj, MetropRcpp, McmcObj) {

  ###Set data object
  K <- DatObj$K

  ###Set MCMC object
  NSims <- McmcObj$NSims

  ###Set Metropolis objects
  MetropPhi <- MetropRcpp$MetropPhi
  AcceptancePhi <- MetropRcpp$AcceptancePhi
  MetropT <- MetropRcpp$MetropT
  AcceptanceT <- MetropRcpp$AcceptanceT
  OriginalTuners <- MetrObj$OriginalTuners
  AcceptancePct <- c(AcceptancePhi, AcceptanceT) / NSims
  MetrSummary <- cbind(AcceptancePct, c(MetropPhi, MetropT), OriginalTuners)
  SigmaInd <- which(lower.tri(apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("A", 1:K, "_"), x)), diag = TRUE), arr.ind = TRUE)
  ALabs <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("A", 1:K, "_"), x))[SigmaInd[order(SigmaInd[, 1]), ]]
  rownames(MetrSummary) <- c(paste0("Phi", 1:K), ALabs)
  colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")

  ###Summarize and output
  return(MetrSummary)

}



###Verify the class of our regression object------------------------------------------------------------------------
#' is.lmcGP
#'
#' \code{is.lmcGP} is a general test of an object being interpretable as a
#' \code{\link{lmcGP}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{lmcGP}} class is defined as the regression object that
#'  results from the \code{\link{lmcGP}} regression function.
#'
#' @export
is.lmcGP <- function(x) {
  identical(attributes(x)$class, "lmcGP")
}
