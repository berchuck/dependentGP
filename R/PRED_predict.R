###Function to obtain predictions for lmcGP object
#'
#' predict.lmcGP
#'
#' Predicts future observations from the \code{\link{lmcGP}} model.
#'
#' @param object a \code{\link{lmcGP}} model object for which predictions
#'  are desired from.
#'
#' @param NewTimes a numeric vector including desired time(s) points for prediction.
#'  
#' @param type a character string indicating the type of predictions to be returned. Either, 
#'  "response" or "derivative".
#'
#' @param Verbose A logical indicating whether MCMC sampler progress should be printed.
#'
#' @param ... other arguments.
#'
#' @details \code{predict.lmcGP} uses Bayesian krigging to predict the various processes.
#'
#' @return \code{predict.lmcGP} returns a \code{list} for type = "response" and a \code{matrix} for
#'  type = "derivative".
#'
#'   \describe{
#'
#'   \item{type = "response"}{For type = "response", the list contains the predictive samples for gamma, theta, 
#'   and the observed data y. Each of these are arrays with three dimensions. The first dimension
#'   is the number of time points to be predicted, the second the number of latent dimensions, and
#'   finally, the third dimensions is the number of posterior samples specified in the original model.}
#'
#'   \item{type = "derivative"}{An array of the derivative samples. The array has three dimensions. The first dimension
#'   is the number of time points to be predicted, the second the number of latent dimensions, and
#'   finally, the third dimensions is the number of posterior samples specified in the original model.}
#'
#'   }
#'
#' @author Samuel I. Berchuck
#' @export
###Prediction function for lmcGP function
predict.lmcGP <- function(object, NewTimes, type = "response", Verbose = TRUE, ...) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!is.lmcGP(object)) stop('"object" must be of class lmcGP')
  if (missing(NewTimes)) stop('"NewTimes" is missing')
  if (!is.numeric(NewTimes)) stop('NewTimes must be a vector')
  if (any(is.na(NewTimes))) stop("NewTimes may have no missing values")
  if (any(!is.finite(NewTimes))) stop("NewTimes must have strictly finite entries")
  if (!all(NewTimes >= 0)) stop('NewTimes vector has at least one negative entry')
  if (!is.character(type)) stop('response must be a character string')
  if (!type %in% c("response", "derivative")) stop('response must be one of "response" or "derivative"')

  ###Set seed for reproducibility
  set.seed(54)

  ###Set data objects
  DatObj <- object$datobj
  K <- DatObj$K
  T <- DatObj$T
  N <- DatObj$N
  AInd <- DatObj$AInd
  KernInd <- DatObj$KernInd

  ###Update A index
  Mat <- matrix(1:(K ^ 2), nrow = K, ncol = K)
  NewAInd <- NULL
  for (k in 1:(dim(AInd)[1])) {
    NewAInd <- c(NewAInd, Mat[AInd[k, 1] + 1, AInd[k, 2] + 1])
  }
  NewAInd <- NewAInd - 1
  
  ###Response prediction
  if (type == "response") {

    ###Create updated distance matrix
    TimeFixed <- DatObj$Time
    WhichOriginalReturnLong <- which(NewTimes %in% TimeFixed)
    WhichOriginalReturn <- which(TimeFixed %in% NewTimes)
    WhichNewReturn <- which(!(NewTimes %in% TimeFixed))
    NNewVisitsReturn <- length(WhichNewReturn)
    NOriginalReturn <- length(WhichOriginalReturn)
    NVisitsReturn <- sum(NNewVisitsReturn + NOriginalReturn)
    Time <- sort(c(TimeFixed, NewTimes))
    Time <- unique(Time)
    WhichOriginal <- which(Time %in% TimeFixed)
    WhichNew <- which(!Time %in% TimeFixed)
    NNewVisits <- length(WhichNew)
    NOriginal <- length(WhichOriginal)
    NVisits <- sum(NNewVisits + NOriginal)
    if (KernInd == 0) TimeDist <- abs(outer(Time, Time, "-"))
    if (KernInd == 1) TimeDist <- 0.5 * (outer(Time, Time, "-"))^2
    if (KernInd == 2) TimeDist <- abs(outer(Time, Time, "-"))
    NewVisits <- WhichNew
    OriginalVisits <- WhichOriginal

    ###Set mcmc object
    NKeep <- dim(object$phi)[1]

    ###Set parameter object
    A <- object$a
    Gamma <- object$gamma
    Sigma2 <- object$sigma2
    Phi <- object$phi

    ###Predict using C++ function
    out <- PredY(NKeep, NVisitsReturn, K, Gamma, T, Phi, 
      TimeDist, KernInd, OriginalVisits - 1, NewVisits - 1,
      NNewVisitsReturn, NOriginalReturn, WhichOriginalReturnLong - 1,
      WhichOriginalReturn - 1, WhichNewReturn - 1, A,
      NewAInd, Sigma2, NVisits, Verbose)
    
    ###Return formated samples
    return(out)

  ###End response prediction
  }

  ###Response prediction
  if (type == "derivative") {

    ###Data objects
    TimeDist <- DatObj$TimeDist
    NNewVisits <- length(NewTimes)
    TimesSort <- sort(NewTimes)
    EyeT <- DatObj$EyeT
    Y <- DatObj$Y

    ###Set mcmc object
    NKeep <- dim(object$phi)[1]

    ###Set parameter object
    A <- object$a
    Sigma2 <- object$sigma2
    Phi <- object$phi

    ###Covariance indeces
    IndecesOrg <- as.matrix(expand.grid(1:T, 1:K))
    IndecesNew <- as.matrix(expand.grid(1:NNewVisits, 1:K))

    ###Predict using C++ function
    Derivative <- PredDerv(NKeep, NNewVisits, K, Phi, Sigma2, A, 
                        NewAInd, N, IndecesOrg - 1, Time, IndecesNew - 1,
                        Y, TimesSort, Verbose)
    
    ###Return formated samples
    return(Derivative$derivative)

  ###End response prediction
  }

}
