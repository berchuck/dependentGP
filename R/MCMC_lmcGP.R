#' MCMC sampler for multivariate Gaussian process model using the linear model of coregionalization
#'
#' \code{lmcGP} is a Markov chain Monte Carlo (MCMC) sampler for a Gaussian process model
#'  using the Bayesian hierarchical framework.
#'
#' @param Y An \code{N} dimensional vector containing the observed outcome data.
#'  Here, \code{N = K * T}, where \code{T} represents the number of temporal observations
#'  and \code{K} the number of latent dimensions. The observations in \code{Y} must be first
#'  ordered temporally and then by latent dimension, meaning the first \code{T} observations
#'  in \code{Y} should come from the first latent dimension.
#'
#' @param Time A \code{T} dimensional vector containing the observed time points for each
#'  vector of outcomes in increasing order.
#'
#' @param Kernel Character string indicating the kernel of the Gaussian process. Options
#'  include: \code{"sqexp"}, \code{"exp"}, \code{"ar1"}.
#'
#' @param Starting Either \code{NULL} or a \code{list} containing starting values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the starting values may be specified.
#'
#'  When \code{NULL} is chosen then default starting values are automatically generated.
#'  Otherwise a \code{list} must be provided with names \code{Gamma}, \code{Sigma2}, \code{Phi}, or
#'  \code{T} containing appropriate objects. \code{Gamma} must be a \code{N} dimensional
#'  vector or a scalar that populates the entire vector, \code{Phi} and \code{Sigma2} \code{K} dimensional vectors, 
#'  or a scalar populating each entire vector, and \code{T} a \code{K x K} dimensional matrix.
#'
#' @param Hypers Either \code{NULL} or a \code{list} containing hyperparameter values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the hyperparameter values may be specified.
#'
#'  When \code{NULL} is chosen then default hyperparameter values are automatically
#'  generated. These default hyperparameters are described in detail in (Berchuck et al.).
#'  Otherwise a \code{list} must be provided with names \code{Phi}, \code{Sigma2} or
#'  \code{T} containing further hyperparameter information. These objects are themselves
#'  \code{lists} and may be constructed as follows.
#'
#'  \code{Sigma2} is a \code{list} with two objects, \code{Alpha} and \code{Beta}.
#'  \code{Alpha} represents the shape of the inverse gamma hyperprior and
#'  must be a scalar, while \code{Beta} represents the rate of the hyperprior
#'  and must be a scalar.
#'
#'  \code{T} is a \code{list} with two objects, \code{Xi} and \code{Psi}. \code{Xi}
#'  represents the degrees of freedom parameter for the inverse-Wishart hyperprior and
#'  must be a real number scalar, while \code{Psi} represents the scale matrix
#'  and must be a \code{K x K} dimensional positive definite matrix.
#'
#'  \code{Phi} is a \code{list} with two objects, \code{APhi} and \code{BPhi}. \code{APhi}
#'  represents the lower bound for the uniform hyperprior, while \code{BPhi} represents
#'  the upper bound. The bounds must be specified carefully. For example, if the exponential
#'  temporal correlation structure is chosen both bounds must be restricted to be non-negative. If a 
#'  scalar is give for these values they will populate the entire vector, otherwise these must be \code{K} 
#'  dimensional.
#'
#' @param Tuning Either \code{NULL} or a \code{list} containing tuning values
#'  to be specified for the MCMC Metropolis steps. If \code{NULL} is not chosen then all
#'  of the tuning values must be specified.
#'
#'  When \code{NULL} is chosen then default tuning values are automatically generated to
#'  \code{1}. Otherwise a \code{list} must be provided with names \code{Phi}, 
#'  and \code{T}. \code{Phi} must be a \code{K} dimensional vector and \code{T} a \code{((K * (K + 1)) / 2)}
#'  vector. Each containing tuning variances for their corresponding Metropolis updates.
#'
#' @param MCMC Either \code{NULL} or a \code{list} containing input values to be used
#'  for implementing the MCMC sampler. If \code{NULL} is not chosen then all
#'  of the MCMC input values must be specified.
#'
#'  \code{NBurn}: The number of sampler scans included in the burn-in phase. (default =
#'  \code{10,000})
#'
#'  \code{NSims}: The number of post-burn-in scans for which to perform the
#'   sampler. (default = \code{100,000})
#'
#'  \code{NThin}: Value such that during the post-burn-in phase, only every
#'  \code{NThin}-th scan is recorded for use in posterior inference (For return values
#'  we define, NKeep = NSims / NThin (default = \code{10}).
#'
#'  \code{NPilot}: The number of times during the burn-in phase that pilot adaptation
#'  is performed (default = \code{20})
#'
#' @param Seed An integer value used to set the seed for the random number generator
#'  (default = 54).
#'  
#' @param Verbose A logical indicating whether MCMC sampler progress should be printed.
#'
#' @details Details of the underlying statistical model can be found in the upcoming paper.
#'
#' @return \code{lmcGP} returns a list containing the following objects
#'
#'   \describe{
#'
#'   \item{\code{gamma}}{\code{NKeep x N} \code{matrix} of posterior samples for \code{gamma}.}
#'
#'   \item{\code{theta}}{\code{NKeep x N} \code{matrix} of posterior samples for \code{theta}.}
#'
#'   \item{\code{sigma2}}{\code{NKeep x K} \code{matrix} of posterior samples for \code{sigma2}.}
#'
#'   \item{\code{phi}}{\code{NKeep x K} \code{matrix} of posterior samples for \code{phi}.}
#'
#'   \item{\code{a}}{\code{NKeep x ((K * (K + 1)) / 2)} \code{matrix} of posterior samples for \code{A}. The
#'   columns have names that describe the samples within them. The row is listed first.}
#'   
#'   \item{\code{t}}{\code{NKeep x ((K * (K + 1)) / 2)} \code{matrix} of posterior samples for \code{T}. The
#'   columns have names that describe the samples within them. The row is listed first.}
#'
#'   \item{\code{metropolis}}{\code{(K + ((K * (K + 1)) / 2)) x 3} \code{matrix} of metropolis
#'   acceptance rates, tuners, a nd original tuners that result from the pilot adaptation. The first \code{K}
#'   correspond to the \code{phi} parameters.}
#'
#'   \item{\code{runtime}}{A \code{character} string giving the runtime of the MCMC sampler.}
#'
#'   \item{\code{datobj}}{A \code{list} of data objects that are used in future \code{lmcGP} functions
#'   and should be ignored by the user.}
#'
#'   }
#'
#' @author Samuel I. Berchuck
#' @references VAE paper forthcoming.
#' @export
lmcGP <- function(Y,
                  Time,
                  Kernel = "sqexp",
                  Hypers = NULL,
                  Starting = NULL,
                  Tuning = NULL,
                  MCMC = NULL,
                  Seed = 54,
                  Verbose = TRUE) {

  ###Function Inputs
  # Y = Y
  # Time = Time
  # Kernel = "sqexp"
  # Hypers = NULL
  # Starting = NULL
  # Tuning = NULL
  # MCMC = NULL
  # Seed = 54
  # Verbose = TRUE

  ###Check for missing objects
  if (missing(Y)) stop("Y: missing")
  if (missing(Time)) stop("Time: missing")

  ###Check model inputs
  CheckInputs(Y, Time, Kernel, Starting, Hypers, Tuning, MCMC, Seed, Verbose)

  ####Set seed for reproducibility
  set.seed(Seed)

  ###Check to see if the job is interactive
  Interactive <- interactive()

  ###Create objects for use in sampler
  DatObj <- CreateDatObj(Y, Time, Kernel)
  HyPara <- CreateHyPara(Hypers, DatObj)
  MetrObj <- CreateMetrObj(Tuning, DatObj)
  Para <- CreatePara(Starting, DatObj, HyPara)
  McmcObj <- CreateMcmc(MCMC, DatObj)
  RawSamples <- CreateStorage(DatObj, McmcObj)

  ###Time MCMC sampler
  BeginTime <- Sys.time()

  ###Run MCMC sampler in Rcpp
  RegObj <- lmcGP_Rcpp(DatObj, HyPara, MetrObj, Para, McmcObj, RawSamples, Interactive, Verbose)
  
  ###Set regression objects
  RawSamples <- RegObj$rawsamples
  MetropRcpp <- RegObj$metropolis

  ###End time
  FinishTime <- Sys.time()
  RunTime <- FinishTime - BeginTime

  ###Collect output to be returned
  DatObjOut <- OutputDatObj(DatObj)
  Metropolis <- SummarizeMetropolis(DatObj, MetrObj, MetropRcpp, McmcObj)
  Samples <- FormatSamples(DatObj, RawSamples)

  ###Return spBFA object
  lmcGP <- list(gamma = Samples$Gamma,
                theta = Samples$Theta,
                sigma2 = Samples$Sigma2,
                phi = Samples$Phi,
                a = Samples$A,
                t = Samples$T,
                metropolis = Metropolis,
                datobj = DatObjOut,
                runtime = paste0("Model runtime: ",round(RunTime, 2), " ", attr(RunTime, "units")))
  lmcGP <- structure(lmcGP, class = "lmcGP")
  return(lmcGP)

  ###End sampler
}
