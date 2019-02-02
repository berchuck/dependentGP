CheckInputs <- function(Y, Time, Kernel, Starting, Hypers, Tuning, MCMC, Seed, Verbose) {
  
  ###Data dimensions
  N <- length(as.numeric(Y))
  T <- length(as.numeric(Time))
  K <- N / T
  
  ###GP Kernel
  if (!Kernel %in% c("exp", "ar1", "sqexp")) stop('Kernel: must be one of "exp", "sqexp" or "ar1"')
  
  ###Seed
  if (!is.wholenumber(Seed)) stop('Seed must be an integer')
  if (is.na(Seed)) stop('Seed cannot be NA')
  if (!is.finite(Seed)) stop('Seed cannot be infinite')

  ###Verbose
  if (!is.logical(Verbose)) stop('Verbose must be a logical')

  ###Data checks for Y
  if (!is.matrix(Y)) stop('Y must be a vector')
  if (length(Y) != N) stop(paste0('Y must have exactly ', N, 'values'))
  if (any(is.na(Y))) stop("Y may have no missing values")
  if (any(!is.finite(Y))) stop("Y must have strictly finite entries")

  ###Data checks for Time
  if (!is.numeric(Time)) stop('Time must be a vector')
  if (length(Time) != T) stop(paste0('Time must have length ', T))
  if (any(is.na(Time))) stop("Time may have no missing values")
  if (any(!is.finite(Time))) stop("Time must have strictly finite entries")
  if (is.unsorted(Time)) stop('Time vector is not in increasing order')
  if (!all(Time >= 0)) stop('Time vector has at least one negative point')

  ###Hypers
  if (!is.null(Hypers)) {
    if (!is.list(Hypers)) stop('Hypers must be a list')
    if (!all(names(Hypers) %in% c("Phi", "Sigma2", "T"))) stop('Hypers: Can only contain lists with names "Phi", "Sigma2" or "T"')

    ###If Sigma2 hyperparameters are provided
    if ("Sigma2" %in% names(Hypers)) {
      if (!is.list(Hypers$Sigma2)) stop('Hypers: "Sigma2" must be a list')
      if (!"Alpha" %in% names(Hypers$Sigma2)) stop('Hypers: "Alpha" value missing')
      if (!is.scalar(Hypers$Sigma2$Alpha)) stop('Hypers: "Alpha" must be a scalar')
      if (is.na(Hypers$Sigma2$Alpha)) stop('Hypers: "Alpha" cannot be NA')
      if (!is.finite(Hypers$Sigma2$Alpha)) stop('Hypers: "Alpha" cannot be infinite')
      if (Hypers$Sigma2$Alpha <= 0) stop('Hypers: "Alpha" must be strictly positive')
      if (!"Beta" %in% names(Hypers$Sigma2)) stop('Hypers: "Beta" value missing')
      if (!is.scalar(Hypers$Sigma2$Beta)) stop('Hypers: "Beta" must be a scalar')
      if (is.na(Hypers$Sigma2$Beta)) stop('Hypers: "Beta" cannot be NA')
      if (!is.finite(Hypers$Sigma2$Beta)) stop('Hypers: "Beta" cannot be infinite')
      if (Hypers$Sigma2$Beta <= 0) stop('Hypers: "Beta" must be strictly positive')
    }
    
    ###If T hyperparameters are provided
    if ("T" %in% names(Hypers)) {
      if (!is.list(Hypers$T)) stop('Hypers: "T" must be a list')
      if (!"Xi" %in% names(Hypers$T)) stop('Hypers: "T" value missing for Xi')
      if (!is.scalar(Hypers$T$Xi)) stop('Hypers: "Xi" must be a scalar')
      if (is.na(Hypers$T$Xi)) stop('Hypers: "Xi" cannot be NA')
      if (!is.finite(Hypers$T$Xi)) stop('Hypers: "T" cannot be infinite')
      if (Hypers$T$Xi <= 0) stop('Hypers: "Xi" must be strictly positive')
      if (!"Psi" %in% names(Hypers$T)) stop('Hypers: "Psi" value missing')
      if (!is.matrix(Hypers$T$Psi)) stop('Hypers: "Psi" must be a matrix')
      if (!dim(Hypers$T$Psi)[1] == O) stop('Hypers: "Psi" must be O dimensional')
      if (!all(!is.na(Hypers$T$Psi))) stop('Hypers: "Psi" cannot have missing values')
      if (!all(is.finite(Hypers$T$Psi))) stop('Hypers: "Psi" cannot have infinite values')
      if (!dim(Hypers$T$Psi)[2] == O) stop('Hypers: "Psi" must be square')
      if (sum( !( (Hypers$T$Psi) == t(Hypers$T$Psi) ) ) > 0) stop('Hypers: "Psi" must be symmetric')
      if ((det(Hypers$T$Psi) - 0) < 0.00001) stop('Hypers: "Psi" is close to singular')
      
    }

    ###If Phi hyperparameters are provided
    if ("Phi" %in% names(Hypers)) {
      if (!is.list(Hypers$Psi)) stop('Hypers: "Phi" must be a list')
      if ((Kernel == "exp") | (Kernel == "sqexp")) {
        if (!"APhi" %in% names(Hypers$Phi)) stop('Hypers: "APhi" value missing')
        if (is.numeric(Hypers$Phi$APhi)) {
          if (length(Hypers$Phi$APhi) != K) stop('Hypers: "APhi" must be length K')
          if (!all(!is.na(Hypers$Phi$APhi))) stop('Hypers: "APhi" cannot have missing values')
          if (!all(is.finite(Hypers$Phi$APhi))) stop('Hypers: "APhi" cannot have infinite values')
          if (any(Hypers$Phi$APhi < 0)) stop('Hypers: All entries in "APhi" must be non-negative')
        }
        if (is.scalar(Hypers$Phi$APhi)) {
          if (is.na(Hypers$Phi$APhi)) stop('Hypers: "APhi" cannot be NA')
          if (!is.finite(Hypers$Phi$APhi)) stop('Hypers: "APhi" cannot be infinite')
          if (Hypers$Phi$APhi < 0) stop('Hypers: "APhi" must be non-negative')
        }
        if ((!is.numeric(Hypers$Phi$APhi)) | (!is.scalar(Hypers$Phi$APhi))) stop('Hypers: "APhi" must be a scalar or a vector')
        if (!"BPhi" %in% names(Hypers$Phi)) stop('Hypers: "BPhi" value missing')
        if (is.numeric(Hypers$Phi$BPhi)) {
          if (length(Hypers$Phi$BPhi) != K) stop('Hypers: "BPhi" must be length K')
          if (!all(!is.na(Hypers$Phi$BPhi))) stop('Hypers: "BPhi" cannot have missing values')
          if (!all(is.finite(Hypers$Phi$BPhi))) stop('Hypers: "BPhi" cannot have infinite values')
          if (any(Hypers$Phi$BPhi < 0)) stop('Hypers: All entries in "BPhi" must be non-negative')
          if (any(Hypers$Phi$BPhi < Hypers$Phi$APhi)) stop('Hypers: "BPhi" must be greater than "APhi"')
        }
        if (is.scalar(Hypers$Phi$BPhi)) {
          if (is.na(Hypers$Phi$BPhi)) stop('Hypers: "BPhi" cannot be NA')
          if (!is.finite(Hypers$Phi$BPhi)) stop('Hypers: "BPhi" cannot be infinite')
          if (Hypers$Phi$BPhi < 0) stop('Hypers: "BPhi" must be non-negative')
          if (Hypers$Phi$BPhi < Hypers$Phi$APhi) stop('Hypers: "BPhi" must be greater than "APhi"')
        }
        if ((!is.numeric(Hypers$Phi$BPhi)) | (!is.scalar(Hypers$Phi$BPhi))) stop('Hypers: "BPhi" must be a scalar or a vector')
      }
      if (Kernel == "ar1") {
        if (!"APhi" %in% names(Hypers$Phi)) stop('Hypers: "APhi" value missing')
        if (is.numeric(Hypers$Phi$APhi)) {
          if (length(Hypers$Phi$APhi) != K) stop('Hypers: "APhi" must be length K')
          if (!all(!is.na(Hypers$Phi$APhi))) stop('Hypers: "APhi" cannot have missing values')
          if (!all(is.finite(Hypers$Phi$APhi))) stop('Hypers: "APhi" cannot have infinite values')
          if (any(Hypers$Phi$APhi < -1) | any(Hypers$Phi$APhi > 1)) stop('Hypers: "APhi" must be in the range (-1, 1)')
        }
        if (is.scalar(Hypers$Phi$APhi)) {
          if (is.na(Hypers$Phi$APhi)) stop('Hypers: "APhi" cannot be NA')
          if (!is.finite(Hypers$Phi$APhi)) stop('Hypers: "APhi" cannot be infinite')
          if (Hypers$Phi$APhi < 0) stop('Hypers: "APhi" must be non-negative')
          if ((Hypers$Phi$APhi < -1) | (Hypers$Phi$APhi > 1)) stop('Hypers: "APhi" must be in the range (-1, 1)')
        }
        if ((!is.numeric(Hypers$Phi$APhi)) | (!is.scalar(Hypers$Phi$APhi))) stop('Hypers: "APhi" must be a scalar or a vector')
        if (!"BPhi" %in% names(Hypers$Phi)) stop('Hypers: "BPhi" value missing')
        if (is.numeric(Hypers$Phi$BPhi)) {
          if (length(Hypers$Phi$BPhi) != K) stop('Hypers: "BPhi" must be length K')
          if (!all(!is.na(Hypers$Phi$BPhi))) stop('Hypers: "BPhi" cannot have missing values')
          if (!all(is.finite(Hypers$Phi$BPhi))) stop('Hypers: "BPhi" cannot have infinite values')
          if (any(Hypers$Phi$BPhi < 0)) stop('Hypers: All entries in "BPhi" must be non-negative')
          if (any(Hypers$Phi$BPhi < -1) | any(Hypers$Phi$BPhi > 1)) stop('Hypers: "BPhi" must be in the range (-1, 1)')
          if (any(Hypers$Phi$BPhi < Hypers$Phi$APhi)) stop('Hypers: "BPhi" must be greater than "APhi"')
        }
        if (is.scalar(Hypers$Phi$BPhi)) {
          if (is.na(Hypers$Phi$BPhi)) stop('Hypers: "BPhi" cannot be NA')
          if (!is.finite(Hypers$Phi$BPhi)) stop('Hypers: "BPhi" cannot be infinite')
          if (Hypers$Phi$BPhi < 0) stop('Hypers: "BPhi" must be non-negative')
          if ((Hypers$Phi$BPhi < -1) | (Hypers$Phi$BPhi > 1)) stop('Hypers: "BPhi" must be in the range (-1, 1)')
          if (Hypers$Phi$BPhi < Hypers$Phi$APhi) stop('Hypers: "BPhi" must be greater than "APhi"')
        }
        if ((!is.numeric(Hypers$Phi$BPhi)) | (!is.scalar(Hypers$Phi$BPhi))) stop('Hypers: "BPhi" must be a scalar or a vector')
      }
    }

  ###End Hyperparameters
  }

  ###Starting Values
  if (!is.null(Starting)) {
    if (!is.list(Starting)) stop('Starting must be a list')
    if (!all(names(Starting) %in% c("Sigma2", "Gamma", "Phi", "T"))) stop('Starting: Can only contain objects with names "Sigma2", "Gamma", "Phi", and "T"')

    ###If Sigma2 starting values is provided
    if ("Sigma2" %in% names(Starting)) {
      if (is.numeric(Starting$Sigma2)) {
        if (!is.numeric(Starting$Sigma2)) stop('Starting: "Sigma2" must be a vector')
        if (length(Starting$Sigma2) != K) stop('Starting: "Sigma2" must be length K')
        if (!all(!is.na(Starting$Sigma2))) stop('Starting: "Sigma2" cannot have missing values')
        if (!all(is.finite(Starting$Sigma2))) stop('Starting: "Sigma2" cannot have infinite values')
        if (any(Starting$Sigma2 <= 0)) stop('Starting: Entries in "Sigma2" must be strictly positive')
      }
      if (is.scalar(Starting$Sigma2)) {
        if (is.na(Starting$Sigma2)) stop('Starting: "Sigma2" cannot be NA')
        if (!is.finite(Starting$Sigma2)) stop('Starting: "Sigma2" cannot be infinite')
        if (Starting$Sigma2 <= 0) stop('Starting: "Sigma2" must be strictly positive')
      }
      if ((!is.numeric(Starting$Sigma2)) | (!is.scalar(Starting$Sigma2))) stop('Starting: "Sigma2" must be a scalar or a vector')
    }
    
    ###If Gamma starting values is provided
    if ("Gamma" %in% names(Starting)) {
      if (is.numeric(Starting$Gamma)) {
        if (!is.numeric(Starting$Gamma)) stop('Starting: "Gamma" must be a vector')
        if (length(Starting$Gamma) != N) stop('Starting: "Gamma" must be length N')
        if (!all(!is.na(Starting$Gamma))) stop('Starting: "Gamma" cannot have missing values')
        if (!all(is.finite(Starting$Gamma))) stop('Starting: "Gamma" cannot have infinite values')
      }
      if (is.scalar(Starting$Gamma)) {
        if (is.na(Starting$Gamma)) stop('Starting: "Gamma" cannot be NA')
        if (!is.finite(Starting$Gamma)) stop('Starting: "Gamma" cannot be infinite')
      }
      if ((!is.numeric(Starting$Gamma)) | (!is.scalar(Starting$Gamma))) stop('Starting: "Gamma" must be a scalar or a vector')
    }
    
    ###If T starting values is provided
    if ("T" %in% names(Starting)) {
      if (!is.matrix(Starting$T)) stop('Starting: "T" must be a matrix')
      if (!dim(Starting$T)[1] == K) stop('Starting: "T" must be K dimensional')
      if (!dim(Starting$T)[2] == K) stop('Starting: "T" must be square')
      if (!all(!is.na(Starting$T))) stop('Starting: "T" cannot have missing values')
      if (!all(is.finite(Starting$T))) stop('Starting: "T" cannot have infinite values')
      if (sum( !( (Starting$T) == t(Starting$T) ) ) > 0) stop('Starting: "T" must be symmetric')
      if ((det(Starting$T) - 0) < 0.0000000001) stop('Starting: "T" is close to singular')
    }

    ###If Phi starting values is provided
    if ("Phi" %in% names(Starting)) {
      if (Kernel == "exp" | Kernel == "sqexp") {
        if (is.numeric(Starting$Phi)) {
          if (length(Starting$Phi) != K) stop('Starting: "Phi" must be length K')
          if (any(is.na(Starting$Phi))) stop('Starting: "Phi" cannot contain NA')
          if (any(!is.finite(Starting$Phi))) stop('Starting: "Phi" cannot be infinite')
          if (any(Starting$Phi <= 0)) stop('Starting: "Phi" must be non-negative')
        }
        if (is.scalar(Starting$Phi)) {
          if (is.na(Starting$Phi)) stop('Starting: "Phi" cannot be NA')
          if (!is.finite(Starting$Phi)) stop('Starting: "Phi" cannot be infinite')
          if (Starting$Phi <= 0) stop('Starting: "Phi" must be non-negative')
        }
        if ((!is.numeric(Starting$Phi)) | (!is.scalar(Starting$Phi))) stop('Starting: "Phi" must be a scalar or a vector')
      }
      if (Kernel == "ar1") {
        if (is.numeric(Starting$Phi)) {
          if (length(Starting$Phi) != K) stop('Starting: "Phi" must be length K')
          if (any(is.na(Starting$Phi))) stop('Starting: "Phi" cannot contain NA')
          if (any(!is.finite(Starting$Phi))) stop('Starting: "Phi" cannot be infinite')
          if (any(Starting$Phi <= -1) | any(Starting$Rho >= 1)) stop('Starting: "Phi" must be in (-1, 1)')
        }
        if (is.scalar(Starting$Phi)) {
          if (is.na(Starting$Phi)) stop('Starting: "Phi" cannot be NA')
          if (!is.finite(Starting$Phi)) stop('Starting: "Phi" cannot be infinite')
          if ((Starting$Phi <= -1) | (Starting$Rho >= 1)) stop('Starting: "Phi" must be in (-1, 1)')
        }
        if ((!is.numeric(Starting$Phi)) | (!is.scalar(Starting$Phi))) stop('Starting: "Phi" must be a scalar or a vector')        # I make sure that Rho is in (ARho, BRho) in CreatePara();
      }
    }
    # I make sure Phi is in the lower and upper bound in CreatePara()
    
  ###End Starting Values
  }

  ###Tuning Values
  if (!is.null(Tuning)) {
    if (!is.list(Tuning)) stop('Tuning must be a list')
    if (!all(names(Tuning) %in% c("Phi", "T"))) stop('Tuning: Can only contain objects with names "Phi", "T"')

    ###If Phi tuning value is provided
    if ("Phi" %in% names(Tuning)) {
      if (is.numeric(Tuning$Phi)) {
        if (length(Tuning$Phi) != K) stop('Tuning: "Phi" must be length K')
        if (any(is.na(Tuning$Phi))) stop('Tuning: "Phi" cannot be NA')
        if (any(!is.finite(Tuning$Phi))) stop('Tuning: "Phi" cannot be infinite')
        if (any(Tuning$Phi < 0)) stop('Tuning: "Phi" must be non-negative')
      }
      if (is.scalar(Tuning$Phi)) {
        if (is.na(Tuning$Phi)) stop('Tuning: "Phi" cannot be NA')
        if (!is.finite(Tuning$Phi)) stop('Tuning: "Phi" cannot be infinite')
        if (Tuning$Phi < 0) stop('Tuning: "Phi" must be non-negative')
      }
      if ((!is.numeric(Tuning$Phi)) | (!is.scalar(Tuning$Phi))) stop('Tuning: "Phi" must be a scalar or a vector')
    }

    ###If T tuning value is provided
    if ("T" %in% names(Tuning)) {
      if (is.numeric(Tuning$T)) {
        if (length(Tuning$T) != ((K * (K + 1)) / 2)) stop('Tuning: "T" must be length ((K * (K + 1)) / 2)')
        if (any(is.na(Tuning$T))) stop('Tuning: "T" cannot be NA')
        if (any(!is.finite(Tuning$T))) stop('Tuning: "T" cannot be infinite')
        if (any(Tuning$T < 0)) stop('Tuning: "T" must be non-negative')
      }
      if (is.scalar(Tuning$T)) {
        if (is.na(Tuning$T)) stop('Tuning: "T" cannot be NA')
        if (!is.finite(Tuning$T)) stop('Tuning: "T" cannot be infinite')
        if (Tuning$T < 0) stop('Tuning: "T" must be non-negative')
      }
      if ((!is.numeric(Tuning$T)) | (!is.scalar(Tuning$T))) stop('Tuning: "T" must be a scalar or a vector')
    }

  ###End Tuning Values
  }

  ###MCMC Values
  if (!is.null(MCMC)) {
    if (!is.list(MCMC)) stop('MCMC must be a list')
    if (!all(names(MCMC) %in% c("NBurn", "NSims", "NThin", "NPilot"))) stop('MCMC: Can only contain objects with names "NBurn", "NSims", "NThin" and "NPilot"')

    ###If NBurn is provided
    if ("NBurn" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NBurn)) stop('MCMC: "NBurn" must be a scalar')
      if (is.na(MCMC$NBurn)) stop('MCMC: "NBurn" cannot be NA')
      if (!is.finite(MCMC$NBurn)) stop('MCMC: "NBurn" cannot be infinite')
      if (!is.wholenumber(MCMC$NBurn) | MCMC$NBurn < 0) stop('MCMC: "NBurn" must be a non-negative integer')
      if (MCMC$NBurn < 100) stop('MCMC: "NBurn" must be at least 100')
    }

    ###If NSims is provided
    if ("NSims" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NSims)) stop('MCMC: "NSims" must be a scalar')
      if (is.na(MCMC$NSims)) stop('MCMC: "NSims" cannot be NA')
      if (!is.finite(MCMC$NSims)) stop('MCMC: "NSims" cannot be infinite')
      if (!is.wholenumber(MCMC$NSims) | MCMC$NSims <= 0) stop('MCMC: "NSims" must be a positive integer')
      if (MCMC$NSims < 100) stop('MCMC: "NSims" must be at least 100')
    }

    ###If NThin is provided
    if ("NThin" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NThin)) stop('MCMC: "NThin" must be a scalar')
      if (is.na(MCMC$NThin)) stop('MCMC: "NThin" cannot be NA')
      if (!is.finite(MCMC$NThin)) stop('MCMC: "NThin" cannot be infinite')
      if (!is.wholenumber(MCMC$NThin) | MCMC$NThin <= 0) stop('MCMC: "NThin" must be a positive integer')
      # if (!is.wholenumber(MCMC$NSims / MCMC$NThin)) stop('MCMC: "NThin" must be a factor of "NSims"') enforced in CreateMCMC();
    }

    ###If NPilot is provided
    if ("NPilot" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NPilot)) stop('MCMC: "NPilot" must be a scalar')
      if (is.na(MCMC$NPilot)) stop('MCMC: "NPilot" cannot be NA')
      if (!is.finite(MCMC$NPilot)) stop('MCMC: "NPilot" cannot be infinite')
      if (!is.wholenumber(MCMC$NPilot) | MCMC$NPilot < 0) stop('MCMC: "NPilot" must be a positive integer')
      # if (!is.wholenumber(MCMC$NBurn / MCMC$NPilot)) stop('MCMC: "NPilot" must be a factor of "NBurn"') enforced in CreateMCMC();
    }

  ###End MCMC Values
  }
 
###End CheckInputs function 
}

###Helper Functions
is.scalar <- function(x) ((is.numeric(x)) & (length(x) == 1))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
