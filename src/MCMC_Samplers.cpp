#define ARMA_DONT_PRINT_ERRORS //So the cholesky warning is suppressed
#include <RcppArmadillo.h>
#include "MCMC_lmcGP.h"

//Function to sample Sigma2 using a Gibbs sampler step---------------------------------------------------------------
para SampleSigma2(datobj DatObj, para Para, hypara HyPara) {

  //Set data objects
  arma::mat YWide = DatObj.YWide;
  int K = DatObj.K;
  int T = DatObj.T;
  
  //Set parameter objects
  arma::colvec Theta = Para.Theta;
  arma::mat ThetaWide = arma::reshape(Theta, T, K);
  
  //Set hyperparameter objects
  double Alpha = HyPara.Alpha;
  double Beta = HyPara.Beta;

  //Shape is constant over locations
  double Shape = Alpha + 0.5 * T;

  //Declarations
  double Rate;
  arma::colvec Sigma2(K), ResidsK(K);
  
  //Loop over latent dimensions
  for (arma::uword k = 0; k < K; k++) {
    
  
    //Calculate rate
    ResidsK = YWide.col(k) - ThetaWide.col(k);
    Rate = Beta + 0.5 * arma::as_scalar(arma::trans(ResidsK) * ResidsK);
    
    //Sample Sigma2k
    Sigma2(k) = rigammaRcpp(Shape, Rate);
    
  //End loop over latent dimensions
  }

  //Update dependent parameters
  arma::mat Sigma = arma::diagmat(Sigma2);
  arma::mat SigmaInv = arma::diagmat(1 / Sigma2);
  arma::mat CholSigma = arma::diagmat(sqrt(Sigma2));

  //Update parameters object
  Para.Sigma2 = Sigma2;
  Para.Sigma = Sigma;
  Para.SigmaInv = SigmaInv;
  Para.CholSigma = CholSigma;
  return Para;

}



//Function to sample Gamma using a Gibbs sampler step---------------------------------------------------------------
para SampleGamma(datobj DatObj, para Para) {

  //Set data objects
  arma::colvec Y = DatObj.Y;
  int K = DatObj.K;
  int T = DatObj.T;
  arma::mat EyeT = DatObj.EyeT;
  arma::uvec Seq = arma::linspace<arma::uvec>(0, K - 1, K);
  
  //Set parameters
  arma::mat A = Para.A;
  arma::mat SigmaInv = Para.SigmaInv;
  arma::cube HPhiInv = Para.HPhiInv;
  arma::mat GammaWide = Para.GammaWide;
  
  //Loop over latent dimensions
  for (arma::uword k = 0; k < K; k++){
    
    //Calculate moments
    arma::colvec AK = A.col(k);
    arma::uvec Index = find(Seq != k);
    arma::mat AMinusK = A.cols(Index);
    arma::mat GammaMinusK = GammaWide.rows(Index);    
    arma::colvec YK = Y - arma::kron(AMinusK, EyeT) * arma::trans(GammaMinusK);
    arma::rowvec tAKSigma = arma::trans(AK) * SigmaInv;
    arma::mat CovGammaK = CholInv(arma::as_scalar(tAKSigma * AK) * EyeT + HPhiInv.slice(k));
    arma::colvec MeanGammaK = CovGammaK * (arma::kron(tAKSigma, EyeT) * YK);
    
    //Sample Gammak
    GammaWide.row(k) = arma::trans(rmvnormRcpp(1, MeanGammaK, CovGammaK));
    
  //End loop over latent dimensions
  }
  
  //Update parameters dependent on Gamma
  arma::colvec Gamma = arma::reshape(arma::trans(GammaWide), T * K, 1);
  arma::colvec Theta = arma::kron(A, EyeT) * Gamma;

  //Update parameters object
  Para.GammaWide = GammaWide;
  Para.Gamma = Gamma;
  Para.Theta = Theta;
  return Para;

}



//Function to sample new value of Phi using a Metropolis sampler step-----------------------------------------------
std::pair<para, metrobj> SamplePhi(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj) {

  //Set data objects
  int K = DatObj.K;
  int T = DatObj.T;
  arma::mat TimeDist = DatObj.TimeDist;
  arma::mat EyeT = DatObj.EyeT;
  arma::colvec ZeroT = DatObj.ZeroT;
  int KernInd = DatObj.KernInd;
    
  //Set parameter objects
  arma::colvec Phi = Para.Phi;
  arma::cube HPhi = Para.HPhi;
  arma::cube HPhiInv = Para.HPhiInv;
  arma::cube CholHPhi = Para.CholHPhi;
  arma::mat GammaWide = Para.GammaWide;
      
  //Set hyperparameters
  arma::colvec APhi = HyPara.APhi;
  arma::colvec BPhi = HyPara.BPhi;
        
  //Set metropolis objects
  arma::colvec MetropPhi = sqrt(MetrObj.MetropPhi);
  arma::Col<int> AcceptancePhi = MetrObj.AcceptancePhi;
  
  //Loop over latent dimension
  for (arma::uword k = 0; k < K; k++) {
    
    //Dimension specific objects
    double PhiK = Phi(k);
    double APhiK = APhi(k);
    double BPhiK = BPhi(k);
    arma::mat CholHPhiK = CholHPhi.slice(k);
    arma::colvec GammaK = arma::trans(GammaWide.row(k));
      
    //Transform current state to real line
    double BigDelta = log((PhiK - APhiK) / (BPhiK - PhiK));
        
    //Numerical fix for when the propopsal cholesky doesn't exist
    double PhiKProposal, BigDeltaProposal;
    arma::mat CholHPhiKProposal(T, T), HPhiKProposal(T, T);
    bool Cholesky = false;
    while (!Cholesky) {
        
      //Sample a new Proposal
      BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropPhi(k)));

      //Compute PhiK Proposal
      PhiKProposal = (BPhiK * exp(BigDeltaProposal) + APhiK) / (1 + exp(BigDeltaProposal));
                
      //Proposal temporal correlation
      HPhiKProposal = H(PhiKProposal, TimeDist, KernInd, T);
      Cholesky = arma::chol(CholHPhiKProposal, HPhiKProposal);

    }
          
    //ThetaI structure components
    arma::mat RootiGammaK = arma::solve(arma::trimatu(CholHPhiK), EyeT);
    arma::mat RootiGammaKProposal = arma::solve(arma::trimatu(CholHPhiKProposal), EyeT);
    double Component1A = lndMvn(GammaK, ZeroT, RootiGammaKProposal);
    double Component1B = lndMvn(GammaK, ZeroT, RootiGammaK);
    double Component1 = Component1A - Component1B;
            
    //Prior components
    double Component2A = BigDeltaProposal;
    double Component2B = BigDelta;
    double Component2 = Component2A - Component2B;
              
    //Jacobian component 1
    double Component3 = 2 * log((1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)));
                
    //Compute log acceptance ratio
    double LogR = Component1 + Component2 + Component3;
                  
    //Metropolis update
    double RandU = randuRcpp();
    if (log(RandU) < LogR) {
                    
      //Keep Count of Acceptances
      AcceptancePhi(k)++;
              
      //Update dependent parameters
      arma::mat P = arma::solve(arma::trimatu(CholHPhiKProposal), EyeT);
              
      //Update parameters object
      Phi(k) = PhiKProposal;
      HPhi.slice(k) = HPhiKProposal;
      CholHPhi.slice(k) = CholHPhiKProposal;
      HPhiInv.slice(k) = P * arma::trans(P);
                      
    }
                  
  //End loop over latent dimension
  }
  
  //Update Metropolis object
  MetrObj.AcceptancePhi = AcceptancePhi;
    
  //Update Parameter object
  Para.Phi = Phi;
  Para.HPhi = HPhi;
  Para.CholHPhi = CholHPhi;
  Para.HPhiInv = HPhiInv;
      
  //Return output object
  return std::pair<para, metrobj>(Para, MetrObj);
  
}



//Function to sample new value of Phi using a Metropolis sampler step-----------------------------------------------
std::pair<para, metrobj> SampleA(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj) {
  
  //Set data objects
  arma::colvec Y = DatObj.Y;
  int K = DatObj.K;
  arma::mat EyeT = DatObj.EyeT;
  arma::mat EyeK = DatObj.EyeK;
  arma::mat EyeTK = DatObj.EyeTK;
  arma::umat AInd = DatObj.AInd;
  
  //Set parameter objects
  arma::colvec Gamma = Para.Gamma;
  arma::mat A = Para.A;
  arma::mat CholSigma = Para.CholSigma;
  arma::mat RootiY = arma::solve(arma::trimatu(arma::kron(CholSigma, EyeT)), EyeTK);
  arma::mat TMat = Para.TMat;
  double logDetT = Para.logDetT;
  arma::mat TInv = Para.TInv;
  arma::colvec Theta = Para.Theta;
  
  //Set hyperparameters
  double Xi = HyPara.Xi;
  arma::mat Psi = HyPara.Psi;
  
  //Set metropolis objects
  arma::colvec MetropT = sqrt(MetrObj.MetropT);
  arma::Col<int> AcceptanceT = MetrObj.AcceptanceT;

  //Loop over components of A
  for (arma::uword k = 0; k < (((K + 1) * K) / 2); k++) {
    
    //Diagonal or off?
    arma::uvec Indeces = arma::trans(AInd.row(k));
    arma::uword Row = Indeces(0);
    arma::uword Col = Indeces(1);
    bool Diag = (Row == Col);
    
    //Dimension specific objects
    double AComp = A(Row, Col);
    arma::mat AProposal = A;
      
    //Transform current state to real line
    double BigDelta;
    if (Diag) BigDelta = log(AComp);
    if (!Diag) BigDelta = AComp;
    
    //Sample a new Proposal
    double BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropT(k)));
            
    //Compute proposal
    double ACompProposal;
    if (Diag) ACompProposal = exp(BigDeltaProposal);
    if (!Diag) ACompProposal = BigDeltaProposal;
                
    //Proposal values
    AProposal(Row, Col) = ACompProposal;
    arma::mat TMatProposal = AProposal * arma::trans(AProposal);
    arma::colvec ThetaProposal = arma::kron(AProposal, EyeT) * Gamma;
    double val, sign;
    arma::log_det(val, sign, TMatProposal);
    double logDetTProposal = val + log(sign);
    arma::mat AInvProposal = arma::solve(arma::trimatl(AProposal), EyeK);
    arma::mat TInvProposal = arma::trans(AInvProposal) * AInvProposal;
    
    //Likelihood components
    double Component1A = lndMvn(Y, ThetaProposal, RootiY);
    double Component1B = lndMvn(Y, Theta, RootiY);
    double Component1 = Component1A - Component1B;
                  
    //Prior components for T
    double Component2A = (-(Xi + K + 1) / 2) * logDetTProposal - 0.5 * arma::trace(Psi * TInvProposal);
    double Component2B = (-(Xi + K + 1) / 2) * logDetT - 0.5 * arma::trace(Psi * TInv);
    double Component2 = Component2A - Component2B;
                    
    //Jacobian component
    double JacobianProposal = 0;
    double Jacobian = 0;
    for (arma::uword j = 0; j < K; j++) {
      JacobianProposal += (K + 1 - j) * log(AProposal(j, j));
      Jacobian += (K + 1 - j) * log(A(j, j));
    }
    double Component3A = JacobianProposal;
    double Component3B = Jacobian;
    double Component3 = Component3A - Component3B;
                      
    //Transformation Jacobian
    double Component4 = 0;
    if (Diag) Component4 += (BigDeltaProposal - BigDelta);
                        
    //Compute log acceptance ratio
    double LogR = Component1 + Component2 + Component3 + Component4;
                          
    //Metropolis update
    double RandU = randuRcpp();
    if (log(RandU) < LogR) {
                            
      //Keep Count of Acceptances
      AcceptanceT(k)++;
      
      //Update parameters object
      A = AProposal;
      TMat = TMatProposal;
      logDetT = logDetTProposal;
      TInv = TInvProposal;
      Theta = ThetaProposal;
                              
    }
                          
  //End loop over components of A
  }
  
  //Update Metropolis object
  MetrObj.AcceptanceT = AcceptanceT;
    
  //Update Parameter object
  Para.A = A;
  Para.TMat = TMat;
  Para.logDetT = logDetT;
  Para.TInv = TInv;
  Para.Theta = Theta;
  
  //Return output object
  return std::pair<para, metrobj>(Para, MetrObj);
  
}