#include <RcppArmadillo.h>
#include "MCMC_lmcGP.h"

//Prediction of new observations-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List PredY(int NKeep, int NVisitsReturn, int K, arma::mat Gamma, int T, arma::mat Phi, 
                 arma::mat TimeDist, int KernInd, arma::uvec OriginalVisits, arma::uvec NewVisits,
                 int NNewVisitsReturn, int NOriginalReturn, arma::uvec WhichOriginalReturnLong,
                 arma::uvec WhichOriginalReturn, arma::uvec WhichNewReturn, arma::mat A,
                 arma::uvec AInd, arma::mat Sigma2, int NVisits, bool Verbose) {

  //Verbose output
  arma::vec VerboseSeq;
  VerboseSeq << 0.10 << 0.20 << 0.3 << 0.4 << 0.50 << 0.6 << 0.70 << 0.80 << 0.90;
  VerboseSeq *= NKeep;
  if (Verbose) Rcpp::Rcout << std::fixed << "Krigging: 0%.. ";
  
  //Declarations
  // arma::cube GammaKrigged(NVisitsReturn, NKeep, K), ThetaKrigged(NVisitsReturn, NKeep, K), YKrigged(NVisitsReturn, NKeep, K); 
  arma::cube GammaKrigged(NVisitsReturn, K, NKeep), ThetaKrigged(NVisitsReturn, K, NKeep), YKrigged(NVisitsReturn, K, NKeep); 
  
  //Loop over scans
  for (arma::uword s = 0; s < NKeep; s++) {
    
    //Scan specific objects
    arma::mat GammaWideS = arma::reshape(Gamma.row(s), T, K);
    arma::colvec PhiS = arma::trans(Phi.row(s));
    
    //Loop over latent dimensions
    for (arma::uword k = 0; k < K; k++) {
      
      //Dimension specific objects
      arma::colvec GammaK = GammaWideS.col(k);
      double PhiK = PhiS(k);
      arma::mat HPhiK = H(PhiK, TimeDist, KernInd, NVisits);
      
      //Block covariance matrices
      arma::mat HPhiKOrgOrg = HPhiK(OriginalVisits, OriginalVisits);
      arma::mat HPhiKOrgOrgInv = CholInv(HPhiKOrgOrg);
      arma::mat NewToOriginal = HPhiK(NewVisits, OriginalVisits) * HPhiKOrgOrgInv;
      arma::mat NewToOriginalToNew = NewToOriginal * HPhiK(OriginalVisits, NewVisits);
      arma::mat NewToNew = HPhiK(NewVisits, NewVisits);
      
      //Moment calculations
      arma::mat CovK = NewToNew - NewToOriginalToNew;
      arma::colvec MeanK = NewToOriginal * GammaK;
      
      //Sample posterior predictive
      arma::colvec NewSamp;
      // if (NNewVisitsReturn > 0) NewSamp = arma::trans(arma::chol(CovK)) * rnormSNRcpp(NNewVisitsReturn) + MeanK;
      if (NNewVisitsReturn > 0) NewSamp = rmvnormRcppRobust(MeanK, CovK);
      
      //Store samples
      arma::colvec OriginalSamp = GammaK;
      arma::colvec gammaKrigged(NVisitsReturn);
      if (NOriginalReturn > 0) gammaKrigged(WhichOriginalReturnLong) = OriginalSamp(WhichOriginalReturn);
      if (NNewVisitsReturn > 0) gammaKrigged(WhichNewReturn) = NewSamp;
      GammaKrigged.slice(s).col(k) = gammaKrigged;
    
    //End loop over latent dimensions
    }
    
    //Add a new percentage
    if (Verbose) Rcpp::Rcout.precision(0);
    if (Verbose) if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
      Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
  }
  
  //Output final percentage
  if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;
  
  //Obtain predictive distributions of theta and y
  for (arma::uword s = 0; s < NKeep; s++) {
    arma::mat GammaS = arma::trans(GammaKrigged.slice(s));
    arma::rowvec AVec = A.row(s);
    arma::mat AS(K, K, arma::fill::zeros);
    AS(AInd) = AVec;
    arma::mat ThetaS = AS * GammaS;
    ThetaKrigged.slice(s) = arma::trans(ThetaS);
    arma::colvec Sigma2S = arma::trans(Sigma2.row(s));
    arma::mat YS(K, NVisitsReturn);
    for (arma::uword t = 0; t < NVisitsReturn; t++) YS.col(t) = rnormVecRcpp(ThetaS.col(t), arma::sqrt(Sigma2S));
    YKrigged.slice(s) = arma::trans(YS);
  }
  
  //Return list object
  return Rcpp::List::create(Rcpp::Named("gamma") = GammaKrigged,
                            Rcpp::Named("theta") = ThetaKrigged,
                            Rcpp::Named("y") = YKrigged);
}



//Prediction of derivatives-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List PredDerv(int NKeep, int NNewVisits, int K, arma::mat Phi, arma::mat Sigma2, arma::mat A, 
                    arma::uvec AInd, int N, arma::mat IndecesOrg, arma::colvec Time, arma::mat IndecesNew,
                    arma::colvec Y, arma::colvec TimesSort, bool Verbose) {
  
  //Verbose output
  arma::vec VerboseSeq;
  VerboseSeq << 0.10 << 0.20 << 0.3 << 0.4 << 0.50 << 0.60 << 0.70 << 0.80 << 0.90;
  VerboseSeq *= NKeep;
  if (Verbose) Rcpp::Rcout << std::fixed << "Krigging: 0%.. ";
  
  //Declarations
  arma::cube Derivative(NNewVisits, K, NKeep);
  arma::mat Sigma11(N, N), Sigma22(NNewVisits * K, NNewVisits * K), Sigma21(NNewVisits * K, N), AS(K, K);
  arma::mat SigmaStar((NNewVisits * K), N), Cov((NNewVisits * K), (NNewVisits * K));
  arma::colvec Mean(NNewVisits * K), NewSamp(NNewVisits * K);
  arma::rowvec Sigma2S(K), PhiS(K), AVec(K);
  arma::uword RowT, RowK, ColT, ColK, s, i, j, p;
  int KStar;
  double Sum, PhiSp, TimeRowT, TimeColT, TimeDiff, TimeSortRowT, TimeSortColT, TimeSortDiff;

  //Loop over scans
  for (s = 0; s < NKeep; s++) {

    //Scan objects
    PhiS = Phi.row(s);
    Sigma2S = Sigma2.row(s);
    AVec = A.row(s);
    AS.zeros();
    AS(AInd) = AVec;
    
    //Covariance for original outcomes
    for (i = 0; i < N; i++) {
      for (j = 0; j <= i; j++)  {
        RowT = IndecesOrg(i, 0);
        RowK = IndecesOrg(i, 1);
        ColT = IndecesOrg(j, 0);
        ColK = IndecesOrg(j, 1);
        KStar = std::min(ColK, RowK);
        TimeRowT = Time(RowT);
        TimeColT = Time(ColT);
        TimeDiff = (TimeRowT - TimeColT);
        Sum = 0;
        for (p = 0; p <= KStar; p++) {
          PhiSp = PhiS(p);
          Sum += AS(RowK, p) * AS(ColK, p) * exp(-0.5 * (1 / (PhiSp * PhiSp)) * TimeDiff * TimeDiff);
        }
        if (i == j) Sum += Sigma2S(RowK);
        Sigma11(i, j) = Sum;
      }
    }
    Sigma11 = arma::symmatl(Sigma11);

    //Covariance for the new derivatives
    for (i = 0; i < (NNewVisits * K); i++) {
      for (j = 0; j <= i; j++)  {
        RowT = IndecesNew(i, 0);
        RowK = IndecesNew(i, 1);
        ColT = IndecesNew(j, 0);
        ColK = IndecesNew(j, 1);
        KStar = std::min(ColK, RowK);
        TimeSortRowT = TimesSort(RowT);
        TimeSortColT = TimesSort(ColT);
        TimeSortDiff = (TimeSortRowT - TimeSortColT);
        Sum = 0;
        for (p = 0; p <= KStar; p++) {
          PhiSp = PhiS(p);
          Sum += AS(RowK, p) * AS(ColK, p) * (1 / (PhiSp * PhiSp * PhiSp * PhiSp)) * (PhiSp * PhiSp - TimeSortDiff * TimeSortDiff) * exp(-0.5 * (1 / (PhiSp * PhiSp)) * TimeSortDiff * TimeSortDiff);
        }
        Sigma22(i, j) = Sum;
      }
    }
    Sigma22 = arma::symmatl(Sigma22);
    
    //Off-diagonal covariance
    for (i = 0; i < (NNewVisits * K); i++) {
      RowT = IndecesNew(i, 0);
      RowK = IndecesNew(i, 1);
      for (j = 0; j < N; j++) {
        ColT = IndecesOrg(j, 0);
        ColK = IndecesOrg(j, 1);
        KStar = std::min(ColK, RowK);
        TimeSortRowT = TimesSort(RowT);
        TimeColT = Time(ColT);
        TimeDiff = (TimeSortRowT - TimeColT);
        Sum = 0;
        for (p = 0; p <= KStar; p++) {
          PhiSp = PhiS(p);
          Sum += (-1) * AS(RowK, p) * AS(ColK, p) * (1 / (PhiSp * PhiSp)) * TimeDiff * exp(-0.5 * (1 / (PhiSp * PhiSp)) * TimeDiff * TimeDiff);
        }
        Sigma21(i, j) = Sum;
      }
    }
    
    //Calculate moments
    SigmaStar = Sigma21 * CholInv(Sigma11);
    Cov = Sigma22 - SigmaStar * arma::trans(Sigma21);
    Mean = SigmaStar * Y;
    
    //Sample derivative
    NewSamp = rmvnormRcppRobust(Mean, Cov);
    Derivative.slice(s) = arma::reshape(NewSamp, NNewVisits, K);

    //Reset matrices
    Sigma11.zeros();
    Sigma21.zeros();
    Sigma22.zeros();
    
    //Add a new percentage
    if (Verbose) Rcpp::Rcout.precision(0);
    if (Verbose) if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
      Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
    
  //End loop over scans
  }
  
  //Output final percentage
  if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;
  
  //Return list object
  return Rcpp::List::create(Rcpp::Named("derivative") = Derivative);

//End function
}


