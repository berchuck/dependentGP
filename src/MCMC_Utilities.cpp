#include <RcppArmadillo.h>
#include "MCMC_lmcGP.h"

//Initiate burn-in progress bar--------------------------------------------------------------------------
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive) {

  //Set MCMC object
  int BarLength = McmcObj.BarLength;

  //Initialize burn-in bar
  if (Interactive) {
    Rcpp::Rcout << std::fixed << "Burn-in progress:  |";
    for (int i = 0; i < BarLength - 1; i++) Rcpp::Rcout << std::fixed << " ";
    Rcpp::Rcout << std::fixed <<  "|" << std::fixed;
  }
  if (!Interactive) {
    Rcpp::Rcout << std::fixed << "Burn-in progress:  0%..  ";
  }

}



//Function to pilot adapt tuning parameter--------------------------------------------------------------
double PilotAdaptFunc(double TuningParameter, double AcceptancePct) {

  //Adjust tuning parameter using scaling based on size of acceptance rate
  if (AcceptancePct >= 0.90) TuningParameter *= 1.3;
  if ( (AcceptancePct >= 0.75 ) & (AcceptancePct < 0.90 ) ) TuningParameter *= 1.2;
  if ( (AcceptancePct >= 0.45 ) & (AcceptancePct < 0.75 ) ) TuningParameter *= 1.1;
  if ( (AcceptancePct <= 0.25 ) & (AcceptancePct > 0.15 ) ) TuningParameter *= 0.9;
  if ( (AcceptancePct <= 0.15 ) & (AcceptancePct > 0.10 ) ) TuningParameter *= 0.8;
  if (AcceptancePct <= 0.10) TuningParameter *= 0.7;
  return TuningParameter;

}



//Function for implementing pilot adaptation in MCMC sampler--------------------------------------------
metrobj PilotAdaptation(metrobj MetrObj, mcmcobj McmcObj, datobj DatObj) {

  //Set Data objects
  int K = DatObj.K;
  
  //Set Metropolis objects
  arma::colvec MetropPhi = MetrObj.MetropPhi;
  arma::Col<int> AcceptancePhi = MetrObj.AcceptancePhi;
  arma::colvec MetropT = MetrObj.MetropT;
  arma::Col<int> AcceptanceT = MetrObj.AcceptanceT;
  
  //Set MCMC objects
  int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;

  //Update Phi
  double PctPhi, PctT, Update, Current;
  for (arma::uword k = 0; k < K; k++) {
    Current = MetropPhi(k);
    PctPhi = AcceptancePhi(k) / double(PilotAdaptDenominator);
    Update = PilotAdaptFunc(Current, PctPhi);
    MetropPhi(k) = Update;
    AcceptancePhi(k) = 0;
  }
  
  //Update T
  for (arma::uword k = 0; k < (((K + 1) * K) / 2); k++) {
    Current = MetropT(k);
    PctT = AcceptanceT(k) / double(PilotAdaptDenominator);
    Update = PilotAdaptFunc(Current, PctT);
    MetropT(k) = Update;
    AcceptanceT(k) = 0;
  }
  
  //Update Metropolis object
  MetrObj.MetropPhi = MetropPhi;
  MetrObj.AcceptancePhi = AcceptancePhi;
  MetrObj.MetropT = MetropT;
  MetrObj.AcceptanceT = AcceptanceT;
  return MetrObj;

}



//Output Metropolis object for summary-------------------------------------------------------------------
Rcpp::List OutputMetrObj(metrobj MetrObj, datobj DatObj) {

  Rcpp::List Out;
  Out = Rcpp::List::create(Rcpp::Named("AcceptanceT") = MetrObj.AcceptanceT,
                           Rcpp::Named("MetropT") = MetrObj.MetropT,
                           Rcpp::Named("AcceptancePhi") = MetrObj.AcceptancePhi,
                           Rcpp::Named("MetropPhi") = MetrObj.MetropPhi,
                           Rcpp::Named("OriginalTuners") = MetrObj.OriginalTuners);
  return Out;

}



//Initiate burn-in progress bar-------------------------------------------------------------------------------------
void SamplerProgress(int s, mcmcobj McmcObj) {

  //Set MCMC object
  int NSims = McmcObj.NSims;
  int NBurn = McmcObj.NBurn;

  //Add a new percentage
  Rcpp::Rcout.precision(0);
  if (s < NSims + NBurn)Rcpp::Rcout << std::fixed << 100 * (s - NBurn) / NSims << "%.. ";
  if (s == NSims + NBurn) Rcpp::Rcout << std::fixed << 100 * (s - NBurn) / NSims << "%!";

}



//Function for storing raw MCMC samples to to an object in memory
arma::colvec StoreSamples(datobj DatObj, para Para) {

  //Set data object
  int N = DatObj.N;
  int K = DatObj.K;

  //Set parameter objects
  arma::colvec Gamma = Para.Gamma;
  arma::colvec Theta = Para.Theta;
  arma::colvec Sigma2 = Para.Sigma2;
  arma::colvec Phi = Para.Phi;
  arma::mat A = Para.A;
  arma::mat TMat = Para.TMat;
  
  //Save raw samples
  int counter = 0;
  arma::colvec col(N + N + K + K + ((K + 1) * K) / 2 + ((K + 1) * K) / 2);
  for (arma::uword i = 0; i < N; i++) {
    col(counter) = Gamma(i);
    counter++;
  }
  for (arma::uword i = 0; i < N; i++) {
    col(counter) = Theta(i);
    counter++;
  }
  for (arma::uword i = 0; i < K; i++) {
    col(counter) = Sigma2(i);
    counter++;
  }
  for (arma::uword i = 0; i < K; i++) {
    col(counter) = Phi(i);
    counter++;
  }
  for (arma::uword i = 0; i < K; i++) {
    for (arma::uword j = 0; j <= i; j++) {
      col(counter) = A(i, j);
      counter++;
    }
  }
  for (arma::uword i = 0; i < K; i++) {
    for (arma::uword j = 0; j <= i; j++) {
      col(counter) = TMat(i, j);
      counter++;
    }
  }
  return col;
}



// //Update burn-in progress bar----------------------------------------------------------------------------
// void UpdateBurnInBar(int s, mcmcobj McmcObj) {
//
//   //Set MCMC object
//   arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
//   int BarLength = McmcObj.BarLength;
//
//   //Add a new star
//   arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
//   arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
//   int NewStar = NewStarBooleanVec(0);
//   for (int i = 0; i < (BarLength + 1 - NewStar); i++) Rcpp::Rcout << std::fixed << "\b";
//   Rcpp::Rcout << std::fixed << "*";
//   for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
//   Rcpp::Rcout << std::fixed << "|";
//
// }




//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateBurnInBarInt(int s, mcmcobj McmcObj) {

  //Set MCMC object
  arma::vec WhichBurnInProgressInt = McmcObj.WhichBurnInProgressInt;
  arma::uvec NewStarBoolean = find(s == WhichBurnInProgressInt);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);

  //Add percentage to submited job mode
  Rcpp::Rcout.precision(0);
  if (NewStar == 0) Rcpp::Rcout << std::fixed << "10%.. ";
  if (NewStar == 1) Rcpp::Rcout << std::fixed << "20%.. ";
  if (NewStar == 2) Rcpp::Rcout << std::fixed << "30%.. ";
  if (NewStar == 3) Rcpp::Rcout << std::fixed << "40%.. ";
  if (NewStar == 4) Rcpp::Rcout << std::fixed << "50%.. ";
  if (NewStar == 5) Rcpp::Rcout << std::fixed << "60%.. ";
  if (NewStar == 6) Rcpp::Rcout << std::fixed << "70%.. ";
  if (NewStar == 7) Rcpp::Rcout << std::fixed << "80%.. ";
  if (NewStar == 8) Rcpp::Rcout << std::fixed << "90%.. ";
  if (NewStar == 9) Rcpp::Rcout << std::fixed << "100%!";

}



// //Update burn-in progress bar----------------------------------------------------------------------------
// void UpdateBurnInBar(int s, mcmcobj McmcObj) {
//
//   //Set MCMC object
//   arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
//   int BarLength = McmcObj.BarLength;
//
//   //Number of new star in interactive mode
//   arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
//   arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
//   int NewStar = NewStarBooleanVec(0);
//   for (int i = 0; i < (BarLength + 1 - NewStar); i++) Rcpp::Rcout << std::fixed << "\b";
//   Rcpp::Rcout << std::fixed << "*";
//   for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
//   Rcpp::Rcout << std::fixed << "|";
//
// }



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateBurnInBar(int s, mcmcobj McmcObj) {

  //Set MCMC object
  arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
  int BarLength = McmcObj.BarLength;

  //Add a new star
  arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);
  Rcpp::Rcout << std::fixed << "\rBurn-in progress:  |";
  for (int i = 0; i < NewStar; i++) Rcpp::Rcout << std::fixed << "*";
  for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
  Rcpp::Rcout << std::fixed << "|";

}

