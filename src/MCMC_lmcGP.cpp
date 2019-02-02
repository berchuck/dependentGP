#include <RcppArmadillo.h>
#include "MCMC_lmcGP.h"

//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List lmcGP_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                      Rcpp::List MetrObj_List, Rcpp::List Para_List,  
                      Rcpp::List McmcObj_List, arma::mat RawSamples, 
                      bool Interactive, bool Verbose) {
  
  //Convet Rcpp::Lists to C++ structs
  datobj DatObj = ConvertDatObj(DatObj_List);
  hypara HyPara = ConvertHyPara(HyPara_List);
  metrobj MetrObj = ConvertMetrObj(MetrObj_List);
  para Para = ConvertPara(Para_List);
  mcmcobj McmcObj = ConvertMcmcObj(McmcObj_List);

  //Set objects to be used in MCMC sampler
  int NTotal = McmcObj.NTotal;
  int NBurn = McmcObj.NBurn;
  arma::vec WhichPilotAdapt = McmcObj.WhichPilotAdapt;
  arma::vec WhichKeep = McmcObj.WhichKeep;
  arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
  arma::vec WhichBurnInProgressInt = McmcObj.WhichBurnInProgressInt;
  arma::vec WhichSamplerProgress = McmcObj.WhichSamplerProgress;
  std::pair<para, metrobj> Update;

  //User output
  if (Verbose) BeginBurnInProgress(McmcObj, Interactive);

  //Begin MCMC Sampler
  for (int s = 1; s < NTotal + 1; s++) {

    //Check for user interrupt every 500 iterations
    if (s % 500 == 0) Rcpp::checkUserInterrupt();

    //Gibbs step for Sigma2
    Para = SampleSigma2(DatObj, Para, HyPara);

    //Gibbs step for Gamma
    Para = SampleGamma(DatObj, Para);

    // Metropolis step for Phi
    Update = SamplePhi(DatObj, Para, HyPara, MetrObj);
    Para = Update.first;
    MetrObj = Update.second;

    // Metropolis step for A
    Update = SampleA(DatObj, Para, HyPara, MetrObj);
    Para = Update.first;
    MetrObj = Update.second;
    
    //Pilot adaptation
    if (std::find(WhichPilotAdapt.begin(), WhichPilotAdapt.end(), s) != WhichPilotAdapt.end())
      MetrObj = PilotAdaptation(MetrObj, McmcObj, DatObj);

    //Store raw samples
    if (std::find(WhichKeep.begin(), WhichKeep.end(), s) != WhichKeep.end())
      RawSamples.cols(find(s == WhichKeep)) = StoreSamples(DatObj, Para);

    //Update burn-in progress bar
    if (Verbose) if (Interactive) if (std::find(WhichBurnInProgress.begin(), WhichBurnInProgress.end(), s) != WhichBurnInProgress.end())
      UpdateBurnInBar(s, McmcObj);
    if (Verbose) if (!Interactive) if (std::find(WhichBurnInProgressInt.begin(), WhichBurnInProgressInt.end(), s) != WhichBurnInProgressInt.end())
      UpdateBurnInBarInt(s, McmcObj);

    //Post burn-in progress
    if (Verbose) if (s == NBurn) Rcpp::Rcout << std::fixed << "\nSampler progress:  0%.. ";
    if (Verbose) if (std::find(WhichSamplerProgress.begin(), WhichSamplerProgress.end(), s) != WhichSamplerProgress.end())
       SamplerProgress(s, McmcObj);

  //End MCMC Sampler
  }

  //Output Metropolis object for summary
  Rcpp::List Metropolis = OutputMetrObj(MetrObj, DatObj);

  //Return raw samples
  return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                            Rcpp::Named("metropolis") = Metropolis);

  return DatObj_List;
  
//End MCMC sampler function
}
