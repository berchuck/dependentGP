#include <RcppArmadillo.h>
#include "MCMC_lmcGP.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobj ConvertDatObj(Rcpp::List DatObj_List) {

  //Set objects from List
  arma::colvec Y = DatObj_List["Y"];
  arma::mat YWide = DatObj_List["YWide"];
  int T = DatObj_List["T"];
  int N = DatObj_List["N"];
  int K = DatObj_List["K"];
  arma::mat EyeT = DatObj_List["EyeT"];
  arma::mat EyeK = DatObj_List["EyeK"];
  arma::mat EyeTK = DatObj_List["EyeTK"];
  arma::colvec ZeroT = DatObj_List["ZeroT"];
  arma::mat TimeDist = DatObj_List["TimeDist"];
  arma::colvec Time = DatObj_List["Time"];
  arma::umat AInd = DatObj_List["AInd"];
  int KernInd = DatObj_List["KernInd"];
  
  //Convert to C++ struct
  datobj DatObj;
  DatObj.Y = Y;
  DatObj.YWide = YWide;
  DatObj.T = T;
  DatObj.N = N;
  DatObj.K = K;
  DatObj.EyeT = EyeT;
  DatObj.EyeK = EyeK;
  DatObj.EyeTK = EyeTK;
  DatObj.ZeroT = ZeroT;
  DatObj.TimeDist = TimeDist;
  DatObj.Time = Time;
  DatObj.AInd = AInd;
  DatObj.KernInd = KernInd;
  return DatObj;

}



//Function to convert Rcpp::List HyPara to a custom C++ struct hypara--------------------------------------------------
hypara ConvertHyPara(Rcpp::List HyPara_List) {

  //Set objects from List
  double Alpha = HyPara_List["Alpha"];
  double Beta = HyPara_List["Beta"];
  arma::colvec APhi = HyPara_List["APhi"];
  arma::colvec BPhi = HyPara_List["BPhi"];
  double Xi = HyPara_List["Xi"];
  arma::mat Psi = HyPara_List["Psi"];
  
  //Convert to C++ struct
  hypara HyPara;
  HyPara.Alpha = Alpha;
  HyPara.Beta = Beta;
  HyPara.APhi = APhi;
  HyPara.BPhi = BPhi;
  HyPara.Xi = Xi;
  HyPara.Psi = Psi;
  return HyPara;

}



//Function to convert Rcpp::List MetrObj to a custom C++ struct metrobj-----------------------------------------------
metrobj ConvertMetrObj(Rcpp::List MetrObj_List) {

  //Set objects from List
  arma::colvec MetropT = MetrObj_List["MetropT"];
  arma::colvec MetropPhi = MetrObj_List["MetropPhi"];
  arma::Col<int> AcceptanceT = MetrObj_List["AcceptanceT"];
  arma::Col<int> AcceptancePhi = MetrObj_List["AcceptancePhi"];
  arma::colvec OriginalTuners = MetrObj_List["OriginalTuners"];
  
  //Convert to C++ struct
  metrobj MetrObj;
  MetrObj.MetropT = MetropT;
  MetrObj.MetropPhi = MetropPhi;
  MetrObj.AcceptanceT = AcceptanceT;
  MetrObj.AcceptancePhi = AcceptancePhi;
  MetrObj.OriginalTuners = OriginalTuners;
  return MetrObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
para ConvertPara(Rcpp::List Para_List) {
 
  //Set objects from List
  arma::colvec Gamma = Para_List["Gamma"];
  arma::mat GammaWide = Para_List["GammaWide"];
  arma::colvec Theta = Para_List["Theta"];
  arma::colvec Sigma2 = Para_List["Sigma2"];
  arma::mat Sigma = Para_List["Sigma"];
  arma::mat SigmaInv = Para_List["SigmaInv"];
  arma::mat CholSigma = Para_List["CholSigma"];
  arma::colvec Phi = Para_List["Phi"];
  arma::cube HPhi = Para_List["HPhi"];
  arma::cube HPhiInv = Para_List["HPhiInv"];
  arma::cube CholHPhi = Para_List["CholHPhi"];
  arma::mat TMat = Para_List["TMat"];
  double logDetT = Para_List["logDetT"];
  arma::mat TInv = Para_List["TInv"];
  arma::mat A = Para_List["A"];
  
  //Convert to C++ struct
  para Para;
  Para.Gamma = Gamma;
  Para.GammaWide = GammaWide;
  Para.Theta = Theta;
  Para.Sigma2 = Sigma2;
  Para.Sigma = Sigma;
  Para.SigmaInv = SigmaInv;
  Para.CholSigma = CholSigma;
  Para.Phi = Phi;
  Para.HPhi = HPhi;
  Para.HPhiInv = HPhiInv;
  Para.CholHPhi = CholHPhi;
  Para.TMat = TMat;
  Para.logDetT = logDetT;
  Para.TInv = TInv;
  Para.A = A;
  return Para;
}



//Function to convert Rcpp::List McmcObj to a custom C++ struct mcmcmobj-----------------------------------------------------
mcmcobj ConvertMcmcObj(Rcpp::List McmcObj_List) {

  //Set objects from List
  int NBurn = McmcObj_List["NBurn"];
  int NSims = McmcObj_List["NSims"];
  int NThin = McmcObj_List["NThin"];
  int NPilot = McmcObj_List["NPilot"];
  int NTotal = McmcObj_List["NTotal"];
  int NKeep = McmcObj_List["NKeep"];
  arma::vec WhichKeep = McmcObj_List["WhichKeep"];
  arma::vec WhichPilotAdapt = McmcObj_List["WhichPilotAdapt"];
  arma::vec WhichBurnInProgress = McmcObj_List["WhichBurnInProgress"];
  arma::vec WhichBurnInProgressInt = McmcObj_List["WhichBurnInProgressInt"];
  arma::vec WhichSamplerProgress = McmcObj_List["WhichSamplerProgress"];
  arma::vec BurnInProgress = McmcObj_List["BurnInProgress"];
  int BarLength = McmcObj_List["BarLength"];
  int PilotAdaptDenominator = McmcObj_List["PilotAdaptDenominator"];

  //Convert to C++ struct
  mcmcobj McmcObj;
  McmcObj.NBurn = NBurn;
  McmcObj.NSims = NSims;
  McmcObj.NThin = NThin;
  McmcObj.NPilot = NPilot;
  McmcObj.NTotal = NTotal;
  McmcObj.NKeep = NKeep;
  McmcObj.WhichKeep = WhichKeep;
  McmcObj.WhichPilotAdapt = WhichPilotAdapt;
  McmcObj.WhichBurnInProgress = WhichBurnInProgress;
  McmcObj.WhichBurnInProgressInt = WhichBurnInProgressInt;
  McmcObj.WhichSamplerProgress = WhichSamplerProgress;
  McmcObj.BurnInProgress = BurnInProgress;
  McmcObj.PilotAdaptDenominator = PilotAdaptDenominator;
  McmcObj.BarLength = BarLength;
  return McmcObj;
}



