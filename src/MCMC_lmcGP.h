//
//  functions used in dependentGP package
//

#ifndef __dependentGP__
#define __dependentGP__

//MCMC Sampler
Rcpp::List lmcGP_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                      Rcpp::List MetrObj_List, Rcpp::List Para_List,  
                      Rcpp::List McmcObj_List, arma::mat RawSamples, 
                      bool Interactive, bool Verbose);
  
//STRUCT DEFINITIONS
struct datobj {
  arma::colvec Y;
  arma::mat YWide;
  int T;
  int N;
  int K;
  arma::mat EyeT;
  arma::mat EyeK;
  arma::mat EyeTK;
  arma::colvec ZeroT;
  arma::mat TimeDist;
  arma::colvec Time;
  arma::umat AInd;
  int KernInd;
};
struct hypara {
  double Alpha;
  double Beta;
  arma::colvec APhi;
  arma::colvec BPhi;
  double Xi;
  arma::mat Psi;
};
struct metrobj {
  arma::colvec MetropT;
  arma::colvec MetropPhi;
  arma::Col<int> AcceptanceT;
  arma::Col<int> AcceptancePhi;
  arma::colvec OriginalTuners;
};
struct para {
  arma::colvec Gamma;
  arma::mat GammaWide;
  arma::colvec Theta;
  arma::colvec Sigma2;
  arma::mat Sigma;
  arma::mat SigmaInv;
  arma::mat CholSigma;
  arma::colvec Phi;
  arma::cube HPhi;
  arma::cube HPhiInv;
  arma::cube CholHPhi;
  arma::mat TMat;
  double logDetT;
  arma::mat TInv;
  arma::mat A;
};
struct mcmcobj {
  int NBurn;
  int NSims;
  int NThin;
  int NPilot;
  int NTotal;
  int NKeep;
  arma::vec WhichKeep;
  arma::vec WhichPilotAdapt;
  arma::vec WhichBurnInProgress;
  arma::vec WhichBurnInProgressInt;
  arma::vec WhichSamplerProgress;
  arma::vec BurnInProgress;
  int BarLength;
  int PilotAdaptDenominator;
};

//COVARIANCE FUNCTIONS
arma::mat H(double Phi, arma::mat const& Dist, int KernInd, int T);
  
//DISTRIBUTION FUNCTIONS
arma::vec rnormRcpp(int n, double mean, double sd);
// arma::vec sampleRcpp(arma::Col<int> const& x, int size, bool replace, arma::vec const& prob);
// double rtnormRcppMSM(double mean, double sd, double lower, double upper);
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
// double pnormRcpp(double q);
// double lpnormRcpp(double q);
// double UpperpnormRcpp(double q);
// double lUpperpnormRcpp(double q);
double rigammaRcpp(double Alpha, double Theta);
// double rgammaRcpp(double Alpha, double Theta);
// arma::mat rwishRcpp(double n, arma::mat const& V);
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti);
double randuRcpp();
// double rtnormRcpp(double mean, double sd, bool Above);
// arma::vec rtnormRcppMSM(int N, arma::vec const& mean, arma::vec const& sd, double lower, double upper);
arma::vec rnormSNRcpp(int n);
arma::vec rnormVecRcpp(arma::vec const& mean, arma::vec const& sd);
arma::colvec rmvnormRcppRobust(arma::colvec const& Mu, arma::mat const& Sigma);
  
//MCMC CONVERSION FUNCTIONS
datobj ConvertDatObj(Rcpp::List DatObj_List);
hypara ConvertHyPara(Rcpp::List HyPara_List);
metrobj ConvertMetrObj(Rcpp::List MetrObj_List);
para ConvertPara(Rcpp::List Para_List);
mcmcobj ConvertMcmcObj(Rcpp::List McmcObj_List);

//MCMC SAMPLER FUNCTIONS
para SampleSigma2(datobj DatObj, para Para, hypara HyPara);
para SampleGamma(datobj DatObj, para Para);
std::pair<para, metrobj> SamplePhi(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj);
std::pair<para, metrobj> SampleA(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj);

//MCMC UTILITY FUNCTIONS
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive);
Rcpp::List OutputMetrObj(metrobj MetrObj, datobj DatObj);
metrobj PilotAdaptation(metrobj MetrObj, mcmcobj McmcObj, datobj DatObj);
void SamplerProgress(int s, mcmcobj McmcObj);
arma::colvec StoreSamples(datobj DatObj, para Para);
void UpdateBurnInBar(int s, mcmcobj McmcObj);
void UpdateBurnInBarInt(int s, mcmcobj McmcObj);

//UTILITY FUNCTIONS
arma::mat CholInv(arma::mat const& Cov);
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol);
arma::mat makeSymm(arma::mat const& A);

#endif // __dependentGP__
