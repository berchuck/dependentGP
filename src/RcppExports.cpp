// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// H
arma::mat H(double Phi, arma::mat const& Dist, int KernInd, int T);
RcppExport SEXP _dependentGP_H(SEXP PhiSEXP, SEXP DistSEXP, SEXP KernIndSEXP, SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Dist(DistSEXP);
    Rcpp::traits::input_parameter< int >::type KernInd(KernIndSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(H(Phi, Dist, KernInd, T));
    return rcpp_result_gen;
END_RCPP
}
// lmcGP_Rcpp
Rcpp::List lmcGP_Rcpp(Rcpp::List DatObj_List, Rcpp::List HyPara_List, Rcpp::List MetrObj_List, Rcpp::List Para_List, Rcpp::List McmcObj_List, arma::mat RawSamples, bool Interactive);
RcppExport SEXP _dependentGP_lmcGP_Rcpp(SEXP DatObj_ListSEXP, SEXP HyPara_ListSEXP, SEXP MetrObj_ListSEXP, SEXP Para_ListSEXP, SEXP McmcObj_ListSEXP, SEXP RawSamplesSEXP, SEXP InteractiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type HyPara_List(HyPara_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MetrObj_List(MetrObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type McmcObj_List(McmcObj_ListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type RawSamples(RawSamplesSEXP);
    Rcpp::traits::input_parameter< bool >::type Interactive(InteractiveSEXP);
    rcpp_result_gen = Rcpp::wrap(lmcGP_Rcpp(DatObj_List, HyPara_List, MetrObj_List, Para_List, McmcObj_List, RawSamples, Interactive));
    return rcpp_result_gen;
END_RCPP
}
// PredY
Rcpp::List PredY(int NKeep, int NVisitsReturn, int K, arma::mat Gamma, int T, arma::mat Phi, arma::mat TimeDist, int KernInd, arma::uvec OriginalVisits, arma::uvec NewVisits, int NNewVisitsReturn, int NOriginalReturn, arma::uvec WhichOriginalReturnLong, arma::uvec WhichOriginalReturn, arma::uvec WhichNewReturn, arma::mat A, arma::uvec AInd, arma::mat Sigma2, int NVisits);
RcppExport SEXP _dependentGP_PredY(SEXP NKeepSEXP, SEXP NVisitsReturnSEXP, SEXP KSEXP, SEXP GammaSEXP, SEXP TSEXP, SEXP PhiSEXP, SEXP TimeDistSEXP, SEXP KernIndSEXP, SEXP OriginalVisitsSEXP, SEXP NewVisitsSEXP, SEXP NNewVisitsReturnSEXP, SEXP NOriginalReturnSEXP, SEXP WhichOriginalReturnLongSEXP, SEXP WhichOriginalReturnSEXP, SEXP WhichNewReturnSEXP, SEXP ASEXP, SEXP AIndSEXP, SEXP Sigma2SEXP, SEXP NVisitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< int >::type NVisitsReturn(NVisitsReturnSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeDist(TimeDistSEXP);
    Rcpp::traits::input_parameter< int >::type KernInd(KernIndSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type OriginalVisits(OriginalVisitsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type NewVisits(NewVisitsSEXP);
    Rcpp::traits::input_parameter< int >::type NNewVisitsReturn(NNewVisitsReturnSEXP);
    Rcpp::traits::input_parameter< int >::type NOriginalReturn(NOriginalReturnSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type WhichOriginalReturnLong(WhichOriginalReturnLongSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type WhichOriginalReturn(WhichOriginalReturnSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type WhichNewReturn(WhichNewReturnSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type AInd(AIndSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma2(Sigma2SEXP);
    Rcpp::traits::input_parameter< int >::type NVisits(NVisitsSEXP);
    rcpp_result_gen = Rcpp::wrap(PredY(NKeep, NVisitsReturn, K, Gamma, T, Phi, TimeDist, KernInd, OriginalVisits, NewVisits, NNewVisitsReturn, NOriginalReturn, WhichOriginalReturnLong, WhichOriginalReturn, WhichNewReturn, A, AInd, Sigma2, NVisits));
    return rcpp_result_gen;
END_RCPP
}
// PredDerv
Rcpp::List PredDerv(int NKeep, int NNewVisits, int K, arma::mat Phi, arma::mat Sigma2, arma::mat A, arma::uvec AInd, int N, arma::mat IndecesOrg, arma::colvec Time, arma::mat IndecesNew, arma::colvec Y, arma::colvec TimesSort);
RcppExport SEXP _dependentGP_PredDerv(SEXP NKeepSEXP, SEXP NNewVisitsSEXP, SEXP KSEXP, SEXP PhiSEXP, SEXP Sigma2SEXP, SEXP ASEXP, SEXP AIndSEXP, SEXP NSEXP, SEXP IndecesOrgSEXP, SEXP TimeSEXP, SEXP IndecesNewSEXP, SEXP YSEXP, SEXP TimesSortSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< int >::type NNewVisits(NNewVisitsSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma2(Sigma2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type AInd(AIndSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type IndecesOrg(IndecesOrgSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Time(TimeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type IndecesNew(IndecesNewSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type TimesSort(TimesSortSEXP);
    rcpp_result_gen = Rcpp::wrap(PredDerv(NKeep, NNewVisits, K, Phi, Sigma2, A, AInd, N, IndecesOrg, Time, IndecesNew, Y, TimesSort));
    return rcpp_result_gen;
END_RCPP
}
// Play
void Play();
RcppExport SEXP _dependentGP_Play() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Play();
    return R_NilValue;
END_RCPP
}
// CholInv
arma::mat CholInv(arma::mat const& Cov);
RcppExport SEXP _dependentGP_CholInv(SEXP CovSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type Cov(CovSEXP);
    rcpp_result_gen = Rcpp::wrap(CholInv(Cov));
    return rcpp_result_gen;
END_RCPP
}
// makeSymm
arma::mat makeSymm(arma::mat const& A);
RcppExport SEXP _dependentGP_makeSymm(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(makeSymm(A));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dependentGP_H", (DL_FUNC) &_dependentGP_H, 4},
    {"_dependentGP_lmcGP_Rcpp", (DL_FUNC) &_dependentGP_lmcGP_Rcpp, 7},
    {"_dependentGP_PredY", (DL_FUNC) &_dependentGP_PredY, 19},
    {"_dependentGP_PredDerv", (DL_FUNC) &_dependentGP_PredDerv, 13},
    {"_dependentGP_Play", (DL_FUNC) &_dependentGP_Play, 0},
    {"_dependentGP_CholInv", (DL_FUNC) &_dependentGP_CholInv, 1},
    {"_dependentGP_makeSymm", (DL_FUNC) &_dependentGP_makeSymm, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_dependentGP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
