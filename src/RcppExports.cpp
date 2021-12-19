// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// D8slope
arma::dmat D8slope(arma::imat& im_dir, arma::dmat& nm_dem, arma::imat& im_fDo, double ns_fpl, int is_ths);
RcppExport SEXP _RPhosFate_D8slope(SEXP im_dirSEXP, SEXP nm_demSEXP, SEXP im_fDoSEXP, SEXP ns_fplSEXP, SEXP is_thsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type im_dir(im_dirSEXP);
    Rcpp::traits::input_parameter< arma::dmat& >::type nm_dem(nm_demSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type im_fDo(im_fDoSEXP);
    Rcpp::traits::input_parameter< double >::type ns_fpl(ns_fplSEXP);
    Rcpp::traits::input_parameter< int >::type is_ths(is_thsSEXP);
    rcpp_result_gen = Rcpp::wrap(D8slope(im_dir, nm_dem, im_fDo, ns_fpl, is_ths));
    return rcpp_result_gen;
END_RCPP
}
// dir_sth
arma::imat dir_sth(arma::imat& im_dir, arma::imat& im_sth, arma::imat& im_fDo, int is_ths);
RcppExport SEXP _RPhosFate_dir_sth(SEXP im_dirSEXP, SEXP im_sthSEXP, SEXP im_fDoSEXP, SEXP is_thsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type im_dir(im_dirSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type im_sth(im_sthSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type im_fDo(im_fDoSEXP);
    Rcpp::traits::input_parameter< int >::type is_ths(is_thsSEXP);
    rcpp_result_gen = Rcpp::wrap(dir_sth(im_dir, im_sth, im_fDo, is_ths));
    return rcpp_result_gen;
END_RCPP
}
// transportCpp
List transportCpp(S4 parameters, double ns_dep_ovl, double ns_tfc_inl, S4 helpers, S4 order, IntegerMatrix im_cha, IntegerMatrix im_dir, IntegerMatrix im_inl, IntegerMatrix im_rip, NumericMatrix nm_man, NumericMatrix nm_xxe, NumericMatrix nm_rhy, NumericMatrix nm_slp);
RcppExport SEXP _RPhosFate_transportCpp(SEXP parametersSEXP, SEXP ns_dep_ovlSEXP, SEXP ns_tfc_inlSEXP, SEXP helpersSEXP, SEXP orderSEXP, SEXP im_chaSEXP, SEXP im_dirSEXP, SEXP im_inlSEXP, SEXP im_ripSEXP, SEXP nm_manSEXP, SEXP nm_xxeSEXP, SEXP nm_rhySEXP, SEXP nm_slpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< double >::type ns_dep_ovl(ns_dep_ovlSEXP);
    Rcpp::traits::input_parameter< double >::type ns_tfc_inl(ns_tfc_inlSEXP);
    Rcpp::traits::input_parameter< S4 >::type helpers(helpersSEXP);
    Rcpp::traits::input_parameter< S4 >::type order(orderSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type im_cha(im_chaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type im_dir(im_dirSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type im_inl(im_inlSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type im_rip(im_ripSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nm_man(nm_manSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nm_xxe(nm_xxeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nm_rhy(nm_rhySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nm_slp(nm_slpSEXP);
    rcpp_result_gen = Rcpp::wrap(transportCpp(parameters, ns_dep_ovl, ns_tfc_inl, helpers, order, im_cha, im_dir, im_inl, im_rip, nm_man, nm_xxe, nm_rhy, nm_slp));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RPhosFate_D8slope", (DL_FUNC) &_RPhosFate_D8slope, 5},
    {"_RPhosFate_dir_sth", (DL_FUNC) &_RPhosFate_dir_sth, 4},
    {"_RPhosFate_transportCpp", (DL_FUNC) &_RPhosFate_transportCpp, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_RPhosFate(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
