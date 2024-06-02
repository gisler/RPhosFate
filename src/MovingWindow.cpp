#include <RcppArmadillo.h>
#include "MovingWindow.h"

arma::uword MovingWindow::determineFacet(const double& ns_dir_inf) {
  arma::uword us_fct{};

  if        (ns_dir_inf <   90.0 && ns_dir_inf >=  45.0) {
    us_fct = 1;
  } else if (ns_dir_inf <   45.0 && ns_dir_inf >=   0.0) {
    us_fct = 2;
  } else if (ns_dir_inf <= 360.0 && ns_dir_inf >= 315.0) {
    us_fct = 3;
  } else if (ns_dir_inf <  315.0 && ns_dir_inf >= 270.0) {
    us_fct = 4;
  } else if (ns_dir_inf <  270.0 && ns_dir_inf >= 225.0) {
    us_fct = 5;
  } else if (ns_dir_inf <  225.0 && ns_dir_inf >= 180.0) {
    us_fct = 6;
  } else if (ns_dir_inf <  180.0 && ns_dir_inf >= 135.0) {
    us_fct = 7;
  } else if (ns_dir_inf <  135.0 && ns_dir_inf >=  90.0) {
    us_fct = 8;
  } else{
    Rcpp::stop("\"dir_inf\" out of range.");
  }

  return us_fct;
};

X1X2 MovingWindow::determine_x1x2(
  const double& ns_dir_inf,
  const arma::uword& i,
  const arma::uword& j,
  const arma::dmat& nm_xxx
) {
  arma::sword is_row{static_cast<arma::sword>(i)};
  arma::sword is_col{static_cast<arma::sword>(j)};
  arma::uword us_fct{determineFacet(ns_dir_inf)};
  arma::sword is_row_x1{is_row + drdc.iv_dr_x1[us_fct]};
  arma::sword is_col_x1{is_col + drdc.iv_dc_x1[us_fct]};
  arma::sword is_row_x2{is_row + drdc.iv_dr_x2[us_fct]};
  arma::sword is_col_x2{is_col + drdc.iv_dc_x2[us_fct]};

  double x1{}, x2{};
  if (is_row_x1 == -1 || is_row_x1 == is_rws || is_col_x1 == -1 || is_col_x1 == is_cls) {
    x1 = NA_REAL;
  } else {
    x1 = nm_xxx.at(is_row_x1, is_col_x1);
  };
  if (is_row_x2 == -1 || is_row_x2 == is_rws || is_col_x2 == -1 || is_col_x2 == is_cls) {
    x2 = NA_REAL;
  } else {
    x2 = nm_xxx.at(is_row_x2, is_col_x2);
  };

  return X1X2{x1, x2};
};
