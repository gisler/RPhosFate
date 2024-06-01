#include <RcppArmadillo.h>
#include "MovingWindow.h"

arma::uword MovingWindow::determineFacet(const double ns_dir_inf) {
  if        (ns_dir_inf <   90.0 && ns_dir_inf >=  45.0) {
    return 1;
  } else if (ns_dir_inf <   45.0 && ns_dir_inf >=   0.0) {
    return 2;
  } else if (ns_dir_inf <= 360.0 && ns_dir_inf >= 315.0) {
    return 3;
  } else if (ns_dir_inf <  315.0 && ns_dir_inf >= 270.0) {
    return 4;
  } else if (ns_dir_inf <  270.0 && ns_dir_inf >= 225.0) {
    return 5;
  } else if (ns_dir_inf <  225.0 && ns_dir_inf >= 180.0) {
    return 6;
  } else if (ns_dir_inf <  180.0 && ns_dir_inf >= 135.0) {
    return 7;
  } else if (ns_dir_inf <  135.0 && ns_dir_inf >=  90.0) {
    return 8;
  }
};

X1X2 MovingWindow::determine_x1x2(
  const double ns_dir_inf,
  const arma::uword i,
  const arma::uword j,
  const arma::dmat& nm_xxx
) {
  arma::uword is_fct{determineFacet(ns_dir_inf)};

  arma::sword is_row{static_cast<arma::sword>(i)};
  arma::sword is_col{static_cast<arma::sword>(j)};
  arma::sword is_row_x1{is_row + drdc.iv_dr_x1[is_fct]};
  arma::sword is_col_x1{is_col + drdc.iv_dc_x1[is_fct]};
  arma::sword is_row_x2{is_row + drdc.iv_dr_x2[is_fct]};
  arma::sword is_col_x2{is_col + drdc.iv_dc_x2[is_fct]};

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
