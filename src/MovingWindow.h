#ifndef MOVINGWINDOW_H
#define MOVINGWINDOW_H

#include <RcppArmadillo.h>

// Facet elevation modifiers
const struct {
  arma::ivec8 iv_dr_x1{ 0, -1, -1,  0,  0,  1, 1, 0};
  arma::ivec8 iv_dc_x1{ 1,  0,  0, -1, -1,  0, 0, 1};

  arma::ivec8 iv_dr_x2{-1, -1, -1, -1,  1,  1, 1, 1};
  arma::ivec8 iv_dc_x2{ 1,  1, -1, -1, -1, -1, 1, 1};
} drdc;

template <typename T>
struct X1X2 {
  T x1{};
  T x2{};
};

class MovingWindow {
public:
  const arma::sword is_rws{};
  const arma::sword is_cls{};

  arma::uword determineFacet(const double& ns_dir_inf);

  template <typename T>
  X1X2<T> determine_x1x2(
    const double& ns_dir_inf,
    const arma::uword& i,
    const arma::uword& j,
    const arma::Mat<T>& xm_xxx,
    const T NA
  );
};

#endif

arma::uword MovingWindow::determineFacet(const double& ns_dir_inf) {
  arma::uword us_fct{};

  if        (ns_dir_inf <   90.0 && ns_dir_inf >=  45.0) {
    us_fct = 0;
  } else if (ns_dir_inf <   45.0 && ns_dir_inf >=   0.0) {
    us_fct = 1;
  } else if (ns_dir_inf <= 360.0 && ns_dir_inf >= 315.0) {
    us_fct = 2;
  } else if (ns_dir_inf <  315.0 && ns_dir_inf >= 270.0) {
    us_fct = 3;
  } else if (ns_dir_inf <  270.0 && ns_dir_inf >= 225.0) {
    us_fct = 4;
  } else if (ns_dir_inf <  225.0 && ns_dir_inf >= 180.0) {
    us_fct = 5;
  } else if (ns_dir_inf <  180.0 && ns_dir_inf >= 135.0) {
    us_fct = 6;
  } else if (ns_dir_inf <  135.0 && ns_dir_inf >=  90.0) {
    us_fct = 7;
  } else{
    // The first two crash R (v4.4.0 with Rcpp v1.0.12), but would be better
    // choices:
    // Rcpp::stop("\"dir_inf\" out of range.");
    // throw Rcpp::exception("\"dir_inf\" out of range.");
    Rcpp::Rcerr << "Warning: \"dir_inf\" out of range." << std::endl;
  }

  return us_fct;
};

template <typename T>
X1X2<T> MovingWindow::determine_x1x2(
    const double& ns_dir_inf,
    const arma::uword& i,
    const arma::uword& j,
    const arma::Mat<T>& xm_xxx,
    const T NA
) {
  arma::uword us_fct{determineFacet(ns_dir_inf)};
  arma::sword is_row{static_cast<arma::sword>(i)};
  arma::sword is_col{static_cast<arma::sword>(j)};
  arma::sword is_row_x1{is_row + drdc.iv_dr_x1[us_fct]};
  arma::sword is_col_x1{is_col + drdc.iv_dc_x1[us_fct]};
  arma::sword is_row_x2{is_row + drdc.iv_dr_x2[us_fct]};
  arma::sword is_col_x2{is_col + drdc.iv_dc_x2[us_fct]};

  T x1{}, x2{};
  if (is_row_x1 == -1 || is_row_x1 == is_rws ||
      is_col_x1 == -1 || is_col_x1 == is_cls) {
    x1 = NA;
  } else {
    x1 = xm_xxx.at(is_row_x1, is_col_x1);
  };
  if (is_row_x2 == -1 || is_row_x2 == is_rws ||
      is_col_x2 == -1 || is_col_x2 == is_cls) {
    x2 = NA;
  } else {
    x2 = xm_xxx.at(is_row_x2, is_col_x2);
  };

  return X1X2<T>{x1, x2};
};
