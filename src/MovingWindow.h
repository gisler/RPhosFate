#ifndef MOVINGWINDOW_H
#define MOVINGWINDOW_H

#include <RcppArmadillo.h>

const struct {
  arma::ivec8 iv_x1_dr{ 0, -1, -1,  0,  0,  1, 1, 0};
  arma::ivec8 iv_x1_dc{ 1,  0,  0, -1, -1,  0, 0, 1};

  arma::ivec8 iv_x2_dr{-1, -1, -1, -1,  1,  1, 1, 1};
  arma::ivec8 iv_x2_dc{ 1,  1, -1, -1, -1, -1, 1, 1};
} drdc; // Facets: deltas of rows and cols of x1 and x2

struct FacetRowsCols {
  bool ls_x1_oob{false}; // Is row or col of x1 out of bounds?
  arma::uword us_x1_r{}; // Row of x1
  arma::uword us_x1_c{}; // Col of x1

  bool ls_x2_oob{false}; // Is row or col of x2 out of bounds?
  arma::uword us_x2_r{}; // Row of x2
  arma::uword us_x2_c{}; // Col of x2
};

template <typename T>
struct X1X2 {
  T x1{};
  T x2{};
};

class MovingWindow {
public:
  const arma::sword is_rws{};
  const arma::sword is_cls{};

  FacetRowsCols determineFacetRowsCols(
    const double& ns_dir_inf,
    const arma::uword& us_row,
    const arma::uword& us_col
  );

  template <typename T>
  X1X2<T> get_x1x2(
    const double& ns_dir_inf,
    const arma::uword& i,
    const arma::uword& j,
    const arma::Mat<T>& xm_xxx,
    const T NA
  );
};

#endif

FacetRowsCols MovingWindow::determineFacetRowsCols(
  const double& ns_dir_inf,
  const arma::uword& us_row,
  const arma::uword& us_col
) {
  // Determine facet
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
  } else {
    // The first two crash R (v4.4.0 with Rcpp v1.0.12), but would be better
    // choices:
    // Rcpp::stop("\"dir_inf\" out of range.");
    // throw Rcpp::exception("\"dir_inf\" out of range.");
    Rcpp::Rcerr << "Warning: \"dir_inf\" out of range." << std::endl;
  }

  // Determine rows and cols of x1 and x2
  arma::sword is_row{static_cast<arma::sword>(us_row)};
  arma::sword is_col{static_cast<arma::sword>(us_col)};
  arma::sword is_x1_row{is_row + drdc.iv_x1_dr[us_fct]};
  arma::sword is_x1_col{is_col + drdc.iv_x1_dc[us_fct]};
  arma::sword is_x2_row{is_row + drdc.iv_x2_dr[us_fct]};
  arma::sword is_x2_col{is_col + drdc.iv_x2_dc[us_fct]};

  FacetRowsCols fctRC{};
  if (is_x1_row == -1 || is_x1_row == is_rws ||
      is_x1_col == -1 || is_x1_col == is_cls) {
    fctRC.ls_x1_oob = true;
  } else {
    fctRC.us_x1_r = us_row + drdc.iv_x1_dr[us_fct];
    fctRC.us_x1_c = us_col + drdc.iv_x1_dc[us_fct];
  }
  if (is_x2_row == -1 || is_x2_row == is_rws ||
      is_x2_col == -1 || is_x2_col == is_cls) {
    fctRC.ls_x2_oob = true;
  } else {
    fctRC.us_x2_r = us_row + drdc.iv_x2_dr[us_fct];
    fctRC.us_x2_c = us_col + drdc.iv_x2_dc[us_fct];
  }

  return fctRC;
}

template <typename T>
X1X2<T> MovingWindow::get_x1x2(
    const double& ns_dir_inf,
    const arma::uword& i,
    const arma::uword& j,
    const arma::Mat<T>& xm_xxx,
    const T NA
) {
  FacetRowsCols fctRC{determineFacetRowsCols(ns_dir_inf, i, j)};

  T x1{}, x2{};
  if (fctRC.ls_x1_oob) {
    x1 = NA;
  } else {
    x1 = xm_xxx.at(fctRC.us_x1_r, fctRC.us_x1_c);
  }
  if (fctRC.ls_x2_oob) {
    x2 = NA;
  } else {
    x2 = xm_xxx.at(fctRC.us_x2_r, fctRC.us_x2_c);
  }

  return X1X2<T>{x1, x2};
}
