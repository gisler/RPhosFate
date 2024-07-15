#ifndef MOVINGWINDOW_H
#define MOVINGWINDOW_H

#include <RcppArmadillo.h>

const struct {
  arma::ivec8 iv_dr{-1, -1, -1,
                     0,      0,
                     1,  1,  1}; // Deltas of rows
  arma::ivec8 iv_dc{-1,  0,  1,
                    -1,      1,
                    -1,  0,  1}; // Deltas of cols

  arma::dvec8 nv_dir_min{ 90.0, 135.0, 180.0,
                          45.0,        225.0,
                           0.0, 315.0, 270.0}; // Lower bounds of inflow directions
  arma::dvec8 nv_dir_mid{135.0, 180.0, 225.0,
                          90.0,        270.0,
                          45.0, 360.0, 315.0}; // Midpoint of inflow directions
  arma::dvec8 nv_dir_max{180.0, 225.0, 270.0,
                         135.0,        315.0,
                          90.0,  45.0, 360.0}; // Upper bounds of inflow directions

  arma::uvec3 uv_oob_lr{0, 1, 2}; // Indices, which are out of bounds when row == 0
  arma::uvec3 uv_oob_ur{5, 6, 7}; // Indices, which are out of bounds when row == max row
  arma::uvec3 uv_oob_lc{0, 3, 5}; // Indices, which are out of bounds when col == 0
  arma::uvec3 uv_oob_uc{2, 4, 7}; // Indices, which are out of bounds when col == max col

  arma::uvec3 uv_oob{0, 0, 0}; // Vector for setting indices, which are out of bounds to 0
} ifl;

const struct {
  arma::ivec8 iv_x1_dr{ 0, -1, -1,  0,  0,  1, 1, 0};
  arma::ivec8 iv_x1_dc{ 1,  0,  0, -1, -1,  0, 0, 1};

  arma::ivec8 iv_x2_dr{-1, -1, -1, -1,  1,  1, 1, 1};
  arma::ivec8 iv_x2_dc{ 1,  1, -1, -1, -1, -1, 1, 1};
} fct_drdc; // Facets: deltas of rows and cols of x1 and x2

struct FacetProperties {
  bool ls_x1_oob{false}; // Is row or col of x1 out of bounds or proportion of x1 == 0.0?
  double ns_p1{};        // Proportion of x1
  arma::uword us_x1_r{}; // Row of x1
  arma::uword us_x1_c{}; // Col of x1

  bool ls_x2_oob{false}; // Is row or col of x2 out of bounds or proportion of x2 == 0.0?
  double ns_p2{};        // Proportion of x2
  arma::uword us_x2_r{}; // Row of x2
  arma::uword us_x2_c{}; // Col of x2
};

template <typename T>
struct X1X2 {
  T x1{};
  T x2{};
};

struct CalcOrder {
  std::vector<arma::uword> uv_r{};
  std::vector<arma::uword> uv_c{};

  CalcOrder(arma::uword n):
    uv_r{},
    uv_c{}
  {
    uv_r.reserve(n);
    uv_c.reserve(n);
  }
};

class MovingWindow {
public:
  const arma::uword us_rws{};
  const arma::uword us_cls{};
  const arma::sword is_rws{};
  const arma::sword is_cls{};

  MovingWindow(arma::uword us_rws, arma::uword us_cls):
    us_rws{us_rws},
    us_cls{us_cls},
    is_rws{static_cast<arma::sword>(us_rws)},
    is_cls{static_cast<arma::sword>(us_cls)}
  {}

  FacetProperties determineFacetProperties(
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

  arma::dvec8 get_ifl_p(
    const arma::dmat& nm_dir_inf,
    const arma::uword& us_row,
    const arma::uword& us_col
  );

  template <typename T>
  arma::Col<T> get_ifl(
    const arma::dvec8& nv_ifl_p,
    const arma::uword& us_row,
    const arma::uword& us_col,
    const arma::Mat<T>& xm_xxx,
    const T NA
  );
};

#endif

inline FacetProperties MovingWindow::determineFacetProperties(
  const double& ns_dir_inf,
  const arma::uword& us_row,
  const arma::uword& us_col
) {
  arma::uword us_fct{};
  FacetProperties fct{};

  // Determine facet and proportions of x1 and x2
  if        (ns_dir_inf >=  45.0 && ns_dir_inf <   90.0) {
    us_fct = 0;

    fct.ns_p1 = (ns_dir_inf - 45.0) / 45.0;
    fct.ns_p2 = 1.0 - fct.ns_p1;

  } else if (ns_dir_inf >=   0.0 && ns_dir_inf <   45.0) {
    us_fct = 1;

    fct.ns_p2 = ns_dir_inf / 45.0;
    fct.ns_p1 = 1.0 - fct.ns_p2;

  } else if (ns_dir_inf >= 315.0 && ns_dir_inf <= 360.0) {
    us_fct = 2;

    fct.ns_p1 = (ns_dir_inf - 315.0) / 45.0;
    fct.ns_p2 = 1.0 - fct.ns_p1;

  } else if (ns_dir_inf >= 270.0 && ns_dir_inf <  315.0) {
    us_fct = 3;

    fct.ns_p2 = (ns_dir_inf - 270.0) / 45.0;
    fct.ns_p1 = 1.0 - fct.ns_p2;

  } else if (ns_dir_inf >= 225.0 && ns_dir_inf <  270.0) {
    us_fct = 4;

    fct.ns_p1 = (ns_dir_inf - 225.0) / 45.0;
    fct.ns_p2 = 1.0 - fct.ns_p1;

  } else if (ns_dir_inf >= 180.0 && ns_dir_inf <  225.0) {
    us_fct = 5;

    fct.ns_p2 = (ns_dir_inf - 180.0) / 45.0;
    fct.ns_p1 = 1.0 - fct.ns_p2;

  } else if (ns_dir_inf >= 135.0 && ns_dir_inf <  180.0) {
    us_fct = 6;

    fct.ns_p1 = (ns_dir_inf - 135.0) / 45.0;
    fct.ns_p2 = 1.0 - fct.ns_p1;

  } else if (ns_dir_inf >=  90.0 && ns_dir_inf <  135.0) {
    us_fct = 7;

    fct.ns_p2 = (ns_dir_inf - 90.0) / 45.0;
    fct.ns_p1 = 1.0 - fct.ns_p2;

  } else {
    // The first two crash R (v4.4.0 with Rcpp v1.0.12), but would be better
    // choices:
    // Rcpp::stop("\"dir_inf\" out of range.");
    // throw Rcpp::exception("\"dir_inf\" out of range.");
    Rcpp::Rcerr << "Warning: \"dir_inf\" out of range." << std::endl;

    fct.ls_x1_oob = true;
    fct.ls_x2_oob = true;
  }

  // Determine rows and cols of x1 and x2
  arma::sword is_row{static_cast<arma::sword>(us_row)};
  arma::sword is_col{static_cast<arma::sword>(us_col)};
  arma::sword is_x1_row{is_row + fct_drdc.iv_x1_dr[us_fct]};
  arma::sword is_x1_col{is_col + fct_drdc.iv_x1_dc[us_fct]};
  arma::sword is_x2_row{is_row + fct_drdc.iv_x2_dr[us_fct]};
  arma::sword is_x2_col{is_col + fct_drdc.iv_x2_dc[us_fct]};

  if (is_x1_row == -1 || is_x1_row == is_rws ||
      is_x1_col == -1 || is_x1_col == is_cls ||
      fct.ns_p1 == 0.0) {
    fct.ls_x1_oob = true;
  } else {
    fct.us_x1_r = us_row + fct_drdc.iv_x1_dr[us_fct];
    fct.us_x1_c = us_col + fct_drdc.iv_x1_dc[us_fct];
  }
  if (is_x2_row == -1 || is_x2_row == is_rws ||
      is_x2_col == -1 || is_x2_col == is_cls ||
      fct.ns_p2 == 0.0) {
    fct.ls_x2_oob = true;
  } else {
    fct.us_x2_r = us_row + fct_drdc.iv_x2_dr[us_fct];
    fct.us_x2_c = us_col + fct_drdc.iv_x2_dc[us_fct];
  }

  return fct;
}

template <typename T>
inline X1X2<T> MovingWindow::get_x1x2(
  const double& ns_dir_inf,
  const arma::uword& i,
  const arma::uword& j,
  const arma::Mat<T>& xm_xxx,
  const T NA
) {
  FacetProperties fct{determineFacetProperties(ns_dir_inf, i, j)};
  T x1{}, x2{};

  if (fct.ls_x1_oob) {
    x1 = NA;
  } else {
    x1 = xm_xxx.at(fct.us_x1_r, fct.us_x1_c);
  }
  if (fct.ls_x2_oob) {
    x2 = NA;
  } else {
    x2 = xm_xxx.at(fct.us_x2_r, fct.us_x2_c);
  }

  return X1X2<T>{x1, x2};
}

inline arma::dvec8 MovingWindow::get_ifl_p(
  const arma::dmat& nm_dir_inf,
  const arma::uword& us_row,
  const arma::uword& us_col
) {
  arma::uvec8 uv_cll(arma::fill::value(1));

  if (us_row == 0) {
    uv_cll.elem(ifl.uv_oob_lr) = ifl.uv_oob;
  }
  if (us_row == us_rws - 1) {
    uv_cll.elem(ifl.uv_oob_ur) = ifl.uv_oob;
  }
  if (us_col == 0) {
    uv_cll.elem(ifl.uv_oob_lc) = ifl.uv_oob;
  }
  if (us_col == us_cls - 1) {
    uv_cll.elem(ifl.uv_oob_uc) = ifl.uv_oob;
  }

  double ns_dir_inf{};
  arma::dvec8 nv_ifl_p(arma::fill::value(NA_REAL));

  for (arma::uword k = 0; k < uv_cll.n_elem; ++k) {
    if (uv_cll[k] == 1) {
      ns_dir_inf = nm_dir_inf.at(us_row + ifl.iv_dr[k], us_col + ifl.iv_dc[k]);

      if (k == 6) {
        if (ns_dir_inf > ifl.nv_dir_min[k] || ns_dir_inf < ifl.nv_dir_max[k]) {
          if (ns_dir_inf > ifl.nv_dir_min[k]) {
            nv_ifl_p[k] = (ns_dir_inf - ifl.nv_dir_min[k]) / 45.0;
          } else {
            nv_ifl_p[k] = (ifl.nv_dir_max[k] - ns_dir_inf) / 45.0;
          }
        }
      } else {
        if (ns_dir_inf > ifl.nv_dir_min[k] && ns_dir_inf < ifl.nv_dir_max[k]) {
          if (ns_dir_inf <= ifl.nv_dir_mid[k]) {
            nv_ifl_p[k] = (ns_dir_inf - ifl.nv_dir_min[k]) / 45.0;
          } else {
            nv_ifl_p[k] = (ifl.nv_dir_max[k] - ns_dir_inf) / 45.0;
          }
        }
      }
    }
  }

  return nv_ifl_p;
}

template <typename T>
inline arma::Col<T> MovingWindow::get_ifl(
  const arma::dvec8& nv_ifl_p,
  const arma::uword& us_row,
  const arma::uword& us_col,
  const arma::Mat<T>& xm_xxx,
  const T NA
) {
  arma::Col<T> xv_ifl(8, arma::fill::value(NA));

  for (arma::uword k = 0; k < nv_ifl_p.n_elem; ++k) {
    if (nv_ifl_p[k] > 0.0) {
      xv_ifl[k] = xm_xxx.at(us_row + ifl.iv_dr[k], us_col + ifl.iv_dc[k]);
    }
  }

  return xv_ifl;
}
