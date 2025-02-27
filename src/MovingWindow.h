#ifndef FOCALWINDOW_H
#define FOCALWINDOW_H

#include <RcppArmadillo.h>

// Struct required for get_ofl_facetProperties()
const struct {
  arma::ivec8 iv_x1_dr { 0, -1, -1,  0,  0,  1, 1, 0};
  arma::ivec8 iv_x1_dc { 1,  0,  0, -1, -1,  0, 0, 1};

  arma::ivec8 iv_x2_dr {-1, -1, -1, -1,  1,  1, 1, 1};
  arma::ivec8 iv_x2_dc { 1,  1, -1, -1, -1, -1, 1, 1};
} fct_drdc; // DInf facets: focal window vector deltas of the row and column indices of the receiving cells x1 and x2

// Struct returned by get_ofl_facetProperties()
struct FacetProperties {
  bool ls_x1_oob {false}; // Is the row or column index of the receiving cell x1 out of bounds or the proportion for x1 == 0.0?
  double ns_p1 {NA_REAL}; // Proportion for x1
  arma::uword us_x1_r {}; // Row index of x1
  arma::uword us_x1_c {}; // Column index of x1

  bool ls_x2_oob {false}; // Is the row or column index of the receiving cell x2 out of bounds or the proportion for x2 == 0.0?
  double ns_p2 {NA_REAL}; // Proportion for x2
  arma::uword us_x2_r {}; // Row index of x2
  arma::uword us_x2_c {}; // Column index of x2
};

// Struct returned by get_ofl_x1x2()
template <typename T>
struct X1X2 {
  T x1 {}; // Value of the horizontal or vertical receiving cell x1
  T x2 {}; // Value of the diagonal receiving cell x2

  FacetProperties fct {};

  // Constructor
  X1X2(T NA_):
    x1 {NA_},
    x2 {NA_}
  {}
};

// Struct required for the *_ifl_* methods
const struct {
  arma::ivec8 iv_dr {-1, -1, -1,
                      0,      0,
                      1,  1,  1}; // Focal window vector deltas of the row indices
  arma::ivec8 iv_dc {-1,  0,  1,
                     -1,      1,
                     -1,  0,  1}; // Focal window vector deltas of the column indices

  arma::dvec8 nv_dir_min { 90.0, 135.0, 180.0,
                           45.0,        225.0,
                            0.0, 315.0, 270.0}; // Focal window vector lower bounds of the inflow directions
  arma::dvec8 nv_dir_mid {135.0, 180.0, 225.0,
                           90.0,        270.0,
                           45.0, 360.0, 315.0}; // Focal window vector midpoints of the inflow directions
  arma::dvec8 nv_dir_max {180.0, 225.0, 270.0,
                          135.0,        315.0,
                           90.0,  45.0, 360.0}; // Focal window vector upper bounds of the inflow directions

  arma::uvec3 uv_oob_lr {0, 1, 2}; // Focal window vector indices, which are out of bounds when row == 0
  arma::uvec3 uv_oob_ur {5, 6, 7}; // Focal window vector indices, which are out of bounds when row == max row
  arma::uvec3 uv_oob_lc {0, 3, 5}; // Focal window vector indices, which are out of bounds when col == 0
  arma::uvec3 uv_oob_uc {2, 4, 7}; // Focal window vector indices, which are out of bounds when col == max col

  arma::uvec3 uv_oob {0, 0, 0}; // Vector for setting focal window vector indices, which are out of bounds to 0
} ifl;

// Struct holding the transport calculation order
struct CalcOrder {
  std::vector<arma::uword> uv_r {}; // Vector holding row indices
  std::vector<arma::uword> uv_c {}; // vector holding column indices

  // Constructor
  CalcOrder(arma::uword n):
    uv_r {},
    uv_c {}
  {
    uv_r.reserve(n);
    uv_c.reserve(n);
  }
};

// Class handling all issues related to outflow and inflow by focusing on an
// examined cell and its eight neighbours
class FocalWindow {
public:
  const arma::uword us_rws {}; // Total number of rows of the extent of the river catchment as unsigned integer
  const arma::uword us_cls {}; // Total number of columns of the extent of the river catchment as unsigned integer
  const arma::sword is_rws {}; // Total number of rows as signed integer
  const arma::sword is_cls {}; // Total number of columns as signed integer

  // Constructor
  FocalWindow(arma::uword us_rws, arma::uword us_cls):
    us_rws {us_rws},
    us_cls {us_cls},
    is_rws {static_cast<arma::sword>(us_rws)},
    is_cls {static_cast<arma::sword>(us_cls)}
  {}

  // Methods
  FacetProperties get_ofl_facetProperties(
    const double ns_dir_inf,
    const arma::uword us_row,
    const arma::uword us_col
  );

  template <typename T>
  X1X2<T> get_ofl_x1x2(
    const FacetProperties& fct,
    const arma::Mat<T>& xm_xxx,
    const T NA_
  );

  double inc_ofl_x1x2(
    const X1X2<int>& x1x2,
    const double x1,
    const double x2,
    arma::dmat& nm_xxx
  );

  arma::dvec8 get_ifl_p(
    const arma::dmat& nm_dir_inf,
    const arma::uword us_row,
    const arma::uword us_col
  );

  template <typename T>
  arma::dvec8 get_ifl_x(
    const arma::dvec8& nv_ifl_p,
    const arma::uword us_row,
    const arma::uword us_col,
    const arma::Mat<T>& xm_xxx
  );
};

#endif

//' Determines the DInf facet of a cell and returns its properties
//'
//' First, the outflowing proportions are calculated from the DInf flow
//' direction. Then the row and column indices of the receiving cells x1 and x2
//' are determined. In case a receiving cell would lie outside of the extent of
//' the river catchment or has a receiving proportion of 0.0, it is flagged as
//' out of bounds.
//'
//' @param us_row The row index of the examined cell.
//' @param us_col The column index of the examined cell.
//' @param ns_dir_inf The DInf flow direction at the examined cell.
//'
//' @return A FacetProperties struct holding the DInf facet properties.
inline FacetProperties FocalWindow::get_ofl_facetProperties(
  const double ns_dir_inf,
  const arma::uword us_row,
  const arma::uword us_col
) {
  arma::uword us_fct {};
  FacetProperties fct {};

  // Determine the DInf facet and calculate its outflowing proportions
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
    Rcpp::stop("\"dir_inf\" out of range.");
  }

  // Determine the row and column indices of the receiving cells x1 and x2
  arma::sword is_row {static_cast<arma::sword>(us_row)};
  arma::sword is_col {static_cast<arma::sword>(us_col)};
  arma::sword is_x1_row {is_row + fct_drdc.iv_x1_dr[us_fct]};
  arma::sword is_x1_col {is_col + fct_drdc.iv_x1_dc[us_fct]};
  arma::sword is_x2_row {is_row + fct_drdc.iv_x2_dr[us_fct]};
  arma::sword is_x2_col {is_col + fct_drdc.iv_x2_dc[us_fct]};

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

//' Gets the values of the receiving cells x1 and x2
//'
//' In case a receiving cell is not out of bounds, its value is extracted from
//' the provided matrix. Furthermore, the DInf facet properties of the examined
//' cell are copied to the returned struct.
//'
//' @param fct The FacetProperties struct of the examined cell.
//' @param xm_xxx The matrix holding the values of the receiving cells.
//' @param NA_ The NA value corresponding to the matrix's data type.
//'
//' @return An X1X2 struct holding the values of the receiving cells x1 and x2
//'   and the DInf facet properties of the examined cell.
template <typename T>
inline X1X2<T> FocalWindow::get_ofl_x1x2(
  const FacetProperties& fct,
  const arma::Mat<T>& xm_xxx,
  const T NA_
) {
  X1X2<T> x1x2(NA_);

  if (!fct.ls_x1_oob) {
    x1x2.x1 = xm_xxx.at(fct.us_x1_r, fct.us_x1_c);
  }
  if (!fct.ls_x2_oob) {
    x1x2.x2 = xm_xxx.at(fct.us_x2_r, fct.us_x2_c);
  }

  x1x2.fct = fct;

  return x1x2;
}

//' Increases the existing values of the receiving cells x1 and x2
//'
//' In case a receiving cell is not out of bounds or NA_INTEGER in a conditional
//' layer, its existing value is increased by the respective provided value.
//'
//' @param x1x2 An X1X2<int> struct holding the values of the receiving cells x1
//'   and x2 of a conditional layer and the DInf facet properties of the
//'   examined cell.
//' @param x1 The value by which the existing value of the receiving cell x1 in
//'   nm_xxx shall be increased.
//' @param x2 The value by which the existing value of the receiving cell x2 in
//'   nm_xxx shall be increased.
//' @param nm_xxx The numeric matrix holding the values of the receiving cells
//'   x1 and x2, which shall be increased.
//'
//' @return The sum of the values by which the existing values of the receiving
//'   cells x1 and x2 were actually increased, i.e. the total outflowing load.
inline double FocalWindow::inc_ofl_x1x2(
  const X1X2<int>& x1x2,
  const double x1,
  const double x2,
  arma::dmat& nm_xxx
) {
  double ns_xxx {0.0};

  if (!Rcpp::IntegerMatrix::is_na(x1x2.x1)) {
    double ns_xxx_x1 {nm_xxx.at(x1x2.fct.us_x1_r, x1x2.fct.us_x1_c)};
    if (Rcpp::NumericMatrix::is_na(ns_xxx_x1)) {
      ns_xxx_x1 = 0.0;
    }

    nm_xxx.at(x1x2.fct.us_x1_r, x1x2.fct.us_x1_c) = ns_xxx_x1 + x1;
    ns_xxx += x1;
  }

  if (!Rcpp::IntegerMatrix::is_na(x1x2.x2)) {
    double ns_xxx_x2 {nm_xxx.at(x1x2.fct.us_x2_r, x1x2.fct.us_x2_c)};
    if (Rcpp::NumericMatrix::is_na(ns_xxx_x2)) {
      ns_xxx_x2 = 0.0;
    }

    nm_xxx.at(x1x2.fct.us_x2_r, x1x2.fct.us_x2_c) = ns_xxx_x2 + x2;
    ns_xxx += x2;
  }

  return ns_xxx;
}

//' Title
//'
//' Description
//'
//' @param
//'
//' @return
inline arma::dvec8 FocalWindow::get_ifl_p(
  const arma::dmat& nm_dir_inf,
  const arma::uword us_row,
  const arma::uword us_col
) {
  // Determine cells, which are out of bounds
  arma::uvec8 uv_cll(arma::fill::ones);

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

  // Determine proportions
  arma::dvec8 nv_ifl_p(arma::fill::zeros);

  for (arma::uword k = 0; k < uv_cll.n_elem; ++k) {
    if (uv_cll[k] == 1) {
      double ns_dir_inf {
        nm_dir_inf.at(us_row + ifl.iv_dr[k], us_col + ifl.iv_dc[k])
      };

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

//' Title
//'
//' Description
//'
//' @param
//'
//' @return
template <typename T>
inline arma::dvec8 FocalWindow::get_ifl_x(
  const arma::dvec8& nv_ifl_p,
  const arma::uword us_row,
  const arma::uword us_col,
  const arma::Mat<T>& xm_xxx
) {
  arma::dvec8 nv_ifl(arma::fill::zeros);

  for (arma::uword k = 0; k < nv_ifl_p.n_elem; ++k) {
    if (nv_ifl_p[k] > 0.0) {
      double ns_ifl {static_cast<double>(
        xm_xxx.at(us_row + ifl.iv_dr[k], us_col + ifl.iv_dc[k])
      )};

      if (ns_ifl > 0.0) {
        nv_ifl[k] = ns_ifl;
      }
    }
  }

  return nv_ifl;
}
