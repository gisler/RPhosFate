#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "DInfWindow.h"

// [[Rcpp::export]]
arma::dmat DInfSlopeCpp(
  const arma::dmat& nm_dir_inf,
  const arma::dmat& nm_dem,
  const double ns_res,
  const int is_ths = 1
) {
  const double ns_res_dgl {std::sqrt(2.0 * ns_res * ns_res)};

  DinfWindow dinfWindow {nm_dir_inf.n_rows, nm_dir_inf.n_cols};

  arma::dmat nm_slp_inf(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_REAL)
  );

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 0; i < nm_dir_inf.n_rows; ++i) {
    for (arma::uword j = 0; j < nm_dir_inf.n_cols; ++j) {
      double ns_dir_inf {nm_dir_inf.at(i, j)};

      if (Rcpp::NumericMatrix::is_na(ns_dir_inf) || ns_dir_inf == -1.0) {
        continue;
      }

      FacetProperties fct {dinfWindow.get_ofl_facetProperties(ns_dir_inf, i, j)};
      X1X2<double> e1e2 {dinfWindow.get_ofl_x1x2<double>(fct, nm_dem, NA_REAL)};

      if (Rcpp::NumericMatrix::is_na(e1e2.x1) &&
          Rcpp::NumericMatrix::is_na(e1e2.x2)) {
        continue;
      }

      double ns_e0 {nm_dem.at(i, j)};

      if (Rcpp::NumericMatrix::is_na(e1e2.x1)) {
        nm_slp_inf.at(i, j) = (ns_e0 - e1e2.x2) / ns_res_dgl * 100.0;
        continue;
      }

      double ns_s1 {(ns_e0 - e1e2.x1) / ns_res};

      if (Rcpp::NumericMatrix::is_na(e1e2.x2)) {
        nm_slp_inf.at(i, j) = ns_s1 * 100.0;
        continue;
      }

      double ns_s2 {(e1e2.x1 - e1e2.x2) / ns_res};
      nm_slp_inf.at(i, j) = std::sqrt(ns_s1 * ns_s1 + ns_s2 * ns_s2) * 100.0;
    }
  }

  return nm_slp_inf;
}
