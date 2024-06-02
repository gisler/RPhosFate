#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "MovingWindow.h"

// [[Rcpp::export]]
arma::dmat DInfSlope(
  const arma::dmat& nm_dir_inf,
  const arma::dmat& nm_dem,
  const double ns_res,
  const int is_ths = 1
) {
  MovingWindow movingWindow{
    static_cast<arma::sword>(nm_dem.n_rows),
    static_cast<arma::sword>(nm_dem.n_cols)
  };
  X1X2 e1e2{};
  double ns_e0{}, ns_s1{}, ns_s2{}, ns_s{};

  arma::dmat nm_slp_inf(
    arma::size(nm_dem),
    arma::fill::value(NA_REAL)
  );

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 0; i < nm_dem.n_rows; ++i) {
    for (arma::uword j = 0; j < nm_dem.n_cols; ++j) {
      if (Rcpp::NumericMatrix::is_na(nm_dir_inf.at(i, j))) {
        continue;
      }

      ns_e0 = nm_dem.at(i, j);
      e1e2 = movingWindow.determine_x1x2(nm_dir_inf.at(i, j), i, j, nm_dem);
      ns_s1 = (     ns_e0 - e1e2.ns_x1) / ns_res;
      ns_s2 = (e1e2.ns_x1 - e1e2.ns_x2) / ns_res;
      ns_s = std::sqrt(ns_s1 * ns_s1 + ns_s2 * ns_s2);

      nm_slp_inf.at(i, j) = ns_s;
    }
  }

  return nm_slp_inf;
}
