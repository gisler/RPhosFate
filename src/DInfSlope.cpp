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
  const double ns_res_dgl{std::sqrt(2.0 * ns_res * ns_res)};

  double ns_dir_inf{}, ns_e0{}, ns_s1{}, ns_s2{};
  MovingWindow movingWindow{
    static_cast<arma::sword>(nm_dem.n_rows),
    static_cast<arma::sword>(nm_dem.n_cols)
  };
  X1X2<double> e1e2{};

  arma::dmat nm_slp_inf(
    arma::size(nm_dem),
    arma::fill::value(NA_REAL)
  );

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 0; i < nm_dem.n_rows; ++i) {
    for (arma::uword j = 0; j < nm_dem.n_cols; ++j) {
      ns_dir_inf = nm_dir_inf.at(i, j);

      if (Rcpp::NumericMatrix::is_na(ns_dir_inf) || ns_dir_inf == -1.0) {
        continue;
      }

      e1e2 = movingWindow.determine_x1x2<double>(ns_dir_inf, i, j, nm_dem, NA_REAL);

      if (Rcpp::NumericMatrix::is_na(e1e2.x1) &&
          Rcpp::NumericMatrix::is_na(e1e2.x2)) {
        continue;
      }

      ns_e0 = nm_dem.at(i, j);

      if (std::set<double>{45.0, 135.0, 225.0, 315.0}.count(ns_dir_inf) > 0 ||
          Rcpp::NumericMatrix::is_na(e1e2.x1)) {
        nm_slp_inf.at(i, j) = (ns_e0 - e1e2.x2) / ns_res_dgl;
        continue;
      }

      ns_s1 = (ns_e0 - e1e2.x1) / ns_res;

      if (std::set<double>{0.0, 90.0, 180.0, 270.0, 360.0}.count(ns_dir_inf) > 0 ||
          Rcpp::NumericMatrix::is_na(e1e2.x2)) {
        nm_slp_inf.at(i, j) = ns_s1;
        continue;
      }

      ns_s2 = (e1e2.x1 - e1e2.x2) / ns_res;
      nm_slp_inf.at(i, j) = std::sqrt(ns_s1 * ns_s1 + ns_s2 * ns_s2);
    }
  }

  return nm_slp_inf;
}
