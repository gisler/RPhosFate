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
  MovingWindow movingWindow{nm_dem.n_rows, nm_dem.n_cols};
  X1X2 e1e2{};

  arma::dmat nm_slp_inf(
    arma::size(nm_dem),
    arma::fill::value(NA_REAL)
  );

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 0; i < movingWindow.is_rws; ++i) {
    for (arma::uword j = 0; j < movingWindow.is_cls; ++j) {
      if (Rcpp::NumericMatrix::is_na(nm_dir_inf.at(i, j))) {
        continue;
      }

      e1e2 = movingWindow.determine_x1x2(nm_dir_inf.at(i, j), i, j, nm_dem);
    }
  }

  return nm_slp_inf;
}
