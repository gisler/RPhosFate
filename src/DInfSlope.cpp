#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::dmat DInfSlope(
  const arma::dmat& nm_dir_inf,
  const arma::dmat& nm_dem,
  const double ns_res,
  const int is_ths = 1
) {



  arma::dmat nm_slp_inf(
    arma::size(nm_dem),
    arma::fill::value(NA_REAL)
  );

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 0; i < nm_dir_inf.n_rows; ++i) {
    for (arma::uword j = 0; j < nm_dir_inf.n_cols; ++j) {
      // if (Rcpp::IntegerMatrix::is_na(im_dir.at(i - 1, j - 1))) {
      //   continue;
      // }
    }
  }

  return nm_slp_inf;
}
