#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::dmat D8slope(
  arma::imat& im_dir,
  arma::dmat& nm_dem,
  arma::imat& im_fDo,
  double ns_fpl,
  int is_ths = 1
) {
  /* With integers, missing values are stored as the smallest integer
   * (-2.147.483.648). See also https://adv-r.hadley.nz/rcpp.html
   */
  int NA_integer_ = Rcpp::IntegerVector::get_na();

  arma::dmat nm_dem_pad(
    arma::size(nm_dem) + 2,
    arma::fill::value(Rcpp::NumericVector::get_na())
  );
  nm_dem_pad(1, 1, arma::size(nm_dem)) = nm_dem;

  arma::uvec4 indices = {0, 2, 6, 8};
  arma::ivec4 iv_fDo_dgl = {im_fDo.elem(indices)};
  double ns_fpl_dgl = {std::sqrt(2.0 * std::pow(ns_fpl, 2.0))};

  arma::dmat nm_slp(
    arma::size(nm_dem),
    arma::fill::value(Rcpp::NumericVector::get_na())
  );

  double ns_tmp;
  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 1; i < nm_dem_pad.n_rows - 1; ++i) {
    for (arma::uword j = 1; j < nm_dem_pad.n_cols - 1; ++j) {
      if (im_dir.at(i - 1, j - 1) == NA_integer_) {
        continue;
      }

      ns_tmp = nm_dem_pad.at(i, j) - arma::conv_to<double>::from(
        nm_dem_pad.submat(i - 1, j - 1, i + 1, j + 1).eval().elem(
          arma::find(im_dir.at(i - 1, j - 1) == im_fDo)
        )
      );

      if (std::find(
        iv_fDo_dgl.begin(),
        iv_fDo_dgl.end(),
        im_dir.at(i - 1, j - 1)
      ) == iv_fDo_dgl.end()) {
        nm_slp.at(i - 1, j - 1) = ns_tmp / ns_fpl     * 100.0;
      } else {
        nm_slp.at(i - 1, j - 1) = ns_tmp / ns_fpl_dgl * 100.0;
      }
    }
  }

  return nm_slp;
}
