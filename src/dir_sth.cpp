#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::imat dir_sth(
  const arma::imat& im_dir,
  const arma::imat& im_sth,
  const arma::imat& im_fDo,
  const int is_ths = 1
) {
  /* With integers, missing values are stored as the smallest integer
   * (-2.147.483.648). See also https://adv-r.hadley.nz/rcpp.html
   */
  const int NA_integer_ = Rcpp::IntegerVector::get_na();

  arma::imat im_sth_pad(arma::size(im_sth) + 2, arma::fill::value(NA_integer_));
  im_sth_pad(1, 1, arma::size(im_sth)) = im_sth;

  int is_sth = {1};
  arma::imat im_xxx(arma::size(im_sth), arma::fill::value(NA_integer_));

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 1; i < im_sth_pad.n_rows - 1; ++i) {
    for (arma::uword j = 1; j < im_sth_pad.n_cols - 1; ++j) {
      if (im_dir.at(i - 1, j - 1) == NA_integer_ ||
          im_sth.at(i - 1, j - 1) != NA_integer_) {
        continue;
      }

      int is_tmp;
      is_tmp = arma::conv_to<int>::from(
        im_sth_pad.submat(i - 1, j - 1, i + 1, j + 1).eval().elem(
          arma::find(im_dir.at(i - 1, j - 1) == im_fDo)
        )
      );

      if (is_tmp != NA_integer_ && is_tmp != 0) {
        im_xxx.at(i - 1, j - 1) = is_sth;
        ++is_sth;
      }
    }
  }

  return im_xxx;
}
