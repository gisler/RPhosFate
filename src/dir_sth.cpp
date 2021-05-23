#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::export]]
arma::imat dir_sth(
  arma::imat& im_dir,
  arma::imat& im_sth,
  arma::imat& im_fDo,
  int is_ths = 1
) {
  /* With integers, missing values are stored as the smallest integer
   * (-2.147.483.648). See also https://adv-r.hadley.nz/rcpp.html
   */
  int NA_integer_ = Rcpp::IntegerVector::get_na();

  arma::imat im_sth_pad(arma::size(im_sth) + 2);
  im_sth_pad.fill(NA_integer_);
  im_sth_pad(1, 1, arma::size(im_sth)) = im_sth;

  int x = {1};
  arma::imat im_xxx(arma::size(im_sth));
  im_xxx.fill(NA_integer_);

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 1; i < im_sth_pad.n_rows - 1; ++i) {
    for (arma::uword j = 1; j < im_sth_pad.n_cols - 1; ++j) {
      if (im_dir.at(i - 1, j - 1) == NA_integer_ ||
          im_sth.at(i - 1, j - 1) != NA_integer_) {
        continue;
      }

      arma::imat33 im_tmp;
      int tmp;
      im_tmp = im_sth_pad.submat(i - 1, j - 1, i + 1, j + 1);
      tmp = arma::conv_to<int>::from(
        im_tmp.elem(arma::find(im_dir.at(i - 1, j - 1) == im_fDo))
      );
      if (tmp != NA_integer_ && tmp != 0) {
        im_xxx.at(i - 1, j - 1) = x;
        ++x;
      }
    }
  }

  return im_xxx;
}
