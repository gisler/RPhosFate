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
  arma::imat im_sth_pad(arma::size(im_sth) + 2, arma::fill::value(NA_INTEGER));
  im_sth_pad(1, 1, arma::size(im_sth)) = im_sth;

  int is_sth = {1};
  arma::imat im_xxx(arma::size(im_sth), arma::fill::value(NA_INTEGER));

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 1; i < im_sth_pad.n_rows - 1; ++i) {
    for (arma::uword j = 1; j < im_sth_pad.n_cols - 1; ++j) {
      if ( Rcpp::IntegerMatrix::is_na(im_dir.at(i - 1, j - 1)) ||
          !Rcpp::IntegerMatrix::is_na(im_sth.at(i - 1, j - 1))) {
        continue;
      }

      int is_tmp;
      is_tmp = arma::conv_to<int>::from(
        im_sth_pad.submat(i - 1, j - 1, i + 1, j + 1).eval().elem(
          arma::find(im_dir.at(i - 1, j - 1) == im_fDo)
        )
      );

      if (!Rcpp::IntegerMatrix::is_na(is_tmp) && is_tmp != 0) {
        im_xxx.at(i - 1, j - 1) = is_sth;
        ++is_sth;
      }
    }
  }

  return im_xxx;
}
