#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "MovingWindow.h"

// [[Rcpp::export]]
Rcpp::List ripInlCpp(
  const arma::dmat& nm_dir_inf,
  const arma::imat& im_cha,
  const arma::imat& im_rds,
  const int is_ths = 1
) {
  FocalWindow focalWindow {nm_dir_inf.n_rows, nm_dir_inf.n_cols};
  int is_rip {1}, is_inl {1};

  arma::imat im_rip(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_INTEGER)
  );
  arma::imat im_inl(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_INTEGER)
  );

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 0; i < nm_dir_inf.n_rows; ++i) {
    for (arma::uword j = 0; j < nm_dir_inf.n_cols; ++j) {
      double ns_dir_inf {nm_dir_inf.at(i, j)};

      if ( Rcpp::NumericMatrix::is_na(ns_dir_inf) ||
          !Rcpp::IntegerMatrix::is_na(im_cha.at(i, j)) ||
          !Rcpp::IntegerMatrix::is_na(im_rds.at(i, j)) ||
           ns_dir_inf == -1.0) {
        continue;
      }

      X1X2<int> cha1cha2 {focalWindow.get_x1x2<int>(ns_dir_inf, i, j, im_cha, NA_INTEGER)};
      X1X2<int> rds1rds2 {focalWindow.get_x1x2<int>(ns_dir_inf, i, j, im_rds, NA_INTEGER)};

      if (!Rcpp::IntegerMatrix::is_na(cha1cha2.x1) ||
          !Rcpp::IntegerMatrix::is_na(cha1cha2.x2)) {
        im_rip.at(i, j) = is_rip;
        # pragma omp atomic
        ++is_rip;
      }
      if ((rds1rds2.x1 == 1 && Rcpp::IntegerMatrix::is_na(cha1cha2.x1)) ||
          (rds1rds2.x2 == 1 && Rcpp::IntegerMatrix::is_na(cha1cha2.x2))) {
        im_inl.at(i, j) = is_inl;
        # pragma omp atomic
        ++is_inl;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("im_rip") = im_rip,
    Rcpp::Named("im_inl") = im_inl
  );
}
