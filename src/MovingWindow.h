#ifndef MOVINGWINDOW_H
#define MOVINGWINDOW_H

#include <RcppArmadillo.h>

// Facet elevation modifiers
const struct {
  arma::ivec8 iv_dr_x1{ 0, -1, -1,  0,  0,  1, 1, 0};
  arma::ivec8 iv_dc_x1{ 1,  0,  0, -1, -1,  0, 0, 1};

  arma::ivec8 iv_dr_x2{-1, -1, -1, -1,  1,  1, 1, 1};
  arma::ivec8 iv_dc_x2{ 1,  1, -1, -1, -1, -1, 1, 1};
} drdc;

struct X1X2 {
  double x1{};
  double x2{};
};

class MovingWindow {
public:
  const arma::uword is_rws{};
  const arma::uword is_cls{};

  arma::uword determineFacet(const double ns_dir_inf);

  X1X2 determine_x1x2(
    const double ns_dir_inf,
    const arma::uword i,
    const arma::uword j,
    const arma::dmat& nm_xxx
  );
};

#endif
