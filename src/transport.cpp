#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "MovingWindow.h"

// [[Rcpp::export]]
Rcpp::List transportCpp(
  const arma::dmat& nm_acc_inf,
  const arma::dmat& nm_dir_inf,
  const arma::dmat& nm_slp_inf,
  const arma::dmat& nm_man,
  const arma::imat& im_cha,
  const arma::imat& im_rds,
  const arma::imat& im_rip,
  const arma::imat& im_inl,
  const Rcpp::String substance,
  const Rcpp::S4& parameters,
  const Rcpp::S4& helpers,
  const int is_ths = 1
) {
  MovingWindow movingWindow{nm_dir_inf.n_rows, nm_dir_inf.n_cols};

  /* Transport calculation order
   * ===========================
   */

  /* Determine number of inflowing cells
   * -----------------------------------
   */

  arma::imat im_ifl(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_INTEGER)
  );

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 0; i < nm_dir_inf.n_rows; ++i) {
    for (arma::uword j = 0; j < nm_dir_inf.n_cols; ++j) {
      double ns_dir_inf{nm_dir_inf.at(i, j)};

      if (Rcpp::NumericMatrix::is_na(nm_acc_inf.at(i, j)) ||
          Rcpp::NumericMatrix::is_na(ns_dir_inf) ||
          ns_dir_inf == -1.0) {
        continue;
      }

      arma::dvec8 nv_ifl_p{movingWindow.get_ifl_p(nm_dir_inf, i, j)};
      im_ifl.at(i, j) = arma::accu(
        movingWindow.get_ifl_x<double>(nv_ifl_p, i, j, nm_acc_inf, NA_REAL) > 0.0
      );
    }
  }

  /* Determine order of rows and cols indices
   * ----------------------------------------
   */

  CalcOrder ord(arma::accu(im_ifl >= 0));

  for (arma::uword i = 0; i < im_ifl.n_rows; ++i) {
    for (arma::uword j = 0; j < im_ifl.n_cols; ++j) {
      if (im_ifl.at(i, j) == 0) {
        ord.uv_r.push_back(i);
        ord.uv_c.push_back(j);
      }
    }
  }

  // arma::imat im_ord(
  //   arma::size(nm_dir_inf),
  //   arma::fill::value(NA_INTEGER)
  // );

  FacetProperties fct{};
  arma::sword x1{}, x2{};

  for (arma::uword n = 0; n < ord.uv_r.capacity(); ++n) {
    if (n == ord.uv_r.size()) {
      Rcpp::stop("Warning: Could not determine a hydrologic consistent transport calculation order.");
    }

    fct = movingWindow.determineFacetProperties(
      nm_dir_inf.at(ord.uv_r[n], ord.uv_c[n]),
      ord.uv_r[n],
      ord.uv_c[n]
    );

    if (fct.ls_x1_oob) {
      x1 = NA_INTEGER;
    } else {
      x1 = im_ifl.at(fct.us_x1_r, fct.us_x1_c);
      if (!Rcpp::IntegerMatrix::is_na(x1)) {
        --x1;
        im_ifl.at(fct.us_x1_r, fct.us_x1_c) = x1;
      }
    }
    if (fct.ls_x2_oob) {
      x2 = NA_INTEGER;
    } else {
      x2 = im_ifl.at(fct.us_x2_r, fct.us_x2_c);
      if (!Rcpp::IntegerMatrix::is_na(x2)) {
        --x2;
        im_ifl.at(fct.us_x2_r, fct.us_x2_c) = x2;
      }
    }

    if (x1 == 0) {
      ord.uv_r.push_back(fct.us_x1_r);
      ord.uv_c.push_back(fct.us_x1_c);
    }
    if (x2 == 0) {
      ord.uv_r.push_back(fct.us_x2_r);
      ord.uv_c.push_back(fct.us_x2_c);
    }

    // im_ord.at(ord.uv_r[n], ord.uv_c[n]) = n;
  }

  /* Retentions and transports
   * =========================
   */

  /* Global variable declarations
   * ----------------------------
   */

  const double ns_rhy_a{parameters.slot("ns_rhy_a")};
  const double ns_rhy_b{parameters.slot("ns_rhy_b")};
  const double ns_tfc_inl{parameters.slot("ns_tfc_inl")};
  const double ns_cha_rto{parameters.slot("ns_cha_rto")};
  const double ns_man_rip{parameters.slot("ns_man_rip")};
  const double ns_man_cha{parameters.slot("ns_man_cha")};

  const Rcpp::NumericVector nv_enr_rto{parameters.slot("nv_enr_rto")};
  const double ns_enr_rto{nv_enr_rto[substance]};

  const double ns_dep_ovl_tmp{parameters.slot("ns_dep_ovl")};
  const double ns_dep_ovl{(substance == "SS") ?
    ns_dep_ovl_tmp : ns_dep_ovl_tmp / ns_enr_rto};
  const double ns_dep_cha{parameters.slot("ns_dep_cha")};

  const double ns_res{helpers.slot("ns_res")};
  const double ns_siz{helpers.slot("ns_siz")};

  // Flow path length of riparian zone cells
  const double ns_fpl_rip{(ns_res - ns_res * ns_cha_rto) / 2.0};

  int is_cha{}, is_rds{}, is_rip{}, is_inl{};
  double ns_dir_inf{}, ns_man{};

  double ns_rhy{}, ns_dir_inf_rad{}, ns_fpl{}, ns_v{}, ns_rtm{}, ns_rtm_rip{},
    ns_tfc_ifl{}, ns_tfc_lcl{}, ns_tfc_rip{};

  arma::uword i{}, j{};

  arma::dmat nm_xxr(
      arma::size(nm_dir_inf),
      arma::fill::value(NA_REAL)
  );
  arma::dmat nm_xxt(
      arma::size(nm_dir_inf),
      arma::fill::value(NA_REAL)
  );
  arma::dmat nm_xxt_out(
      arma::size(nm_dir_inf),
      arma::fill::value(NA_REAL)
  );

  /* Calculation of retentions and transports
   * ----------------------------------------
   */

  for (arma::uword n = 0; n < ord.uv_r.size(); ++n) {
    i = ord.uv_r[n];
    j = ord.uv_c[n];

    is_cha = im_cha.at(i, j);
    is_rds = im_rds.at(i, j);
    is_rip = im_rip.at(i, j);
    is_inl = im_inl.at(i, j);
    ns_dir_inf = nm_dir_inf.at(i, j);

    if (Rcpp::IntegerMatrix::is_na(is_cha)) {
      ns_man = nm_man.at(i, j);
    } else {
      ns_man = ns_man_cha;
    }

    // Hydraulic radius
    ns_rhy = ns_rhy_a * std::pow(nm_acc_inf.at(i, j) * ns_siz * 1e-6, ns_rhy_b);
    // Flow direction in radian
    ns_dir_inf_rad = ns_dir_inf * arma::datum::pi / 180.0;
    // Flow path length
    ns_fpl = ns_res * (std::abs(std::sin(ns_dir_inf_rad)) +
      std::abs(std::cos(ns_dir_inf_rad)));
    // Flow velocity
    ns_v = (1.0 / ns_man) * std::pow(ns_rhy, 2.0 / 3.0) *
      std::pow(nm_slp_inf.at(i, j) / 100.0, 1.0 / 2.0);
    // Residence time
    ns_rtm = ns_fpl / ns_v;
    // Residence time of riparian zone cell
    if (!Rcpp::IntegerMatrix::is_na(is_rip)) {
      ns_rtm_rip = ns_fpl_rip / ns_v;
    }

    // nv_ifl_p = movingWindow.get_ifl_p(nm_dir_inf, i, j);
    // im_ifl.at(i, j) = arma::accu(
    //   movingWindow.get_ifl_x<double>(nv_ifl_p, i, j, nm_acc_inf, NA_REAL) > 0.0
    // );
    // Transfer coefficients
    if (Rcpp::IntegerMatrix::is_na(is_cha)) {
      // Overland cell
      ns_tfc_ifl = std::exp(-ns_dep_ovl * ns_rtm);
      ns_tfc_lcl = std::exp(-ns_dep_ovl * ns_rtm * 0.5);

      // Riparian zone cell
      if (!Rcpp::IntegerMatrix::is_na(is_rip)) {
        ns_tfc_rip = std::exp(-ns_dep_ovl * ns_rtm_rip);
      }
    } else {
      // Channel cell
      ns_tfc_ifl = std::exp(-ns_dep_cha * ns_rtm);
      ns_tfc_lcl = std::exp(-ns_dep_cha * ns_rtm * 0.5);
    }
  }

  /* Cell loads and transfers
   * ========================
   */

  /* Global variable declarations
   * ----------------------------
   */

  for (arma::uword n = ord.uv_r.size() - 1; n >= 0; --n) {
    i = ord.uv_r[n];
    j = ord.uv_c[n];

    is_cha = im_cha.at(i, j);

    if (!Rcpp::IntegerMatrix::is_na(is_cha)) {
      continue;
    }
  }

//   /* Transport
//      =========
//   */
//   NumericMatrix nm_xxr_inp = na_real_matrix(is_rws, is_cls); // Retention of riparian zone and inlets (empty matrix to loop through)
//   NumericMatrix nm_xxt_inp = na_real_matrix(is_rws, is_cls); // Transport of riparian zone and inlets (empty matrix to loop through)

//   div_t code;
//
//   n = iv_ord_row.size();
//   for (int i = 0; i < n; ++i) { // Top-down calculation of retention and transport (overland and channels)
//     im_fDi_foc = im_fDi; // Moving inflow direction window
//     moving_data_window(iv_ord_row[i], iv_ord_col[i], is_rws, is_cls, is_row_min, is_row_max, is_col_min, is_col_max, im_fDi_foc);
//     im_foc     = im_dir(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving potential inflow direction window
//     nm_xxt_foc = nm_xxt(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving transport data window
//     if (IntegerMatrix::is_na(im_cha(iv_ord_row[i], iv_ord_col[i]))) { // Overland
//       nm_xxr(iv_ord_row[i], iv_ord_col[i]) = nm_xxe(iv_ord_row[i], iv_ord_col[i]) * (1.0 - nm_tfc_lcl(iv_ord_row[i], iv_ord_col[i])) + focal_sum_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc) * (1.0 - nm_tfc_tpt(iv_ord_row[i], iv_ord_col[i])); // Retention
//       nm_xxt(iv_ord_row[i], iv_ord_col[i]) = nm_xxe(iv_ord_row[i], iv_ord_col[i]) + focal_sum_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc) - nm_xxr(iv_ord_row[i], iv_ord_col[i]); // Transport
//       if (!IntegerMatrix::is_na(im_rip(iv_ord_row[i], iv_ord_col[i]))) { // Riparian zone
//         nm_xxr_inp(iv_ord_row[i], iv_ord_col[i]) = nm_xxt(iv_ord_row[i], iv_ord_col[i]) * (1.0 - nm_tfc_rip(iv_ord_row[i], iv_ord_col[i])); // Retention of riparian zone
//         nm_xxt_inp(iv_ord_row[i], iv_ord_col[i]) = nm_xxt(iv_ord_row[i], iv_ord_col[i]) - nm_xxr_inp(iv_ord_row[i], iv_ord_col[i]); // Transport of riparian zone
//       } else if (!IntegerMatrix::is_na(im_inl(iv_ord_row[i], iv_ord_col[i]))) { // Inlet
//         nm_xxr_inp(iv_ord_row[i], iv_ord_col[i]) = nm_xxt(iv_ord_row[i], iv_ord_col[i]) * (1.0 - ns_tfc_inl); // Retention of inlet
//         nm_xxt_inp(iv_ord_row[i], iv_ord_col[i]) = nm_xxt(iv_ord_row[i], iv_ord_col[i]) - nm_xxr_inp(iv_ord_row[i], iv_ord_col[i]); // Transport of inlet
//         code = div(im_inl(iv_ord_row[i], iv_ord_col[i]), is_cls);
//         is_row = code.quot - 1; // Outlet row number from integer code (C++ indices start at 0)
//         is_col = code.rem - 1; // Outlet column number from integer code (C++ indices start at 0)
//         if (NumericMatrix::is_na(nm_xxt_out(is_row, is_col))) { // Outlet load
//           nm_xxt_out(is_row, is_col) =  nm_xxt_inp(iv_ord_row[i], iv_ord_col[i]);
//         } else {
//           nm_xxt_out(is_row, is_col) += nm_xxt_inp(iv_ord_row[i], iv_ord_col[i]);
//         }
//       }
//     } else { // Channel
//       nm_xxt_inp_foc = nm_xxt_inp(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving transport of riparian zone and inlet data window
//       im_cha_foc     =     im_cha(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving channel data window
//       nm_xxr(iv_ord_row[i], iv_ord_col[i]) = focal_sum_matrix_cond(nm_xxt_inp_foc, im_fDi_foc, im_foc) * (1.0 - nm_tfc_lcl(iv_ord_row[i], iv_ord_col[i])) + focal_sum_two_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc, im_cha_foc) * (1.0 - nm_tfc_tpt(iv_ord_row[i], iv_ord_col[i]));
//       if (NumericMatrix::is_na(nm_xxt_out(iv_ord_row[i], iv_ord_col[i]))) {
//         nm_xxt(iv_ord_row[i], iv_ord_col[i]) = focal_sum_matrix_cond(nm_xxt_inp_foc, im_fDi_foc, im_foc) + focal_sum_two_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc, im_cha_foc) - nm_xxr(iv_ord_row[i], iv_ord_col[i]);
//       } else {
//         nm_xxt(iv_ord_row[i], iv_ord_col[i]) = focal_sum_matrix_cond(nm_xxt_inp_foc, im_fDi_foc, im_foc) + focal_sum_two_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc, im_cha_foc) + nm_xxt_out(iv_ord_row[i], iv_ord_col[i]) - nm_xxr(iv_ord_row[i], iv_ord_col[i]);
//       }
//     }
//   }
//
//   /* Cell loads & transfers
//      ----------------------
//   */
//   double ns_xxt_cld;
//   double ns_xxt_ctf;
//   NumericMatrix nm_xxt_ctf_foc;
//
//   NumericMatrix nm_xxe_net = matrix_minus_matrix_elementwise(nm_xxe, nm_xxr); // Net emission
//   NumericMatrix nm_xxt_cld = na_real_matrix(is_rws, is_cls); // Cell loads (empty matrix to loop through)
//   NumericMatrix nm_xxt_ctf = na_real_matrix(is_rws, is_cls); // Cell transfers (empty matrix to loop through)
//
//   n = iv_ord_ovl_row_rev.size();
//   for (int i = 0; i < n; ++i) { // Bottom-up calculation of cell loads and transfers (overland)
//     im_fDi_foc = im_fDi;
//     moving_data_window(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i], is_rws, is_cls, is_row_min, is_row_max, is_col_min, is_col_max, im_fDi_foc);
//     im_foc     = im_dir(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving potential inflow direction window
//     nm_xxt_foc = nm_xxt(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving transport data window
//     if (!NumericMatrix::is_na(nm_xxt_inp(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]))) { // Riparian zone or inlet cell
//       ns_xxt_cld = nm_xxt_inp(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]); // Intermediate cell load equals input into channel
//     } else if (!NumericMatrix::is_na(nm_xxt_ctf(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]))) { // Cell upstream of a riparian zone or inlet cell
//       ns_xxt_cld = nm_xxt_ctf(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]); // Intermediate cell load equals cell transfer from downstream
//     } else { // Other cell (cell at road embankment without subsurface drainage)
//       ns_xxt_cld = 0.0;
//     }
//     ns_xxt_ctf = ns_xxt_cld; // Intermediate cell transfer
//     ns_xxt_cld = nm_xxt(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]) * ns_xxt_cld / (nm_xxt(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]) + focal_sum_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc)); // Calculated intermediate cell load
//     if (nm_xxe_net(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]) < 0.0 || !std::isfinite(ns_xxt_cld)) { // Max possible cell load
//       ns_xxt_cld = 0.0;
//     } else {
//       ns_xxt_cld = std::min(ns_xxt_cld, nm_xxe_net(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]));
//     }
//     nm_xxt_cld(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]) = ns_xxt_cld; // Cell load
//     ns_xxt_ctf = std::max(ns_xxt_ctf - ns_xxt_cld, 0.0); // Cell load carry equals cell transfer
//     nm_xxt_ctf(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]) = ns_xxt_ctf; // Cell transfer
//     nm_xxt_ctf_foc = nm_xxt_ctf(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving cell transfer data window
//     focal_apportionment(ns_xxt_ctf, nm_xxt_ctf_foc, nm_xxt_foc, im_fDi_foc, im_foc); // Cell load/transfer apportionment
//     replace_submatrix(nm_xxt_ctf, nm_xxt_ctf_foc, is_row_min, is_col_min); // Inserts cell load/transfer apportionment into cell transfer
//   }

  return Rcpp::List::create(
    // Rcpp::Named("nm_xxr"    ) = nm_xxr    ,
    // Rcpp::Named("nm_xxt"    ) = nm_xxt    ,
    // Rcpp::Named("nm_xxt_inp") = nm_xxt_inp,
    // Rcpp::Named("nm_xxt_out") = nm_xxt_out,
    // Rcpp::Named("nm_xxt_cld") = nm_xxt_cld,
    // Rcpp::Named("nm_xxt_ctf") = nm_xxt_ctf
    Rcpp::Named("im_ifl") = im_ifl//,
    // Rcpp::Named("im_ord") = im_ord
  );
}
