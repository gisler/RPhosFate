#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "MovingWindow.h"

// [[Rcpp::export]]
Rcpp::List transportCpp(
  arma::dmat& nm_acc_inf,
  arma::dmat& nm_dir_inf,
  arma::imat& im_cha,
  const int is_ths = 1
) {
  MovingWindow movingWindow{nm_dir_inf.n_rows, nm_dir_inf.n_cols};

  /* Transport calculation order
     ===========================
  */

  // Determine number of inflowing cells
  arma::imat im_ifl(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_INTEGER)
  );

  #pragma omp parallel for num_threads(is_ths) collapse(2)
  for (arma::uword i = 0; i < nm_dir_inf.n_rows; ++i) {
    for (arma::uword j = 0; j < nm_dir_inf.n_cols; ++j) {
      double ns_dir_inf{};
      arma::dvec8 nv_ifl_p{};

      ns_dir_inf = nm_dir_inf.at(i, j);

      if (Rcpp::NumericMatrix::is_na(nm_acc_inf.at(i, j)) ||
          Rcpp::NumericMatrix::is_na(ns_dir_inf) ||
          ns_dir_inf == -1.0) {
        continue;
      }

      nv_ifl_p = movingWindow.get_ifl_p(nm_dir_inf, i, j);
      im_ifl.at(i, j) = arma::accu(
        movingWindow.get_ifl<double>(nv_ifl_p, i, j, nm_acc_inf, NA_REAL) > 0.0
      );
    }
  }

  // Determine order of rows and cols indices
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
      Rcpp::Rcerr << "Warning: Could not determine a hydrologic consistent transport calculation order." << std::endl;
      continue;
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

  /* Transport
     =========
  */


//   double ns_cha_rto = parameters.slot("ns_cha_rto");
//   double ns_man_rip = parameters.slot("ns_man_rip");
//   double ns_man_cha = parameters.slot("ns_man_cha");
//   double ns_dep_cha = parameters.slot("ns_dep_cha");
//   double ns_res = helpers.slot("ns_res");
//   int is_rws = helpers.slot("is_rws");
//   int is_cls = helpers.slot("is_cls");
//   IntegerMatrix im_fDo = helpers.slot("im_fDo");
//   IntegerMatrix im_fDi = helpers.slot("im_fDi");
//   IntegerVector iv_ord_row = order.slot("iv_ord_row");
//   IntegerVector iv_ord_col = order.slot("iv_ord_col");
//   IntegerVector iv_ord_ovl_row_rev = order.slot("iv_ord_ovl_row_rev");
//   IntegerVector iv_ord_ovl_col_rev = order.slot("iv_ord_ovl_col_rev");
//
//   int n;
//
//   /* Preparations
//      ============
//   */
//   NumericMatrix nm_man_pcd = replace_matrix_na_cond(nm_man, im_cha, ns_man_cha); // Manning n combined with channel manning n
//   // Add manning n considering preferential flow path manning n here
//   IntegerVector iv_fDo_dgl = IntegerVector::create(im_fDo[0], im_fDo[2], im_fDo[6], im_fDo[8]); // Diagonal outflow direction vector
//
//   NumericMatrix nm_rtm     = na_real_matrix(is_rws, is_cls); // Residence time (empty matrix to loop through)
//   NumericMatrix nm_rtm_rip = na_real_matrix(is_rws, is_cls); // Riparian zone residence time (empty matrix to loop through)
//   NumericMatrix nm_tfc_lcl = na_real_matrix(is_rws, is_cls); // Local transfer coefficient (empty matrix to loop through)
//   NumericMatrix nm_tfc_tpt = na_real_matrix(is_rws, is_cls); // Transport transfer coefficient (empty matrix to loop through)
//   NumericMatrix nm_tfc_rip = na_real_matrix(is_rws, is_cls); // Riparian zone transfer coefficient (empty matrix to loop through)
//
//   double ns_fpl     = ns_res; // Flow path length
//   double ns_fpl_dgl = sqrt(2.0 * pow(ns_fpl, 2.0)); // Diagonal flow path length
//   double ns_fpl_rip = (ns_fpl - ns_fpl * ns_cha_rto) / 2.0; // Riparian zone flow path length
//
//   n = is_rws * is_cls;
//   for (int i = 0; i < n; ++i) {
//     if (std::find(iv_fDo_dgl.begin(), iv_fDo_dgl.end(), im_dir[i]) == iv_fDo_dgl.end()) { // Residence time
//       nm_rtm[i] = ns_fpl / ((1.0 / nm_man_pcd[i]) * pow(nm_rhy[i], 2.0 / 3.0) * pow(nm_slp[i] / 100.0, 1.0 / 2.0));
//       if (!IntegerMatrix::is_na(im_rip[i])) {
//         nm_rtm_rip[i] = ns_fpl_rip / ((1.0 / ns_man_rip) * pow(nm_rhy[i], 2.0 / 3.0) * pow(nm_slp[i] / 100.0, 1.0 / 2.0));
//       }
//     } else { // Residence time of diagonal flow paths
//       nm_rtm[i] = ns_fpl_dgl / ((1.0 / nm_man_pcd[i]) * pow(nm_rhy[i], 2.0 / 3.0) * pow(nm_slp[i] / 100.0, 1.0 / 2.0));
//       if (!IntegerMatrix::is_na(im_rip[i])) {
//         nm_rtm_rip[i] = ns_fpl_rip / ((1.0 / ns_man_rip) * pow(nm_rhy[i], 2.0 / 3.0) * pow(nm_slp[i] / 100.0, 1.0 / 2.0));
//       }
//     }
//     if (IntegerMatrix::is_na(im_cha[i])) { // Overland
//       nm_tfc_lcl[i] = exp(-ns_dep_ovl * nm_rtm[i] * 0.5);
//       nm_tfc_tpt[i] = exp(-ns_dep_ovl * nm_rtm[i]);
//       if (!IntegerMatrix::is_na(im_rip[i])) { // Riparian zone
//         nm_tfc_rip[i] = exp(-ns_dep_ovl * nm_rtm_rip[i]);
//       }
//     } else { // Channel
//       nm_tfc_lcl[i] = exp(-ns_dep_cha * nm_rtm[i] * 0.5);
//       nm_tfc_tpt[i] = exp(-ns_dep_cha * nm_rtm[i]);
//     }
//   }
//
//   /* Transport
//      =========
//   */
//   NumericMatrix nm_xxr     = na_real_matrix(is_rws, is_cls); // Retention (empty matrix to loop through)
//   NumericMatrix nm_xxt     = na_real_matrix(is_rws, is_cls); // Transport (empty matrix to loop through)
//   NumericMatrix nm_xxr_inp = na_real_matrix(is_rws, is_cls); // Retention of riparian zone and inlets (empty matrix to loop through)
//   NumericMatrix nm_xxt_inp = na_real_matrix(is_rws, is_cls); // Transport of riparian zone and inlets (empty matrix to loop through)
//   NumericMatrix nm_xxt_out = na_real_matrix(is_rws, is_cls); // Outlet loads (empty matrix to loop through)
//   NumericMatrix nm_xxt_cld = na_real_matrix(is_rws, is_cls); // Cell loads (empty matrix to loop through)
//   NumericMatrix nm_xxt_ctf = na_real_matrix(is_rws, is_cls); // Cell transfers (empty matrix to loop through)
//
//   IntegerMatrix im_fDi_foc;
//   int is_row_min;
//   int is_row_max;
//   int is_col_min;
//   int is_col_max;
//   IntegerMatrix im_foc;
//   NumericMatrix nm_xxt_foc;
//
//   /* Transport
//      ---------
//   */
//   div_t code;
//   int is_row;
//   int is_col;
//   NumericMatrix nm_xxt_inp_foc;
//   IntegerMatrix im_cha_foc;
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
