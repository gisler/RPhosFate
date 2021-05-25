#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix na_real_matrix(int rws, int cls) {
  NumericMatrix nm(rws, cls);
  int n = rws * cls;
  for (int i = 0; i < n; ++i) {
    nm[i] = NA_REAL;
  }
  return nm;
}

NumericMatrix replace_matrix_na_cond(NumericMatrix nm_values, IntegerMatrix im_cond, double ns_replace) {
  int rws = nm_values.nrow();
  int cls = nm_values.ncol();
  NumericMatrix nm = na_real_matrix(rws, cls);
  int n = rws * cls;
  for (int i = 0; i < n; ++i) {
    if (IntegerMatrix::is_na(im_cond[i])) {
      nm[i] = nm_values[i];
    } else {
      nm[i] = ns_replace;
    }
  }
  return nm;
}

void moving_data_window(int is_ord_row, int is_ord_col, int is_rws, int is_cls, int &is_row_min, int &is_row_max, int &is_col_min, int &is_col_max, IntegerMatrix &im_fDx_foc) {
  if (is_ord_row == 0) { // Left boundary of moving data window
    is_row_min = 0;
    im_fDx_foc = im_fDx_foc(Range(1, im_fDx_foc.nrow() - 1), Range(0, im_fDx_foc.ncol() - 1)); // Adjusts moving inflow direction window
  } else {
    is_row_min = is_ord_row - 1;
  }
  if (is_ord_row == is_rws - 1) { // Right boundary of moving data window
    is_row_max = is_rws - 1;
    im_fDx_foc = im_fDx_foc(Range(0, im_fDx_foc.nrow() - 2), Range(0, im_fDx_foc.ncol() - 1));
  } else {
    is_row_max = is_ord_row + 1;
  }
  if (is_ord_col == 0) { // Bottom boundary of moving data window
    is_col_min = 0;
    im_fDx_foc = im_fDx_foc(Range(0, im_fDx_foc.nrow() - 1), Range(1, im_fDx_foc.ncol() - 1));
  } else {
    is_col_min = is_ord_col - 1;
  }
  if (is_ord_col == is_cls - 1) { // Top boundary of moving data window
    is_col_max = is_cls - 1;
    im_fDx_foc = im_fDx_foc(Range(0, im_fDx_foc.nrow() - 1), Range(0, im_fDx_foc.ncol() - 2));
  } else {
    is_col_max = is_ord_col + 1;
  }
}

double focal_sum_int_cond(NumericMatrix nm_sum, IntegerMatrix im_foc, int is_cond) {
  double total = 0.0;
  int n = nm_sum.nrow() * nm_sum.ncol();
  for (int i = 0; i < n; ++i) {
    if (!NumericMatrix::is_na(nm_sum[i]) && im_foc[i] == is_cond) {
      total += nm_sum[i];
    }
  }
  return total;
}

double focal_sum_matrix_cond(NumericMatrix nm_sum, IntegerMatrix im_foc, IntegerMatrix im_cond) {
  double total = 0.0;
  int n = nm_sum.nrow() * nm_sum.ncol();
  for (int i = 0; i < n; ++i) {
    if (!NumericMatrix::is_na(nm_sum[i]) && im_foc[i] == im_cond[i]) {
      total += nm_sum[i];
    }
  }
  return total;
}

double focal_sum_two_matrix_cond(NumericMatrix nm_sum, IntegerMatrix im_foc, IntegerMatrix im_cond_a, IntegerMatrix im_cond_b) {
  double total = 0.0;
  int n = nm_sum.nrow() * nm_sum.ncol();
  for (int i = 0; i < n; ++i) {
    if (!NumericMatrix::is_na(nm_sum[i]) && im_foc[i] == im_cond_a[i] && !IntegerMatrix::is_na(im_cond_b[i])) {
      total += nm_sum[i];
    }
  }
  return total;
}

//NumericMatrix matrix_times_matrix_elementwise(NumericMatrix nm_x, NumericMatrix nm_y) {
//  int rws = nm_x.nrow();
//  int cls = nm_x.ncol();
//  NumericMatrix nm = na_real_matrix(rws, cls);
//  int n = rws * cls;
//  for (int i = 0; i < n; ++i) {
//    nm[i] = nm_x[i] * nm_y[i];
//  }
//  return nm;
//}

NumericMatrix matrix_minus_matrix_elementwise(NumericMatrix nm_x, NumericMatrix nm_y) {
  int rws = nm_x.nrow();
  int cls = nm_x.ncol();
  NumericMatrix nm = na_real_matrix(rws, cls);
  int n = rws * cls;
  for (int i = 0; i < n; ++i) {
    nm[i] = nm_x[i] - nm_y[i];
  }
  return nm;
}

void focal_apportionment(double ns_app, NumericMatrix nm_app, NumericMatrix nm_wgt, IntegerMatrix im_foc, IntegerMatrix im_cond) {
  double factor = ns_app / focal_sum_matrix_cond(nm_wgt, im_foc, im_cond);
  if (!std::isfinite(factor)) {
    factor = 0.0;
  }
  int n = nm_app.nrow() * nm_app.ncol();
  for (int i = 0; i < n; ++i) {
    if (!NumericMatrix::is_na(nm_wgt[i]) && im_foc[i] == im_cond[i]) {
      nm_app[i] = nm_wgt[i] * factor;
    }
  }
}

void replace_submatrix(NumericMatrix nm_values, NumericMatrix nm_rpl, int is_row_sub, int is_col_sub) {
  int rws = nm_rpl.nrow();
  int cls = nm_rpl.ncol();
  for (int i = 0; i < rws; ++i) {
    for (int j = 0; j < cls; ++j) {
      nm_values(is_row_sub + i, is_col_sub + j) = nm_rpl(i, j);
    }
  }
}

// [[Rcpp::export]]
List transportCpp(
  S4 parameters,
  double ns_dep_ovl,
  double ns_tfc_inl,
  S4 helpers,
  S4 order,
  IntegerMatrix im_cha,
  IntegerMatrix im_dir,
  IntegerMatrix im_inl,
  IntegerMatrix im_rip,
  NumericMatrix nm_man,
  NumericMatrix nm_xxe,
  NumericMatrix nm_rhy,
  NumericMatrix nm_slp
) {
  double ns_cha_rto = parameters.slot("ns_cha_rto");
  double ns_man_rip = parameters.slot("ns_man_rip");
  double ns_man_cha = parameters.slot("ns_man_cha");
  double ns_dep_cha = parameters.slot("ns_dep_cha");
  int is_res = helpers.slot("is_res");
  int is_rws = helpers.slot("is_rws");
  int is_cls = helpers.slot("is_cls");
  IntegerMatrix im_fDo = helpers.slot("im_fDo");
  IntegerMatrix im_fDi = helpers.slot("im_fDi");
  IntegerVector iv_ord_row = order.slot("iv_ord_row");
  IntegerVector iv_ord_col = order.slot("iv_ord_col");
  IntegerVector iv_ord_ovl_row_rev = order.slot("iv_ord_ovl_row_rev");
  IntegerVector iv_ord_ovl_col_rev = order.slot("iv_ord_ovl_col_rev");

  int n;

  /* Preparations
     ============
  */
  NumericMatrix nm_man_pcd = replace_matrix_na_cond(nm_man, im_cha, ns_man_cha); // Manning n combined with channel manning n
  // Add manning n considering preferential flow path manning n here
  IntegerVector iv_fDo_dgl = IntegerVector::create(im_fDo[0], im_fDo[2], im_fDo[6], im_fDo[8]); // Diagonal outflow direction vector

  NumericMatrix nm_rtm     = na_real_matrix(is_rws, is_cls); // Residence time (empty matrix to loop through)
  NumericMatrix nm_rtm_rip = na_real_matrix(is_rws, is_cls); // Riparian zone residence time (empty matrix to loop through)
  NumericMatrix nm_tfc_lcl = na_real_matrix(is_rws, is_cls); // Local transfer coefficient (empty matrix to loop through)
  NumericMatrix nm_tfc_tpt = na_real_matrix(is_rws, is_cls); // Transport transfer coefficient (empty matrix to loop through)
  NumericMatrix nm_tfc_rip = na_real_matrix(is_rws, is_cls); // Riparian zone transfer coefficient (empty matrix to loop through)

  double ns_fpl     = static_cast<double>(is_res); // Flow path length
  double ns_fpl_dgl = sqrt(2.0 * pow(ns_fpl, 2.0)); // Diagonal flow path length
  double ns_fpl_rip = (ns_fpl - ns_fpl * ns_cha_rto) / 2.0; // Riparian zone flow path length

  n = is_rws * is_cls;
  for (int i = 0; i < n; ++i) {
    if (std::find(iv_fDo_dgl.begin(), iv_fDo_dgl.end(), im_dir[i]) == iv_fDo_dgl.end()) { // Residence time
      nm_rtm[i] = ns_fpl / ((1.0 / nm_man_pcd[i]) * pow(nm_rhy[i], 2.0 / 3.0) * pow(nm_slp[i] / 100.0, 1.0 / 2.0));
      if (!IntegerMatrix::is_na(im_rip[i])) {
        nm_rtm_rip[i] = ns_fpl_rip / ((1.0 / ns_man_rip) * pow(nm_rhy[i], 2.0 / 3.0) * pow(nm_slp[i] / 100.0, 1.0 / 2.0));
      }
    } else { // Residence time of diagonal flow paths
      nm_rtm[i] = ns_fpl_dgl / ((1.0 / nm_man_pcd[i]) * pow(nm_rhy[i], 2.0 / 3.0) * pow(nm_slp[i] / 100.0, 1.0 / 2.0));
      if (!IntegerMatrix::is_na(im_rip[i])) {
        nm_rtm_rip[i] = ns_fpl_rip / ((1.0 / ns_man_rip) * pow(nm_rhy[i], 2.0 / 3.0) * pow(nm_slp[i] / 100.0, 1.0 / 2.0));
      }
    }
    if (IntegerMatrix::is_na(im_cha[i])) { // Overland
      nm_tfc_lcl[i] = exp(-ns_dep_ovl * nm_rtm[i] * 0.5);
      nm_tfc_tpt[i] = exp(-ns_dep_ovl * nm_rtm[i]);
      if (!IntegerMatrix::is_na(im_rip[i])) { // Riparian zone
        nm_tfc_rip[i] = exp(-ns_dep_ovl * nm_rtm_rip[i]);
      }
    } else { // Channel
      nm_tfc_lcl[i] = exp(-ns_dep_cha * nm_rtm[i] * 0.5);
      nm_tfc_tpt[i] = exp(-ns_dep_cha * nm_rtm[i]);
    }
  }

  /* Transport
     =========
  */
  NumericMatrix nm_xxr     = na_real_matrix(is_rws, is_cls); // Retention (empty matrix to loop through)
  NumericMatrix nm_xxt     = na_real_matrix(is_rws, is_cls); // Transport (empty matrix to loop through)
  NumericMatrix nm_xxr_inp = na_real_matrix(is_rws, is_cls); // Retention of riparian zone and inlets (empty matrix to loop through)
  NumericMatrix nm_xxt_inp = na_real_matrix(is_rws, is_cls); // Transport of riparian zone and inlets (empty matrix to loop through)
  NumericMatrix nm_xxt_out = na_real_matrix(is_rws, is_cls); // Outlet loads (empty matrix to loop through)
  NumericMatrix nm_xxt_cld = na_real_matrix(is_rws, is_cls); // Cell loads (empty matrix to loop through)
  NumericMatrix nm_xxt_ctf = na_real_matrix(is_rws, is_cls); // Cell transfers (empty matrix to loop through)

  IntegerMatrix im_fDi_foc;
  int is_row_min;
  int is_row_max;
  int is_col_min;
  int is_col_max;
  IntegerMatrix im_foc;
  NumericMatrix nm_xxt_foc;

  /* Transport
     ---------
  */
  div_t code;
  int is_row;
  int is_col;
  NumericMatrix nm_xxt_inp_foc;
  IntegerMatrix im_cha_foc;

  n = iv_ord_row.size();
  for (int i = 0; i < n; ++i) { // Top-down calculation of retention and transport (overland and channels)
    im_fDi_foc = im_fDi; // Moving inflow direction window
    moving_data_window(iv_ord_row[i], iv_ord_col[i], is_rws, is_cls, is_row_min, is_row_max, is_col_min, is_col_max, im_fDi_foc);
    im_foc     = im_dir(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving potential inflow direction window
    nm_xxt_foc = nm_xxt(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving transport data window
    if (IntegerMatrix::is_na(im_cha(iv_ord_row[i], iv_ord_col[i]))) { // Overland
      nm_xxr(iv_ord_row[i], iv_ord_col[i]) = nm_xxe(iv_ord_row[i], iv_ord_col[i]) * (1.0 - nm_tfc_lcl(iv_ord_row[i], iv_ord_col[i])) + focal_sum_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc) * (1.0 - nm_tfc_tpt(iv_ord_row[i], iv_ord_col[i])); // Retention
      nm_xxt(iv_ord_row[i], iv_ord_col[i]) = nm_xxe(iv_ord_row[i], iv_ord_col[i]) + focal_sum_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc) - nm_xxr(iv_ord_row[i], iv_ord_col[i]); // Transport
      if (!IntegerMatrix::is_na(im_inl(iv_ord_row[i], iv_ord_col[i]))) { // Inlet
        nm_xxr_inp(iv_ord_row[i], iv_ord_col[i]) = nm_xxt(iv_ord_row[i], iv_ord_col[i]) * (1.0 - ns_tfc_inl); // Retention of inlet
        nm_xxt_inp(iv_ord_row[i], iv_ord_col[i]) = nm_xxt(iv_ord_row[i], iv_ord_col[i]) - nm_xxr_inp(iv_ord_row[i], iv_ord_col[i]); // Transport of inlet
        code = div(im_inl(iv_ord_row[i], iv_ord_col[i]), is_cls);
        is_row = code.quot - 1; // Outlet row number from integer code (C++ indices start at 0)
        is_col = code.rem - 1; // Outlet column number from integer code (C++ indices start at 0)
        if (NumericMatrix::is_na(nm_xxt_out(is_row, is_col))) { // Outlet load
          nm_xxt_out(is_row, is_col) =  nm_xxt_inp(iv_ord_row[i], iv_ord_col[i]);
        } else {
          nm_xxt_out(is_row, is_col) += nm_xxt_inp(iv_ord_row[i], iv_ord_col[i]);
        }
      }
      if (!IntegerMatrix::is_na(im_rip(iv_ord_row[i], iv_ord_col[i]))) { // Riparian zone
        nm_xxr_inp(iv_ord_row[i], iv_ord_col[i]) = nm_xxt(iv_ord_row[i], iv_ord_col[i]) * (1.0 - nm_tfc_rip(iv_ord_row[i], iv_ord_col[i])); // Retention of riparian zone
        nm_xxt_inp(iv_ord_row[i], iv_ord_col[i]) = nm_xxt(iv_ord_row[i], iv_ord_col[i]) - nm_xxr_inp(iv_ord_row[i], iv_ord_col[i]); // Transport of riparian zone
      }
    } else { // Channel
      nm_xxt_inp_foc = nm_xxt_inp(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving transport of riparian zone and inlet data window
      im_cha_foc     =     im_cha(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving channel data window
      nm_xxr(iv_ord_row[i], iv_ord_col[i]) = focal_sum_matrix_cond(nm_xxt_inp_foc, im_fDi_foc, im_foc) * (1.0 - nm_tfc_lcl(iv_ord_row[i], iv_ord_col[i])) + focal_sum_two_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc, im_cha_foc) * (1.0 - nm_tfc_tpt(iv_ord_row[i], iv_ord_col[i]));
      if (NumericMatrix::is_na(nm_xxt_out(iv_ord_row[i], iv_ord_col[i]))) {
        nm_xxt(iv_ord_row[i], iv_ord_col[i]) = focal_sum_matrix_cond(nm_xxt_inp_foc, im_fDi_foc, im_foc) + focal_sum_two_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc, im_cha_foc) - nm_xxr(iv_ord_row[i], iv_ord_col[i]);
      } else {
        nm_xxt(iv_ord_row[i], iv_ord_col[i]) = focal_sum_matrix_cond(nm_xxt_inp_foc, im_fDi_foc, im_foc) + focal_sum_two_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc, im_cha_foc) + nm_xxt_out(iv_ord_row[i], iv_ord_col[i]) - nm_xxr(iv_ord_row[i], iv_ord_col[i]);
      }
    }
  }

  /* Cell loads & transfers
     ----------------------
  */
  double ns_xxt_cld;
  double ns_xxt_ctf;
  NumericMatrix nm_xxt_ctf_foc;

  NumericMatrix nm_xxe_net = matrix_minus_matrix_elementwise(nm_xxe, nm_xxr); // Net emission

  n = iv_ord_ovl_row_rev.size();
  for (int i = 0; i < n; ++i) { // Bottom-up calculation of cell loads and transfers (overland)
    im_fDi_foc = im_fDi;
    moving_data_window(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i], is_rws, is_cls, is_row_min, is_row_max, is_col_min, is_col_max, im_fDi_foc);
    im_foc     = im_dir(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving potential inflow direction window
    nm_xxt_foc = nm_xxt(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving transport data window
    if (!NumericMatrix::is_na(nm_xxt_inp(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]))) { // Riparian zone or inlet
      ns_xxt_cld = nm_xxt_inp(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]); // Intermediate cell load equals input into channel
    } else { // Other cell
      ns_xxt_cld = nm_xxt_ctf(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]); // Intermediate cell load equals cell transfer from downstream
    }
    ns_xxt_ctf = ns_xxt_cld; // Intermediate cell transfer
    ns_xxt_cld = nm_xxt(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]) * ns_xxt_cld / (nm_xxt(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]) + focal_sum_matrix_cond(nm_xxt_foc, im_fDi_foc, im_foc)); // Calculated intermediate cell load
    if (nm_xxe_net(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]) < 0.0 || !std::isfinite(ns_xxt_cld)) { // Max possible cell load
      ns_xxt_cld = 0.0;
    } else {
      ns_xxt_cld = std::min(ns_xxt_cld, nm_xxe_net(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]));
    }
    nm_xxt_cld(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]) = ns_xxt_cld; // Cell load
    ns_xxt_ctf = std::max(ns_xxt_ctf - ns_xxt_cld, 0.0); // Cell load carry equals cell transfer
    nm_xxt_ctf(iv_ord_ovl_row_rev[i], iv_ord_ovl_col_rev[i]) = ns_xxt_ctf; // Cell transfer
    nm_xxt_ctf_foc = nm_xxt_ctf(Range(is_row_min, is_row_max), Range(is_col_min, is_col_max)); // Moving cell transfer data window
    focal_apportionment(ns_xxt_ctf, nm_xxt_ctf_foc, nm_xxt_foc, im_fDi_foc, im_foc); // Cell load/transfer apportionment
    replace_submatrix(nm_xxt_ctf, nm_xxt_ctf_foc, is_row_min, is_col_min); // Inserts cell load/transfer apportionment into cell transfer
  }

  List result = List::create(
    Named("nm_xxr"    ) = nm_xxr    ,
    Named("nm_xxt"    ) = nm_xxt    ,
    Named("nm_xxt_inp") = nm_xxt_inp,
    Named("nm_xxt_out") = nm_xxt_out,
    Named("nm_xxt_cld") = nm_xxt_cld,
    Named("nm_xxt_ctf") = nm_xxt_ctf
  );

  return result;
}
