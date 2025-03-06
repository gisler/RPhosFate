#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "DInfWindow.h"

// #define DEBUGGING

// [[Rcpp::export]]
Rcpp::List transportCpp(
  const arma::dmat& nm_acc_inf,
  const arma::dmat& nm_dir_inf,
  const arma::dmat& nm_slp_cap,
  const arma::dmat& nm_man,
  const arma::dmat& nm_xxe,
  const arma::imat& im_cha,
  const arma::imat& im_rds,
  const arma::imat& im_rip,
  const arma::imat& im_inl,
  const Rcpp::String& substance,
  const Rcpp::S4& parameters,
  const Rcpp::S4& helpers,
  const int is_ths = 1
) {
  DinfWindow dinfWindow {nm_dir_inf.n_rows, nm_dir_inf.n_cols};

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
      double ns_dir_inf {nm_dir_inf.at(i, j)};

      if (Rcpp::NumericMatrix::is_na(nm_acc_inf.at(i, j)) ||
          Rcpp::NumericMatrix::is_na(ns_dir_inf) ||
          ns_dir_inf == -1.0) { // A value of -1.0 is used by the WhiteboxTool DInfPointer to designate cells with no flow direction
        continue;
      }

      arma::dvec8 nv_ifl_p {dinfWindow.get_ifl_p(nm_dir_inf, i, j)};
      im_ifl.at(i, j) = arma::accu(
        dinfWindow.get_ifl_xp<double>(nv_ifl_p, i, j, nm_acc_inf) > 0.0
      ); // The flow accumulation can be NA_REAL where the flow direction and therefore the inflow proportion is not (e.g. roads)
    }
  }

  /* Determine order of row and column indices
   * -----------------------------------------
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

  #ifdef DEBUGGING
  Rcpp::Rcout << "Debugging" << std::endl;
  arma::imat im_ord(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_INTEGER)
  );
  #endif

  for (arma::uword n = 0; n < ord.uv_r.capacity(); ++n) {
    if (n == ord.uv_r.size()) {
      Rcpp::stop("Warning: Could not determine a hydrologic consistent transport calculation order.");
    }

    #ifdef DEBUGGING
    im_ord.at(ord.uv_r[n], ord.uv_c[n]) = n;
    #endif

    FacetProperties fct {dinfWindow.get_ofl_facetProperties(
      nm_dir_inf.at(ord.uv_r[n], ord.uv_c[n]),
      ord.uv_r[n],
      ord.uv_c[n]
    )};
    arma::sword is_x1 {NA_INTEGER}, is_x2 {NA_INTEGER};

    if (!fct.ls_x1_oob) {
      is_x1 = im_ifl.at(fct.us_x1_r, fct.us_x1_c);
      if (!Rcpp::IntegerMatrix::is_na(is_x1)) {
        --is_x1;
        im_ifl.at(fct.us_x1_r, fct.us_x1_c) = is_x1;
      }
    }
    if (!fct.ls_x2_oob) {
      is_x2 = im_ifl.at(fct.us_x2_r, fct.us_x2_c);
      if (!Rcpp::IntegerMatrix::is_na(is_x2)) {
        --is_x2;
        im_ifl.at(fct.us_x2_r, fct.us_x2_c) = is_x2;
      }
    }

    if (is_x1 == 0) {
      ord.uv_r.push_back(fct.us_x1_r);
      ord.uv_c.push_back(fct.us_x1_c);
    }
    if (is_x2 == 0) {
      ord.uv_r.push_back(fct.us_x2_r);
      ord.uv_c.push_back(fct.us_x2_c);
    }
  }

  /* Retentions and transports
   * =========================
   */
  const double ns_rhy_a {parameters.slot("ns_rhy_a")};
  const double ns_rhy_b {parameters.slot("ns_rhy_b")};
  const double ns_cha_rto {parameters.slot("ns_cha_rto")};
  const double ns_man_rip {parameters.slot("ns_man_rip")};
  const double ns_man_cha {parameters.slot("ns_man_cha")};

  const Rcpp::NumericVector nv_tfc_inl = parameters.slot("nv_tfc_inl");
  const double ns_tfc_inl {nv_tfc_inl[substance]};
  const Rcpp::NumericVector nv_enr_rto = parameters.slot("nv_enr_rto");
  const double ns_enr_rto {(substance == "SS") ? 1.0 : nv_enr_rto[substance]};
  const double ns_dep_ovl_tmp {parameters.slot("ns_dep_ovl")};
  const double ns_dep_ovl {ns_dep_ovl_tmp / ns_enr_rto};
  const double ns_dep_cha {parameters.slot("ns_dep_cha")};

  const double ns_res {helpers.slot("ns_res")};
  const double ns_siz {helpers.slot("ns_siz")};

  // Flow path length of riparian zone
  const double ns_fpl_rip {(ns_res - ns_res * ns_cha_rto) / 2.0};
  // Strickler coefficient of riparian zone
  const double ns_str_rip {1.0 / ns_man_rip};

  arma::dmat nm_xxr(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_REAL)
  );
  arma::dmat nm_xxt(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_REAL)
  );
  arma::dmat nm_xxt_inp(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_REAL)
  );
  arma::dmat nm_xxt_rip(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_REAL)
  ); // Substance inputs of riparian zone into surface waters
  arma::dmat nm_xxt_out(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_REAL)
  );

  /* Calculation of retentions and transports
   * ----------------------------------------
   */
  for (arma::uword n = 0; n < ord.uv_r.size(); ++n) {
    arma::uword i {ord.uv_r[n]};
    arma::uword j {ord.uv_c[n]};

    int is_cha {im_cha.at(i, j)};
    int is_rip {im_rip.at(i, j)};
    int is_inl {im_inl.at(i, j)};
    double ns_xxe {nm_xxe.at(i, j)};
    double ns_dir_inf {nm_dir_inf.at(i, j)};
    // Capped slope in m / m
    double ns_slp_cap {nm_slp_cap.at(i, j) / 100.0};

    // Flow direction in radian
    double ns_dir_inf_rad {ns_dir_inf * arma::datum::pi / 180.0};
    // Flow path length
    double ns_fpl {ns_res * (std::abs(std::sin(ns_dir_inf_rad)) +
      std::abs(std::cos(ns_dir_inf_rad)))};
    // Strickler coefficient
    double ns_str {(Rcpp::IntegerMatrix::is_na(is_cha)) ?
      1.0 / nm_man.at(i, j) : 1.0 / ns_man_cha};
    // Hydraulic radius
    double ns_rhy {ns_rhy_a * std::pow(
      nm_acc_inf.at(i, j) * ns_siz * 1e-6,
      ns_rhy_b
    )};
    // Hydraulic radius and slope part of flow velocity
    double ns_rhy_slp {std::pow(ns_rhy, 2.0 / 3.0) * std::sqrt(ns_slp_cap)};
    // Residence time
    double ns_rtm {ns_fpl / (ns_str * ns_rhy_slp)};

    // Inflow proportions
    arma::dvec8 nv_ifl_p {dinfWindow.get_ifl_p(nm_dir_inf, i, j)};

    // Overland cell
    if (Rcpp::IntegerMatrix::is_na(is_cha)) {
      // Retention coefficients
      double ns_rtc_ifl {1.0 - std::exp(-ns_dep_ovl * ns_rtm)};
      double ns_rtc_lcl {1.0 - std::exp(-ns_dep_ovl * ns_rtm * 0.5)};

      // Inflowing load
      double ns_xxt_ifl {arma::accu(
        dinfWindow.get_ifl_xp<double>(nv_ifl_p, i, j, nm_xxt)
      )};

      // Retention
      double ns_xxr {ns_xxt_ifl * ns_rtc_ifl + ns_xxe * ns_rtc_lcl};
      nm_xxr.at(i, j) = ns_xxr;
      // Transport
      double ns_xxt {ns_xxt_ifl + ns_xxe - ns_xxr};
      nm_xxt.at(i, j) = ns_xxt;

      // Riparian zone or inlet cell
      if (!Rcpp::IntegerMatrix::is_na(is_rip) ||
          !Rcpp::IntegerMatrix::is_na(is_inl)) {
        // Outflowing instead of inflowing point of view
        FacetProperties fct {dinfWindow.get_ofl_facetProperties(ns_dir_inf, i, j)};

        // Riparian zone cell
        if (!Rcpp::IntegerMatrix::is_na(is_rip)) {
          X1X2<int> cha1cha2 = dinfWindow.get_ofl_x1x2<int>(fct, im_cha, NA_INTEGER);

          // Retention coefficient (0.0 in case there is no riparian zone defined)
          double ns_rtc_rip {};
          if (ns_cha_rto < 1.0) {
            // Residence time
            double ns_rtm_rip {ns_fpl_rip / (ns_str_rip * ns_rhy_slp)};

            ns_rtc_rip = 1.0 - std::exp(-ns_dep_ovl * ns_rtm_rip);
          } else {
            ns_rtc_rip = 0.0;
          }

          // Retention, transport and substance input (of riparian zone) into
          // surface water
          nm_xxt_inp.at(i, j) = dinfWindow.inc_ofl_x1x2<int, INTSXP>(
            cha1cha2,
            ns_xxt - ns_xxt * ns_rtc_rip,
            nm_xxt_rip
          );
        }

        // Inlet cell
        if (!Rcpp::IntegerMatrix::is_na(is_inl)) {
          X1X2<int> rds1rds2 = dinfWindow.get_ofl_x1x2<int>(fct, im_rds, NA_INTEGER);

          // Proportional transport
          double ns_xxt_x1x2 {0.0};
          if (!Rcpp::IntegerMatrix::is_na(rds1rds2.x1)) {
            ns_xxt_x1x2 += ns_xxt * rds1rds2.fct.ns_p1;
          }
          if (!Rcpp::IntegerMatrix::is_na(rds1rds2.x2)) {
            ns_xxt_x1x2 += ns_xxt * rds1rds2.fct.ns_p2;
          }

          // Retention, transport and (additional) substance input into surface
          // water
          ns_xxt_x1x2 *= ns_tfc_inl;
          double ns_xxt_inp {nm_xxt_inp.at(i, j)};
          if (Rcpp::NumericMatrix::is_na(ns_xxt_inp)) {
            ns_xxt_inp = 0.0;
          }
          ns_xxt_inp += ns_xxt_x1x2; // A cell can be a riparian zone and inlet at the same time
          nm_xxt_inp.at(i, j) = ns_xxt_inp;

          // Outlet row and col from inlet code (C++ indices start at 0)
          std::div_t code {std::div(im_inl.at(i, j), dinfWindow.is_cls)};
          arma::uword us_row {static_cast<arma::uword>(code.quot - 1)};
          arma::uword us_col {static_cast<arma::uword>(code.rem  - 1)};

          // Outlet load
          double ns_xxt_out {nm_xxt_out.at(us_row, us_col)};
          if (Rcpp::NumericMatrix::is_na(ns_xxt_out)) {
            ns_xxt_out = 0.0;
          }
          nm_xxt_out(us_row, us_col) = ns_xxt_out + ns_xxt_inp;
        }
      }

    // Channel cell
    } else {
      // Retention coefficients
      double ns_rtc_ifl {1.0 - std::exp(-ns_dep_cha * ns_rtm)};
      double ns_rtc_lcl {1.0 - std::exp(-ns_dep_cha * ns_rtm * 0.5)};

      // Inflowing overland load
      double ns_xxt_rip {nm_xxt_rip.at(i, j)};
      if (Rcpp::NumericMatrix::is_na(ns_xxt_rip)) {
        ns_xxt_rip = 0.0;
      }
      // Inflowing channel load
      arma::dvec8 nv_xxt_cha {
        dinfWindow.get_ifl_xp<double>(nv_ifl_p, i, j, nm_xxt) %
          (dinfWindow.get_ifl_xp<int>(nv_ifl_p, i, j, im_cha) > 0.0)
      };
      double ns_xxt_cha {arma::accu(nv_xxt_cha)};
      // Outlet load
      double ns_xxt_out {nm_xxt_out.at(i, j)};
      if (Rcpp::NumericMatrix::is_na(ns_xxt_out)) {
        ns_xxt_out = 0.0;
      }

      // Retention
      double ns_xxr {ns_xxt_cha * ns_rtc_ifl +
        (ns_xxt_rip + ns_xxt_out) * ns_rtc_lcl};
      nm_xxr.at(i, j) = ns_xxr;
      // Transport
      nm_xxt.at(i, j) = ns_xxt_cha + ns_xxt_rip + ns_xxt_out - ns_xxr;
    }
  }

  /* Cell loads and transfers
   * ========================
   */
  arma::dmat nm_xxt_ctf(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_REAL)
  );
  arma::dmat nm_xxt_cld(
    arma::size(nm_dir_inf),
    arma::fill::value(NA_REAL)
  );

  for (arma::uword n = ord.uv_r.size() - 1; n-- > 0; ) {
    arma::uword i {ord.uv_r[n]};
    arma::uword j {ord.uv_c[n]};

    if (!Rcpp::IntegerMatrix::is_na(im_cha.at(i, j))) {
      continue;
    }

    // Inflow proportions
    arma::dvec8 nv_ifl_p {dinfWindow.get_ifl_p(nm_dir_inf, i, j)};

    // Inflowing overland load
    arma::dvec8 nv_xxt_ifl {
      dinfWindow.get_ifl_xp<double>(nv_ifl_p, i, j, nm_xxt)
    };
    double ns_xxt_ifl {arma::accu(nv_xxt_ifl)};
    // Net emission
    double ns_xxe_net {nm_xxe.at(i, j) - nm_xxr.at(i, j)};

    // Initialise intermediate cell load (step 1)
    double ns_xxt_cld {nm_xxt_ctf.at(i, j)};
    if (Rcpp::NumericMatrix::is_na(ns_xxt_cld)) {
      ns_xxt_cld = 0.0;
    }
    double ns_xxt_inp {nm_xxt_inp.at(i, j)};
    if (!Rcpp::NumericMatrix::is_na(ns_xxt_inp)) {
      ns_xxt_cld += ns_xxt_inp;
    }
    // Initialise intermediate cell transfer (step 2)
    double ns_xxt_ctf {ns_xxt_cld};

    // Cell load (steps 3 and 4)
    if (ns_xxe_net <= 0.0) {
      ns_xxt_cld = 0.0;
    } else {
      double ns_xxt {nm_xxt.at(i, j)};
      ns_xxt_cld = ns_xxt_cld * ns_xxt / (ns_xxt + ns_xxt_ifl);

      ns_xxt_cld = std::min(ns_xxt_cld, ns_xxe_net);
    }
    nm_xxt_cld.at(i, j) = ns_xxt_cld;
    // Cell transfer (steps 5 and 6)
    ns_xxt_ctf = std::max(ns_xxt_ctf - ns_xxt_cld, 0.0);
    nm_xxt_ctf.at(i, j) = ns_xxt_ctf;

    // Weighted apportioning (step 7)
    for (arma::uword k = 0; k < nv_ifl_p.n_elem; ++k) {
      arma::uword us_row {i + ifl.iv_dr[k]};
      arma::uword us_col {j + ifl.iv_dc[k]};

      if (nv_ifl_p[k] > 0.0 &&
          !Rcpp::NumericMatrix::is_na(nm_acc_inf.at(us_row, us_col))) {
        double ns_xxt_ctf_old {nm_xxt_ctf.at(us_row, us_col)};
        if (Rcpp::NumericMatrix::is_na(ns_xxt_ctf_old)) {
          ns_xxt_ctf_old = 0.0;
        }

        nm_xxt_ctf.at(us_row, us_col) = ns_xxt_ctf_old +
          ns_xxt_ctf * nv_xxt_ifl[k] / ns_xxt_ifl;
      }
    }
  }

  return Rcpp::List::create(
    #ifdef DEBUGGING
    Rcpp::Named("im_ifl"    ) = im_ifl    ,
    Rcpp::Named("im_ord"    ) = im_ord    ,
    Rcpp::Named("nm_xxt_rip") = nm_xxt_rip,
    #endif
    Rcpp::Named("nm_xxr"    ) = nm_xxr    ,
    Rcpp::Named("nm_xxt"    ) = nm_xxt    ,
    Rcpp::Named("nm_xxt_inp") = nm_xxt_inp,
    Rcpp::Named("nm_xxt_out") = nm_xxt_out,
    Rcpp::Named("nm_xxt_ctf") = nm_xxt_ctf,
    Rcpp::Named("nm_xxt_cld") = nm_xxt_cld
  );
}
