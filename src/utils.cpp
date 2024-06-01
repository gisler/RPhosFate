#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// void moving_data_window(int is_ord_row, int is_ord_col, int is_rws, int is_cls, int &is_row_min, int &is_row_max, int &is_col_min, int &is_col_max, IntegerMatrix &im_fDx_foc) {
//   if (is_ord_row == 0) { // Left boundary of moving data window
//     is_row_min = 0;
//     im_fDx_foc = im_fDx_foc(Range(1, im_fDx_foc.nrow() - 1), Range(0, im_fDx_foc.ncol() - 1)); // Adjusts moving inflow direction window
//   } else {
//     is_row_min = is_ord_row - 1;
//   }
//   if (is_ord_row == is_rws - 1) { // Right boundary of moving data window
//     is_row_max = is_rws - 1;
//     im_fDx_foc = im_fDx_foc(Range(0, im_fDx_foc.nrow() - 2), Range(0, im_fDx_foc.ncol() - 1));
//   } else {
//     is_row_max = is_ord_row + 1;
//   }
//   if (is_ord_col == 0) { // Top boundary of moving data window
//     is_col_min = 0;
//     im_fDx_foc = im_fDx_foc(Range(0, im_fDx_foc.nrow() - 1), Range(1, im_fDx_foc.ncol() - 1));
//   } else {
//     is_col_min = is_ord_col - 1;
//   }
//   if (is_ord_col == is_cls - 1) { // Bottom boundary of moving data window
//     is_col_max = is_cls - 1;
//     im_fDx_foc = im_fDx_foc(Range(0, im_fDx_foc.nrow() - 1), Range(0, im_fDx_foc.ncol() - 2));
//   } else {
//     is_col_max = is_ord_col + 1;
//   }
// }
