assertdf_cdt <- function(cmt) {
  assertDataFrame(cmt@parameters@df_cdt, min.rows = 1L, min.cols = 4L)
  assertCharacter(
    names(cmt@parameters@df_cdt),
    min.chars = 1L,
    any.missing = FALSE,
    min.len = 4L,
    unique = TRUE
  )
  assertSubset(c("ID", "x", "y"), names(cmt@parameters@df_cdt))
  assertVector(
    cmt@parameters@df_cdt[["ID"]],
    strict = TRUE,
    any.missing = FALSE,
    min.len = 1L,
    unique = TRUE
  )
  qassert(cmt@parameters@df_cdt[["x"]], "N+(0,)")
  qassert(cmt@parameters@df_cdt[["y"]], "N+(0,)")
}

assertCol <- function(cmt, col) {
  qassert(col, "S1")
  assertdf_cdt(cmt)
  assertSubset(col, setdiff(names(cmt@parameters@df_cdt), c("ID", "x", "y")))
  assertNumeric(
    cmt@parameters@df_cdt[[col]],
    lower = 0,
    finite = TRUE,
    all.missing = FALSE,
    min.len = 1L,
  )
}

assertSubstance <- function(cmt, substance) {
  qassert(substance, "S1")
  assertSubset(substance, slotNames(cmt@substance))
}
