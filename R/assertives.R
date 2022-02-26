assertdf_cdt <- function(cmt) {
  assertDataFrame(
    cmt@parameters@df_cdt,
    min.rows = 1L,
    min.cols = 4L,
    .var.name = "df_cdt"
  )
  assertCharacter(
    names(cmt@parameters@df_cdt),
    min.chars = 1L,
    any.missing = FALSE,
    min.len = 4L,
    unique = TRUE,
    .var.name = "names(df_cdt)"
  )
  assertSubset(
    c("ID", "x", "y"),
    names(cmt@parameters@df_cdt),
    .var.name = "names(df_cdt)"
  )
  assertVector(
    cmt@parameters@df_cdt[["ID"]],
    strict = TRUE,
    any.missing = FALSE,
    min.len = 1L,
    unique = TRUE,
    .var.name = 'df_cdt[["ID"]]'
  )
  assertNumeric(
    cmt@parameters@df_cdt[["x"]],
    lower = cmt@helpers@ex_cmt@xmin,
    upper = cmt@helpers@ex_cmt@xmax,
    any.missing = FALSE,
    min.len = 1L,
    .var.name = 'df_cdt[["x"]]'
  )
  assertNumeric(
    cmt@parameters@df_cdt[["y"]],
    lower = cmt@helpers@ex_cmt@ymin,
    upper = cmt@helpers@ex_cmt@ymax,
    any.missing = FALSE,
    min.len = 1L,
    .var.name = 'df_cdt[["y"]]'
  )
}

assertCol <- function(cmt, col) {
  assertdf_cdt(cmt)
  assertChoice(col, setdiff(names(cmt@parameters@df_cdt), c("ID", "x", "y")))
  assertNumeric(
    cmt@parameters@df_cdt[[col]],
    lower = 0,
    finite = TRUE,
    all.missing = FALSE,
    min.len = 1L,
    .var.name = "df_cdt[[col]]"
  )
}
