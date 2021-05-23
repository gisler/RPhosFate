assertCol <- function(cmt, col) {
  qassert(col, "S1")
  assertSubset(col, setdiff(names(cmt@parameters@df_cdt), c("ID", "x", "y")))
}

assertSubstance <- function(cmt, substance) {
  qassert(substance, "S1")
  assertSubset(substance, slotNames(cmt@substance))
}
