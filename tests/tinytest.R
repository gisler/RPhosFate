if (requireNamespace("tinytest", quietly = TRUE)) {
  tinytest::test_package(
    "RPhosFate",
    at_home = isTRUE(as.logical(Sys.getenv("RPhosFate_DEVELOPER")))
  )
}
