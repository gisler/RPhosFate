writeParameters <- function(parameters) {
  parameters <- slots2list(parameters)

  parameters[["nv_enr_rto"]] <- as.list(parameters[["nv_enr_rto"]])

  write_yaml(
    c(list(RPhosFate = as.character(packageVersion("RPhosFate"))), parameters),,
    "parameters.yaml",
    indent.mapping.sequence = TRUE
  )
}

readParameters <- function(arguments) {
  parameters <- read_yaml("parameters.yaml")

  parameters[["nv_enr_rto"]] <- unlist(parameters[["nv_enr_rto"]])
  parameters[["nm_olc"]] <- matrix(parameters[["nm_olc"]], 1L, 2L)
  parameters[["df_cdt"]] <- as.data.frame(parameters[["df_cdt"]])

  modifyList(parameters, arguments)
}
