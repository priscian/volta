#' @export
test_volta <- function()
{
  #cordon(get_peak2_data, path = system.file("extdata", package = "volta"), envir = globalenv(), file_path = paste(dataDir, "peak2-historical.RData", sep = "/"), variables = c("p2", "l"), timestamp... = list(use_seconds = TRUE), action = "run")
  l <- get_peak2_spreadsheet(path = system.file("extdata/Peak2 Overview 122909.xlsx", package = "volta"))

  ## Plot historical Peak 2 data.
  #cordon(plot_peak2_data, x = l, image_dir = "report/report_" %_% make_current_timestamp("%Y%m%d"), save_png = FALSE, envir = globalenv(), file_path = paste(dataDir, "peak2-historical-changepoints.RData", sep = "/"), variables = NULL, timestamp... = list(use_seconds = TRUE), action = "run")
  v <- plot_peak2_data(x = l)

  print(v)

  ### Maecker & Trotter 2006

  e <- new.env()
  load(system.file("extdata/peak2-maecker+trotter2006_2019-09-19+54111.RData", package = "volta"), envir = e)
  attr(e$d, "wavelength") <- c(520, 785, 578, 695) # This attribute is necessary!
  attr(e$d, "filename") <- "Maecker & Trotter 2006" # Not necessary to set this attribute.
  v2006 <- plot_peak2_data(x = list(e$d), image_dir = system.file("images", package = "volta"), save_png = FALSE)

  print(v2006)

  nop()
}
