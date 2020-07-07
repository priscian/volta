#' @export
test_volta <- function()
{
  ### Maecker & Trotter 2006

  e <- new.env()
  load(system.file("extdata/peak2-maecker+trotter2006_2019-09-19+54111.RData", package = "volta"), envir = e)
  attr(e$d, "wavelength") <- c(520, 785, 578, 695) # This attribute is necessary!
  attr(e$d, "experiment_name") <- "Maecker & Trotter 2006" # Not necessary to set this attribute.
  v2006 <- plot_peak2_data(x = list(e$d), image_dir = system.file("images", package = "volta"), save_png = FALSE)

  print(v2006)

  nop()
}
