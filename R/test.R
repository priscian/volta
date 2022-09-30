#' @export
test_volta <- function()
{
  ### Maecker & Trotter 2006

  e <- new.env()
  load(system.file("extdata/peak2-maecker+trotter2006_2019-09-19+54111.RData", package = "volta"), envir = e)
  attr(e$d, "wavelength") <- c(520, 785, 578, 695) # This attribute is necessary!
  attr(e$d, "experiment_name") <- "Maecker & Trotter 2006" # Not necessary to set this attribute.
  v2006 <-
    plot_peak2_data(
      x = list(e$d),
      image_dir = system.file("images", package = "volta"),
      plot_series... = list(x_var = "voltage", lwd = 4, xaxisTicks... = list(nint = 20, force = TRUE), cex.xaxis = 0.8),
      create_smooth_variables... = list(x_var = "voltage"),
      points... = list(cex = 3, lwd = 3, type = "p"),
      save_png = FALSE
    )

  v2006[[1]]$table
}
## Save the resulting plot:
# keystone::dev_print(file = "./Downloads/001 - Maecker & Trotter 2006.png", device = png, width = 9.375, height = 7.3, units = "in", res = 600)
