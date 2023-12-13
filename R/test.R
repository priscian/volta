#' @export
test_volta <- function()
{
  ### Maecker & Trotter 2006

  e <- new.env()
  load(system.file("extdata/peak2-maecker+trotter2006_2019-09-19+54111.RData", package = "volta"),
    envir = e)
  p2 <- list(`Maecker & Trotter 2006` = structure(e$d %>% dplyr::rename(quantity = 1),
    color = keystone::wavelength2col(c(520, 785, 578, 695)), # This attribute is necessary!
    ## N.B. Not necessary to set this attribute:
    experiment_name = "Maecker & Trotter 2006"))
  fip <- p2 %>% find_inflection_point()
  p2 <- sapply(names(p2),
    function(i) { (structure(p2[[i]], plot_data = fip[[i]]) %>% keystone::add_class("volta")) },
    simplify = FALSE)

  plot(p2[[1]])

  attr(p2[[1]], "plot_data")$inflection_points %>%
    dplyr::rename(
      !!names(formals(plot.volta)$x_var_lab %>% keystone::poly_eval()) := "inflection",
      !!names(formals(plot.volta)$y_var_lab %>% keystone::poly_eval()) := "y"
    )
}
## Save the resulting plot:
# keystone::dev_print(file = "./Downloads/001 - Maecker & Trotter 2006.png",
#   device = png, width = 9.375, height = 7.3, units = "in", res = 600)


#' @export
test_volta_2 <- function()
{
}
