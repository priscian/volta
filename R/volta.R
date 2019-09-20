#' @export
get_peak2_data <- function(path, list.files... = list())
{
  list.filesArgs <- list(
    path = path,
    pattern = "\\.xlsx",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  list.filesArgs <- utils::modifyList(list.filesArgs, list.files...)

  if (R.utils::isDirectory(path[1]))
    ff <- do.call(list.files, list.filesArgs)
  else
    ff <- path

  r <- sapply(ff,
    function(a)
    {
      r <- NULL

      nmRe <- "\b?([2-9][0-9][0-9])\b?"
      ccRe <- "(red|orange|green|cyan|blue|indigo|violet)"

      tryCatch({
        numSheets <- length(readxl::excel_sheets(a))
      }, error = function(e) { cat("Error: ", e$message, "\n", sep = ""); cat(basename(a) %_% " cannot be read as spreadsheet.\n") })

      if (!exists("numSheets", where = environment()))
        return (NULL)

      y <- sapply(seq(numSheets),
        function(b)
        {
          Error <- function(e) {
            cat("Error: ", e$message, "\n", sep = "")
            cat(basename(a) %_% ", sheet " %_% b %_% " cannot be imported correctly.\n")
          }

          tryCatch({
            x <- rio::import(a, which = b, n_max = 16, na = c("", "##ERROR##"))
            x <- x %>% janitor::remove_empty(which = "cols") %>% dplyr::filter_all(dplyr::any_vars(!is.na(.)))

            ## Figure out which column is voltage & which other columns to keep.
            for(i in seq_along(x)) {
              testVoltages <- stringr::str_match(x[[i]], nmRe)[, 2] %>% as.numeric
              is.na(testVoltages) <- testVoltages < 200 | testVoltages >= 1000
              if (all(Vectorize(is_invalid)(testVoltages)))
                next

              voltageColumn <- i

              break
            }

            voltageColumnIndex <- grepl("^voltage", colnames(x), ignore.case = TRUE, perl = TRUE)
            if (!is_invalid(which(voltageColumnIndex)))
              voltageColumn <- which(voltageColumnIndex)[1]

            y <- x[, voltageColumn:NCOL(x)]

            ## Find channel columns.
            cci <- grep(ccRe, colnames(y), ignore.case = TRUE, perl = TRUE, value = FALSE)
            if (is_invalid(colnames(y)[cci]))
              stop("No valid channel columns.")

            y[[1]] <- stringr::str_match(y[[1]], nmRe)[, 2] %>% as.numeric
            dev_null <- sapply(2:NCOL(y), function(z) { y[[z]] <<- as.numeric(y[[z]]) })
            y <- y[!is.na(y[[1]]), c(1, cci)]
            y <- y %>% janitor::remove_empty(which = "cols") %>%
              dplyr::filter_all(dplyr::any_vars(!is.na(.))) %>%
              naniar::replace_with_na_all(condition = ~.x <= 0.0)
            names(y)[1] <- "voltage"

            color <- tolower(stringr::str_match(colnames(y)[2:NCOL(y)], stringr::regex(ccRe, ignore_case = TRUE))[, 2])
            wavelength <- stringr::str_match(colnames(y)[2:NCOL(y)], nmRe)[, 2]
            if (all(Vectorize(is_invalid)(wavelength)))
              wavelength <- color_nm_map[color]
            wavelength <- as.numeric(wavelength)

            if (NROW(y) == 0) stop("No rows in data set.")

            # print(y$voltage)
            # print(voltageColumn)
            # print(color)
            # print(wavelength)

            #if (any(y$voltage < 50 | y$voltage > 999)) { print(a); browser() }

            attr(y, "filename") <- a
            attr(y, "sheet") <- b
            attr(y, "color") <- color
            attr(y, "wavelength") <- wavelength

            y
          }, error = Error) # Could add 'warning = Error'.
        }, simplify = FALSE)

      r <- y

      r
    }, simplify = FALSE)

  p2 <- r

  l <- p2 %>% purrr::flatten() %>% # Flatten list to a depth of 1
    purrr::compact() # Remove NULL elements from list

  l
}


#' @export
plot_peak2_data <- function(x, image_dir, save_png = FALSE)
{
  if (save_png && !dir.exists(image_dir))
    dir.create(image_dir, recursive = TRUE)

  pngArgs <- list(
    width = 12.5,
    height = 7.3,
    units = "in",
    res = 600
  )

  indices <- NULL # NULL for all elements
  if (is.null(indices))
    indices <- seq_along(x)
  currentIndex <- head(indices, 1)

  r <- sapply(x[indices],
    function(a)
    {
      color <- wavelength2col(attr(a, "wavelength"))
      sheet <- attr(a, "sheet"); if (is_invalid(sheet)) sheet <- 1
      sheetText <- "(sheet " %_% sheet %_% ")"
      filename <- attr(a, "filename"); if (is_invalid(filename)) filename <- "[no filename]"

      cat("Processing object # ", currentIndex, ", ", basename(filename), " (sheet ", sheet, ")\n", sep = ""); flush.console()

      x <- a %>% dplyr::mutate_at(dplyr::vars(2:NCOL(.)), log10)
      series <- names(x)[2:NCOL(x)]

      Error <- function(e) {
        cat("Fatal error: ", e$message, ". Aborting & continuing next analysis.\n", sep = ""); flush.console()
        if (dev.cur() >= 0L) dev.off()
      }

      ## Use time-series plotting.
      tryCatch({
        if (save_png) do.call(grDevices::png, utils::modifyList(pngArgs, list(filename = paste(image_dir, sprintf("%03d", currentIndex) %_% " - " %_% paste(basename(filename), sheetText) %_% ".png", sep = "/"))))

        r <- plot_series(x, series, x_var = "voltage", log = "", xlab = "PMT Voltage", ylab = expression(paste(log[10], " CV")), main = paste(basename(filename), sheetText), col = color, lwd = 1, trend = FALSE, segmented = TRUE, segmented... = list(breakpoints... = list(h = 3)))

        volts <- sapply(r$segmented$piecewise, function(a) tail(a$sm$psi[, "Est."], 1), simplify = FALSE)
        cv <- sapply(r$segmented$piecewise, function(a) { if (is.null(a$sm)) return (NULL); predict(a$sm, dataframe(voltage = tail(a$sm$psi[, "Est."], 1))) }, simplify = FALSE)

        ## Fit smooth splines & return functions describing it.
        fxs_spline <- sapply(series, function(b) splinefun(x[["voltage"]], x[[b]]))

        changepoints_cv <- mapply(fxs_spline, cv,
          FUN = function(a, b)
          {
            cp <- stats::approx(x = a(x[["voltage"]]), y = x[["voltage"]], xout = b)

            list(x = cp$y, y = cp$x)
          }, SIMPLIFY = FALSE)

        changepoints_volts <- mapply(fxs_spline, volts,
          FUN = function(a, b)
          {
            yval <- a(b)

            list(x = b, y = yval)
          }, SIMPLIFY = FALSE)

        dev_null <- sapply(changepoints_cv, function(a) points(a, col = "black", pch = 4, cex = 1))
        #dev_null <- sapply(changepoints_volts, function(a) points(a, col = "gray50", pch = 3, cex = 1))

        if (save_png) dev.off()
      }, error = Error) # Could add 'warning = Error'.

      rv <- purrr::map_dfr(changepoints_cv, dataframe); rownames(rv) <- names(changepoints_cv); colnames(rv) <- c("PMT_voltage", "log10_CV")

      currentIndex <<- currentIndex + 1

      rv
    }, simplify = FALSE)

  r
}
