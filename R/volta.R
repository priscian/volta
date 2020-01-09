#' @importFrom magrittr %>%
#' @importFrom plinth %_% cordon is_invalid dataframe

#' @export
get_peak2_spreadsheet <- function(path, list.files... = list())
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
      ccRe <- "(red|orange|green|cyan|blue|indigo|violet|yg|uv)"

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

            attr(y, "experiment_name") <- a
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
get_peak2_data <- function(
  path = ".",
  list.dirs... = list(),
  list.files... = list(),
  read.FCS... = list(),
  min_zero = FALSE,
  cv_fun = plinth::bd_rcv,
  cv_multiplier = 100,
  replace_with_na_condition = ~ is.infinite(.x) | is.nan(.x) | .x < 0.0
)
{
  list.dirsArgs <- list(
    path = path
  )
  list.dirsArgs <- utils::modifyList(list.dirsArgs, list.dirs...)
  dirPaths <- do.call(list.dirs, list.dirsArgs)

  list.filesArgs <- list(
    pattern = ".*?\\.fcs$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )

  ff <- sapply(dirPaths,
    function(i)
    {
      list.filesArgs$path = i
      list.filesArgs <- utils::modifyList(list.filesArgs, list.files...)
      filePaths <- do.call(list.files, list.filesArgs)

      filePaths
    }, simplify = FALSE) %>% purrr::compact() # Remove blank elements

  nmRe <- "\b?([2-9][0-9][0-9])\b?"
  ccRe <- "(red|orange|green|cyan|blue|indigo|violet|yg|uv)"

  p2 <- sapply(ff,
    function(i)
    {
      r <- sapply(i,
        function(j)
        {
          read.FCSArgs <- list(
            filename = j
          )
          read.FCSArgs <- utils::modifyList(read.FCSArgs, read.FCS...)
          fcs <- do.call(flowCore::read.FCS, read.FCSArgs)

          allChannelNames <- flowCore::description(fcs)[grep("^\\$P\\d+N", trimws(names(flowCore::description(fcs))), value = TRUE, perl = TRUE)]
          descriptionNames <- stringr::str_match(names(allChannelNames)[unlist(allChannelNames) %in% trimws(colnames(flowCore::description(fcs)$SPILL))], "^(.*?)N$")[, 2]

          channelVoltages <- flowCore::description(fcs)[descriptionNames %_% "V"]
          channelNames <- flowCore::description(fcs)[descriptionNames %_% "N"]

          e <- flowCore::exprs(fcs)[, unlist(channelNames)]
          if (min_zero)
            e <- plyr::aaply(e, 1, function(k) k - min(k, na.rm = TRUE))

          ee <- e %>% plyr::aaply(2, function(k) { cv_fun(k) * cv_multiplier }) %>% t %>% plinth::dataframe()
          attr(ee, "channel_voltages") <- unlist(channelVoltages)
          attr(ee, "experiment_name") <- flowCore::description(fcs)$`EXPERIMENT NAME`

          ee
        }, simplify = FALSE)

      ## It probably won't happen, but this can separate channels that aren't on the same voltage for some reason.
      rr <- sapply(r,
        function(j)
        {
          channel_voltages <- attr(j, "channel_voltages")
          v <- Reduce(plyr::rbind.fill, sapply(unique(channel_voltages),
            function(k)
            {
              cbind(voltage = k, j[, k %in% channel_voltages])
            }, simplify = FALSE)) %>% dplyr::mutate(voltage = as.numeric(plinth::unfactor(voltage)))

          v
        }, simplify = FALSE)

      rrr <- Reduce(plyr::rbind.fill, rr) %>%
        naniar::replace_with_na_all(condition = replace_with_na_condition)

      ## Add some attributes.
      color <- tolower(stringr::str_match(colnames(rrr)[2:NCOL(rrr)], stringr::regex(ccRe, ignore_case = TRUE))[, 2])
      wavelength <- stringr::str_match(colnames(rrr)[2:NCOL(rrr)], nmRe)[, 2]
      if (all(Vectorize(is_invalid)(wavelength)))
        wavelength <- color_nm_map[color]
      wavelength <- as.numeric(wavelength)

      attr(rrr, "experiment_name") <- attr(r[[1]], "experiment_name")
      attr(rrr, "sheet") <- 1
      attr(rrr, "color") <- color
      attr(rrr, "wavelength") <- wavelength

      rrr
    }, simplify = FALSE)

  l <- p2

  l
}


#' @export
plot_peak2_data <- function(
  x,
  report_dir,
  image_dir = report_dir,
  trans_fun = log10, # Also possibly 'identity'
  replace_with_na_condition = ~ is.infinite(.x) | is.nan(.x),
  save_png = FALSE,
  segmented = FALSE,
  zeros_threshold = 0.95,
  plot_derivative = TRUE
)
{
  if (save_png && !dir.exists(image_dir))
    dir.create(image_dir, recursive = TRUE)

  if (!missing(report_dir) && !dir.exists(report_dir))
    dir.create(report_dir, recursive = TRUE)

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
      color <- plinth::wavelength2col(attr(a, "wavelength"))
      sheet <- attr(a, "sheet"); if (is_invalid(sheet)) sheet <- 1
      sheetText <- "(sheet " %_% sheet %_% ")"
      experiment_name <- attr(a, "experiment_name"); if (is_invalid(experiment_name)) experiment_name <- "[no experiment_name]"

      cat("Processing object # ", currentIndex, ", ", basename(experiment_name), " (sheet ", sheet, ")\n", sep = ""); flush.console()

      x <- a %>% dplyr::mutate_at(dplyr::vars(2:NCOL(.)), trans_fun) %>%
        naniar::replace_with_na_all(condition = replace_with_na_condition)
      series <- names(x)[2:NCOL(x)]

      Error <- function(e) {
        cat("Fatal error: ", e$message, ". Aborting & continuing next analysis.\n", sep = ""); flush.console()
        if (dev.cur() >= 0L) dev.off()
      }

      ### Do time-series plotting.

      #tryCatch({
        if (save_png) do.call(grDevices::png, utils::modifyList(pngArgs, list(filename = paste(image_dir, sprintf("%03d", currentIndex) %_% " - " %_% paste(basename(experiment_name), sheetText) %_% ".png", sep = "/"))))

        r <- plinth::plot_series(x, series, x_var = "voltage", log = "", xlab = "PMT Voltage", ylab = expression(paste(log[10], " CV")), main = paste(basename(experiment_name), sheetText), col = color, lwd = 1, trend = FALSE, segmented = segmented, segmented... = list(breakpoints... = list(h = 3)))

        if (segmented) { # Find optimal voltage by segmented linear model
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
        } else { # Find optimal voltage by 2nd derivative.
          z <- plinth::create_smooth_variables(a, series = NULL, x_var = "voltage", pad_by = 3, keep_interpolated = TRUE, deriv = 2, frfast... = list(smooth = "kernel"), interpolated_derivative = FALSE, verbose = TRUE)

          create_smooth_variablesArgs <- as.list(args(plinth::create_smooth_variables))
          dv <- sprintf(create_smooth_variablesArgs$deriv_suffix_template, 2) # 2 => 2nd derivative
          ll <- create_smooth_variablesArgs$lower_ci_suffix
          ul <- create_smooth_variablesArgs$upper_ci_suffix

          zz <- attr(z, "frfast")

          changepoints_cv <- sapply(series,
            function(i)
            {
              zeros <- zz[[i]][[i %_% dv %_% ll]] < 0.0 & zz[[i]][[i %_% dv %_% ul]] > 0.0
              tot_zeros <- sum(zeros, na.rm = TRUE)

              p <- sapply(rev(seq_along(zeros)),
                function (j)
                {
                  sum(zeros[j:length(zeros)], na.rm = TRUE) / tot_zeros
                }, simplify = TRUE)

              cp <- which(rev(p) <= zeros_threshold)[1]

              #list(x = zz[[i]]$voltage[cp], y = trans_fun(drop(plinth::interpNA(zz[[i]][[i]], "linear"))[cp]))
              list(x = zz[[i]]$voltage[cp], y = approx(x = x$voltage, y = x[[i]], xout = zz[[i]]$voltage)$y[cp])
            }, simplify = FALSE)
        }

        dev_null <- sapply(changepoints_cv, function(a) points(a, col = "black", pch = 4, cex = 1))
        #dev_null <- sapply(changepoints_volts, function(a) points(a, col = "gray50", pch = 3, cex = 1))

        if (save_png) dev.off()
      #}, error = Error) # Could add 'warning = Error'.

      rv <- purrr::map_dfr(changepoints_cv, dataframe); rownames(rv) <- names(changepoints_cv); colnames(rv) <- c("PMT_voltage", "log10_CV")

      if (plot_derivative & !segmented) local({
        zzz <- zz #%>%
          #dplyr::mutate_at(dplyr::vars(2:NCOL(.)), dplyr::funs(!! function(x) trans_fun(x)))

        r <- mapply(series %_% dv, color, series,
          FUN = function(i, color, s)
          {
            if (all(Vectorize(is_invalid)(zzz[[s]][[i]])))
              return (nop())

            plinth::plot_series(zzz[[s]], i, x_var = "voltage", log = "", xlab = "PMT Voltage", ylab = "CV", main = "Second Derivative + 95% CI", col = color, lwd = 1, conf_int = TRUE, trend = FALSE, segmented = FALSE, ylim = NULL)

            plinth::vline(sprintf("%.1f", changepoints_cv[[s]]$x), abline... = list(col = scales::alpha("black", 0.4)))
          }, SIMPLIFY = FALSE)
      })

      currentIndex <<- currentIndex + 1

      #browser()
      rv
    }, simplify = FALSE)

  ## Make Peak 2 report.
  if (!missing(report_dir)) {
    rr <- r; names(rr) <- basename(names(r))

    rio::export(rr, paste(report_dir, basename(report_dir) %_% ".xlsx", sep = "/"), rowNames = TRUE)
  }

  r
}
