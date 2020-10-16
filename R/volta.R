
#' Load Peak 2 FCS files at different PMT voltages
#'
#'@description
#' Loads Peak 2 FCS files to gather up channel & voltage data for "knee" analysis.
#'
#' @param path Path to directory holding the Peak 2 FCS files for analysis.
#' @param list.dirs... An optional list of arguments to [base::list.dirs()], for finer selection control of the subdirectories of `path`.
#' @param list.files... An optional list of arguments to [base::list.files()], for finer selection control of Peak 2 FCS files under `path`.
#' @param read.FCS... An optional list of arguments to [flowCore::read.FCS()], for reading in the Peak 2 FCS files.
#' @param min_zero Logical; if `TRUE`, shift expression values for each channel so that the minimum value is zero. The default is `FALSE`.
#' @param cv_fun Function to calculate the coefficient of variation against which PMT voltage is plotted. The default is BD Biosciences robust standard deviation (BD-rSD), `bd_rcv`; alternatives are `cv` & robust coefficient of variation (rCV), `rcv`.
#' @param cv_multiplier Numeric; the default value of 100 pushes all CVs above 1.
#' @param replace_with_na_condition Optional argument `condition` to [naniar::replace_with_na_all()] to replace e.g. negative CV values with `NA`.
#' @param color_nm_map Named vector map of color names to nm wavelength values, the default being [`plinth::color_nm_map`].
#'
#' @return
#' [get_peak2_data()] returns a list (named after the input FCS files) of data matrices, each of whose first column contains all test voltages, with remaining columns containing Peak 2 CV values for each voltage. Each data matrix in the list has the following custom attributes:
#' \item{`experiment_name`}{Usually the name of a directory holding a collection of per-voltage FCS files for a particular instrument.}
#' \item{`color`}{Color names for each channel.}
#' \item{`wavelength`}{Color wavelengths (nm) for each channel.}
#' \item{`raw_cv`}{The CV data matrix for this instument before removal of e.g. negatives, `NaN`s, & so on.}
#'
#' @examples
#' \dontrun{
#' library(volta)
#' pacman::p_load(magrittr)
#' x <- get_peak2_data("QC/Peak 2 and Baseline Calibrations/2020/March 2020/Data",
#'   list.dirs... = list(pattern = "(Waldorf 030220|Z01)$", recursive = TRUE))
#' }

#' @import data.table
#' @importFrom magrittr %>%
#' @importFrom plinth %_% cordon is_invalid dataframe

#' @export
get_peak2_data <- function(
  path = ".",
  list.dirs... = list(),
  list.files... = list(),
  read.FCS... = list(),
  min_zero = FALSE,
  cv_fun = plinth::bd_rcv,
  cv_multiplier = 100,
  replace_with_na_condition = ~ is.infinite(.x) | is.nan(.x) | .x < 0.0,
  color_nm_map = plinth::color_nm_map
)
{
  list.dirsArgs <- list(
    path = path
  )
  list.dirsArgs <- utils::modifyList(list.dirsArgs, list.dirs...)
  list.dirsPattern <- list.dirsArgs$pattern
  list.dirsArgs$pattern <- NULL
  dirPaths <- do.call(list.dirs, list.dirsArgs)

  if (!is.null(list.dirsPattern))
    dirPaths <- dirPaths[stringr::str_detect(dirPaths, list.dirsPattern)]

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
  colorRe <- sprintf("(%s)", paste(names(color_nm_map), collapse = "|"))

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
              cbind(voltage = k, j[, channel_voltages %in% k])
            }, simplify = FALSE)) %>% dplyr::mutate(voltage = as.numeric(plinth::unfactor(voltage)))

          v
        }, simplify = FALSE)

      raw_cv <- Reduce(plyr::rbind.fill, rr)
      rrr <- raw_cv %>%
        naniar::replace_with_na_all(condition = replace_with_na_condition) %>%
        dplyr::arrange(voltage)

      d <- data.table::data.table(rrr)
      ## Remove voltage duplicates by averaging.
      #d <- d[, lapply(.SD, mean, na.rm = TRUE), by = .(voltage), .SDcols = tail(colnames(d), -1)]
      d <- d[, lapply(.SD, function(j) { if (all(is.na(j))) NA_real_ else mean(j, na.rm = TRUE) }), by = .(voltage), .SDcols = tail(colnames(d), -1)]
      rrr <- d %>% as.data.frame()

      ## Add some attributes.
      color <- tolower(stringr::str_match(colnames(rrr)[2:NCOL(rrr)], stringr::regex(colorRe, ignore_case = TRUE))[, 2])
      wavelength <- stringr::str_match(colnames(rrr)[2:NCOL(rrr)], nmRe)[, 2]
      if (all(Vectorize(is_invalid)(wavelength)))
        wavelength <- color_nm_map[color]
      wavelength <- as.numeric(wavelength)

      attr(rrr, "experiment_name") <- attr(r[[1]], "experiment_name")
      attr(rrr, "color") <- color
      attr(rrr, "wavelength") <- wavelength
      attr(rrr, "raw_cv") <- raw_cv

      rrr
    }, simplify = FALSE)

  l <- p2

  l
}


#' Plot Peak 2 CV vs. PMT voltages & find minimum "knee" voltage
#'
#'@description
#' Creates a plot of Peak 2 CV vs. PMT voltages for each channel for each instrument & finds minimum "knee" voltage for optimal resolution sensitivity. Optionally generates a summary report as an Excel workbook.
#'
#' @param x List object generated by [get_peak2_data]
#' @param report_dir \[Coming soon\]
#' @param image_dir ...
#' @param trans_fun ...
#' @param replace_with_na_condition ...
#' @param remove_empty_cols ...
#' @param save_png ...
#' @param plot_series... ...
#' @param create_smooth_variables... ...
#' @param plot_derivative ...
#' @param points... ...
#' @param debug ...
#' @param xlsx_expression ...
#'
#' @return
#'
#' @examples
#' \dontrun{
#' ## There will totally be code here soon.
#' }

#' @export
plot_peak2_data <- function(
  x,
  report_dir,
  image_dir = report_dir,
  trans_fun = log10, # Also possibly 'identity'
  replace_with_na_condition = ~ is.infinite(.x) | is.nan(.x),
  max_cv_threshold = Inf, # Make anomalously high-valued CVs into missings
  keep_longest_continuous = FALSE,
  remove_empty_cols = TRUE,
  save_png = FALSE, png... = list(),
  plot_series... = list(),
  create_smooth_variables... = list(),
  plot_derivative = FALSE,
  points... = list(),
  debug = FALSE,
  round_to_nearest_volts = 1, # Round to nearest multiple of this value
  xlsx_expression = NULL
)
{
  if (save_png && !dir.exists(image_dir))
    dir.create(image_dir, recursive = TRUE)

  if (!missing(report_dir) && !dir.exists(report_dir))
    dir.create(report_dir, recursive = TRUE)

  pngArgs <- list(
    #width = 12.5,
    width = 9.375,
    height = 7.3,
    units = "in",
    res = 600
  )
  pngArgs <- utils::modifyList(pngArgs, png..., keep.null = TRUE)

  indices <- NULL # NULL for all elements
  if (is.null(indices))
    indices <- seq_along(x)
  currentIndex <- head(indices, 1)

  r <- sapply(x[indices],
    function(a)
    {
      color <- plinth::wavelength2col(attr(a, "wavelength"))
      experiment_name <- attr(a, "experiment_name"); if (is_invalid(experiment_name)) experiment_name <- "[no experiment_name]"

      cat("Processing object # ", currentIndex, ", ", basename(experiment_name), "\n", sep = ""); flush.console()

      y <- a %>%
        naniar::replace_with_na_all(condition = eval(substitute(~ .x > MCT, list(MCT = max_cv_threshold)))) %>%
        dplyr::mutate_at(dplyr::vars(2:NCOL(.)), trans_fun) %>%
        naniar::replace_with_na_all(condition = replace_with_na_condition) %>%
        dplyr::mutate_at(dplyr::vars(2:NCOL(.)),
          function(b)
          {
            if (keep_longest_continuous) {
              ## Keep all values in rightmost stretch
              # ii <- b[b %>% plinth::na_unwrap()] %>% is.na %>% `!` %>% rle %>% `$`("lengths") %>% plinth::cum_sum() %>% `c`(1, .)
              # i <- (ii[length(ii) - 1]):(length(b))
              # bb <- rep(NA_real_, length(b))
              # bb[i] <- b[i]

              ## Keep all values including + after longest stretch
              r <- b[b %>% plinth::na_unwrap()] %>% is.na %>% `!` %>% rle
              i <- ((c(0, r$lengths) %>% cumsum)[r$lengths[r$values] %>% which.max] + 1):(length(b))
              bb <- rep(NA_real_, length(b))
              bb[i] <- b[i]
            }
            else {
              bb <- b
            }

            bb
          })

      if (remove_empty_cols) {
        notAllNas <- sapply(y[, 2:NCOL(y)], function(i) !all(is.na(i)))
        if (!all(notAllNas)) {
          color <- color[notAllNas]
          y <- janitor::remove_empty(y, which = c("cols"))

          warning("Channel(s) ", paste(names(notAllNas)[!notAllNas], sep = ","), " has invalid CVs for all voltages.", immediate. = TRUE)
        }
      }

      series <- names(y)[2:NCOL(y)]

      Error <- function(e) {
        cat("Fatal error: ", e$message, ". Aborting & continuing next analysis.\n", sep = ""); flush.console()
        if (dev.cur() >= 0L) dev.off()
      }

      ### Do time-series plotting.

      #tryCatch({
        if (save_png) do.call(grDevices::png, utils::modifyList(pngArgs, list(filename = paste(image_dir, sprintf("%03d", currentIndex) %_% " - " %_% basename(experiment_name) %_% ".png", sep = "/")), keep.null = TRUE))

        plot_seriesArgs <- list(
          x = y,
          series = series,
          x_var = "voltage",
          log = "",
          xlab = "PMT Voltage", ylab = expression(paste(log[10], " CV")),
          main = basename(experiment_name),
          dev.new... = list(width = 9.375, height = 7.3),
          col = color, lwd = 1,
          trend = FALSE,
          segmented = FALSE, segmented... = list(breakpoints... = list(h = 3)),
          legend... = list(x = "topright")
        )
        plot_seriesArgs <- utils::modifyList(plot_seriesArgs, plot_series..., keep.null = TRUE)
        do.call(plinth::plot_series, plot_seriesArgs)

        ## Find smoothed series & derivative(s).
        deriv <- 1:2
        create_smooth_variablesArgs <- list(
          x = y,
          series = NULL,
          x_var = "voltage",
          pad_by = 3,
          keep_interpolated = TRUE,
          deriv = deriv,
          frfast... = list(smooth = "kernel"),
          interpolated_derivative = FALSE,
          loess... = list(span = 1.0),
          verbose = TRUE
        )
        create_smooth_variablesArgs <- utils::modifyList(create_smooth_variablesArgs, create_smooth_variables..., keep.null = TRUE)
        z <- do.call(plinth::create_smooth_variables, create_smooth_variablesArgs)

        create_smooth_variablesArgs <- as.list(args(plinth::create_smooth_variables))
        o <- structure(vector("list", length = length(deriv)), .Names = as.character(deriv))
        for (i in as.character(deriv)) {
          o[[i]] <- list()

          o[[i]]$dv <- sprintf(create_smooth_variablesArgs$deriv_suffix_template, as.numeric(i))
          o[[i]]$ll <- create_smooth_variablesArgs$lower_ci_suffix
          o[[i]]$ul <- create_smooth_variablesArgs$upper_ci_suffix
        }

        zz <- attr(z, "frfast")

        changepoints_cv <- sapply(series,
          function(i)
          {
            if (debug) {
              ## Radius of curvature
              ## V. https://math.stackexchange.com/questions/690677/is-there-a-name-for-the-point-of-a-exponential-curve-where-the-y-axis-significan/1734110#1734110
              yp <- zz[[i]][[i %_% o[["1"]]$dv]]; ypp <- zz[[i]][[i %_% o[["2"]]$dv]]
              R <- ((1 + yp^2)^(3/2)) / abs(ypp)
              dfR <- dataframe(voltage = zz[[i]]$voltage, R = R)

              plinth::plot_series(zz[[i]], i, x_var = "voltage", xlab = "PMT Voltage", ylab = "CV", main = "CV vs. Voltage", col = color[series == i], lwd = 1, conf_int = TRUE, trend = FALSE, segmented = FALSE, ylim = NULL)
              plinth::plot_series(zz[[i]], i %_% sapply(o, `[[`, "dv"), x_var = "voltage", xlab = "PMT Voltage", ylab = "CV", main = "Derivatives + 95% CI", col = color[series == i], lwd = 1, lty = seq(length(o)), legend... = list(lty = seq(length(o))), conf_int = TRUE, trend = FALSE, segmented = FALSE, ylim = c(-0.0001, 0.0001))
              plinth::plot_series(dfR, "R", x_var = "voltage", xlab = "PMT Voltage", ylab = "Radius of Curvature", main = "Radius of Curvature", col = color[series == i], lwd = 1, conf_int = TRUE, trend = FALSE, segmented = FALSE, ylim = NULL)
            }

            x1 <- zz[[i]]$voltage; y1 <- zz[[i]][[i]]
            checkCurve <- inflection::check_curve(x1, y1)
            if (checkCurve$index == 1) { # This works! I should figure out why....
              bese <- inflection::bese(x1, y1, index = checkCurve$index, doparallel = FALSE)
              j <- x1 >= bese$iplast
              x1 <- x1[j]; y1 <- y1[j]
            }
            knee <- inflection::uik(x1, y1)
            cp <- plinth::nearest(x1, knee[1]) # Could just be 'cp <- knee[1]'.

            r <- list(
              x = zz[[i]]$voltage[cp] %>% plyr::round_any(round_to_nearest_volts),
              y = approx(x = y$voltage, y = y[[i]], xout = zz[[i]]$voltage)$y[cp]
            )

            r
          }, simplify = FALSE)

        pointsArgs <- list(
          col = "black",
          pch = 4, cex = 1
        )
        pointsArgs <- utils::modifyList(pointsArgs, points..., keep.null = TRUE)
        dev_null <- sapply(changepoints_cv, function(a) { pointsArgs$x = a; do.call(points, pointsArgs) })

        if (save_png) dev.off()
      #}, error = Error) # Could add 'warning = Error'.

      rv <- purrr::map_dfr(changepoints_cv, dataframe); rownames(rv) <- names(changepoints_cv); colnames(rv) <- c("PMT_voltage", "log10_CV")

      if (plot_derivative) local({
        zzz <- zz #%>% dplyr::mutate_at(dplyr::vars(2:NCOL(.)), dplyr::funs(!! function(x) trans_fun(x)))

        r <- mapply(combine_groups(list(series, sapply(o, `[[`, "dv")), sep = ""),
          rep(color, each = length(o)), rep(series, each = length(o)),
          FUN = function(i, color, s)
          {
            if (all(Vectorize(is_invalid)(zzz[[s]][[i]])))
              return (nop())

            ## Currently not configurable.
            plinth::plot_series(zzz[[s]], i, x_var = "voltage", log = "", xlab = "PMT Voltage", ylab = "CV", main = "First Derivative + 95% CI", col = color, lwd = 1, conf_int = TRUE, trend = FALSE, segmented = FALSE, ylim = NULL)

            plinth::vline(sprintf("%.1f", changepoints_cv[[s]]$x), abline... = list(col = scales::alpha("black", 0.4)))
          }, SIMPLIFY = FALSE)
      })

      currentIndex <<- currentIndex + 1

      rv <- rv %>%
      tibble::rownames_to_column(var = "channel") %>%
      dplyr::mutate(
        PMT_voltage = PMT_voltage,
        log10_CV = log10_CV %>% round(2)
      )

      rv
    }, simplify = FALSE)

  ## Make Peak 2 report.
  if (!missing(report_dir)) {
    rr <- r; names(rr) <- basename(names(r))

    fileName <- paste(report_dir, basename(report_dir) %_% ".xlsx", sep = "/")
    rio::export(rr, fileName, rowNames = FALSE)

    wb <- xlsx::loadWorkbook(fileName)
    plinth::poly_eval(xlsx_expression)

    ## Add plots to report.
    if (save_png) {
      ss <- xlsx::getSheets(wb)
      imageFiles <- list.files(report_dir, "^\\d{3} - .*?\\.png", full.names = TRUE, ignore.case = TRUE)

      dev_null <- sapply(seq_along(ss), function(i) { xlsx::addPicture(imageFiles[i], ss[[i]], scale = 1, startRow = 1, startColumn = 4); nop() })

      xlsx::saveWorkbook(wb, fileName)
    }
  }

  r
}
