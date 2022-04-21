
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
#' @importFrom flowpipe `colnames<-`

## FCS file format: https://en.wikipedia.org/wiki/Flow_Cytometry_Standard

#' @export
get_peak2_data <- function(
  path = ".",
  ## These calcs typically rely on data stored in ad-hoc keywords for parameter n, e.g. PnV (_V_oltage), PnG (_G_ain):
  cv_keyword_fmt = "$P%sV",
  ## 'analysis_channels_re' is one of: 1-element regex/character vector; 2+-element vector of channel names
  analysis_channels_re = stringr::regex("^(?!(fsc|ssc|time).*)", ignore_case = TRUE),
  ## 'use_spillover_channels' overrides 'analysis_channels_re' if TRUE; if so, use channel names from spillover matrix
  use_spillover_channels = FALSE,
  create_pmm_data = FALSE,
  prepare_augmented_fcs_data... = list(),
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
      filePaths <- do.call(plinth::list_files, list.filesArgs)

      filePaths
    }, simplify = FALSE) %>% purrr::compact() # Remove blank elements

  ## N.B. TODO: Allow review of directories & files here w/ interactive go/no-go.
  #browser()
  nmRe <- "\b?([2-9][0-9][0-9])\b?"
  colorRe <- sprintf("(%s)", paste(names(color_nm_map), collapse = "|"))

  p2 <- sapply(ff,
    function(i)
    {
      ii <- i
      if (create_pmm_data) {
        tictoc::tic("Make PMM files")

        ## Make augmented flowCore::flowFrame objects.
        pmm_files <- plinth::psapply(i,
          function(j)
          {
            prepare_augmented_fcs_dataArgs <- list(
              x = j,
              b = 1/150,
              data_dir = sprintf("./data/%s", j %>% dirname %>% basename),
              remove_outliers = FALSE,
              excluded_channels_re = stringr::regex("time|event_length|(width$)", ignore_case = TRUE),
              multisect... = list(max_sample = 5000),
              outfile_suffix = NULL,
              outfile_prefix = expression(outfile_prefix <- paste0(x %>% dirname %>% basename, "-")),
              overwrite = FALSE
            )
            prepare_augmented_fcs_dataArgs <- utils::modifyList(prepare_augmented_fcs_dataArgs, prepare_augmented_fcs_data..., keep.null = TRUE)

            pmm_file <- do.call(flowpipe::prepare_augmented_fcs_data, prepare_augmented_fcs_dataArgs)

            pmm_file
          }, simplify = FALSE)

        ii <- pmm_files

        tictoc::toc()
      }

      r <- sapply(ii,
        function(j)
        {
          fcs <- stain_groups <- NULL
          if (!create_pmm_data) {
            local({
              read.FCSArgs <- list(
                filename = j,
                truncate_max_range = FALSE,
                transformation = FALSE
              )
              read.FCSArgs <- utils::modifyList(read.FCSArgs, read.FCS...)
              fcs <<- do.call(flowCore::read.FCS, read.FCSArgs)
            })
          } else {
            local({
              load(j)
              fcs <<- rlang::duplicate(tff, shallow = FALSE)
              remove(tff)

              exprs_fcs <- flowCore::exprs(fcs)
              pmm <- attr(exprs_fcs, "plus_minus_matrix")

              ## Separate events as stained/unstained for each channel
              stain_groups <<- sapply(colnames(pmm),
                function(a)
                {
                  pm <- pmm %>% dplyr::select(!!a) %>% dplyr::pull()
                  if (!is.factor(pm))
                    return (NULL)
                  pm <- pm %>%
                    forcats::fct_collapse(
                      unstained = c("-", "d"),
                      stained = c("+", "++")
                    )

                  r <- exprs_fcs[, a] %>% split(pm)
                  if (min_zero)
                    r <- sapply(r, function(k) k - min(k, na.rm = TRUE))

                  r
                }, simplify = FALSE) %>% purrr::compact()
            })
          }

          keywords <- fcs %>% flowCore::keyword()
          pData <- fcs %>% flowCore::parameters() %>% flowCore::pData() %>%
            tibble::rownames_to_column(var = "parameter") %>%
            dplyr::mutate(
              quantity_keyword = sprintf(cv_keyword_fmt, stringr::str_extract(parameter, "\\d+$")),
              desc = dplyr::case_when(is.na(desc) ~ name, TRUE ~ desc)
            ) %>%
            dplyr::left_join(
              plinth::dataframe(
                parameter = .$parameter,
                quantity =
                  sapply(.$quantity_keyword,
                    function(a) { r <- keywords[[a]]; if (is_invalid(r)) r <- NA; r })
              ), by = "parameter"
            ) %>%
            dplyr::mutate(quantity = quantity %>% readr::parse_number())

          ## 'pData' should now look something like this:
          #   parameter              name              desc  range minRange maxRange quantity_keyword quantity
          # 1        $P1             FSC-A             FSC-A 262144     0.00   262143             $P1V      450
          # 2        $P2             SSC-A             SSC-A 262144     0.00   262143             $P2V      200
          # 3        $P3   Blue B 515/20-A   Blue B 515/20-A 262144   -63.84   262143             $P3V      250
          # 4        $P4   Blue A 710/50-A   Blue A 710/50-A 262144  -111.00   262143             $P4V      250
          # 5        $P5 Violet H 450/50-A Violet H 450/50-A 262144   -40.18   262143             $P5V      250
          # 6        $P6 Violet G 550/40-A Violet G 550/40-A 262144   -48.02   262143             $P6V      250
          # 7        $P7 Violet F 560/40-A Violet F 560/40-A 262144   -60.76   262143             $P7V      250
          # 8        $P8 Violet E 585/42-A Violet E 585/42-A 262144   -48.02   262143             $P8V      250
          # 9        $P9 Violet D 605/40-A Violet D 605/40-A 262144   -58.80   262143             $P9V      250
          # 10      $P10 Violet C 660/40-A Violet C 660/40-A 262144   -56.84   262143            $P10V      250
          # 11      $P11 Violet B 705/70-A Violet B 705/70-A 262144   -73.50   262143            $P11V      250
          # 12      $P12 Violet A 780/60-A Violet A 780/60-A 262144   -58.80   262143            $P12V      250
          # 13      $P13    Red C 660/20-A    Red C 660/20-A 262144   -34.16   262143            $P13V      250
          # 14      $P14    Red B 710/50-A    Red B 710/50-A 262144   -31.11   262143            $P14V      250
          # 15      $P15    Red A 780/60-A    Red A 780/60-A 262144   -44.53   262143            $P15V      250
          # 16      $P16  Green E 575/25-A  Green E 575/25-A 262144   -47.30   262143            $P16V      250
          # 17      $P17  Green D 610/20-A  Green D 610/20-A 262144   -38.70   262143            $P17V      250
          # 18      $P18  Green C 660/40-A  Green C 660/40-A 262144   -49.88   262143            $P18V      250
          # 19      $P19  Green B 710/50-A  Green B 710/50-A 262144   -48.16   262143            $P19V      250
          # 20      $P20  Green A 780/40-A  Green A 780/40-A 262144   -47.30   262143            $P20V      250
          # 21      $P21              Time              Time 262144     0.00   262143            $P21V       NA

          ## Now pick out the relevant channels.
          ## Often, the spillover matrix uses those channel names:
          # names(keywords) %>% stringr::str_subset(stringr::regex("spill", ignore_case = TRUE)) %>% `[[`(keywords, .) %>% colnames

          ## Check for spillover matrix
          spilloverKey <- names(keywords) %>% stringr::str_subset(stringr::regex("spill", ignore_case = TRUE))
          hasSpilloverMatrix <- spilloverKey %>% `[[`(keywords, .) %>% is.matrix
          spilloverChannels <- NULL
          if (hasSpilloverMatrix)
            spilloverChannels <- spilloverKey %>% `[[`(keywords, .) %>% colnames

          ## Assume that 'pData$desc' has the channel info of interest:
          pData$desc[is.na(pData$desc)] <- pData$name[is.na(pData$desc)]
          description <- fcs %>% flowCore::description()

          if (!is.null(use_spillover_channels) && use_spillover_channels) {
            if (!hasSpilloverMatrix) {
              warning("Spillover matrix not found. Defaulting to 'analysis_channels_re'", immediate. = FALSE)

              cvChannels <- NULL
            } else {
              cvChannels <- spilloverChannels
            }
          }

          if (is.null(use_spillover_channels) || !use_spillover_channels || is.null(cvChannels)) {
            cvChannels <- stringr::str_subset(pData$desc, analysis_channels_re)
          }

          descNameMapFull <- structure(pData$desc, .Names = pData$name)
          descNameMap <- descNameMapFull[descNameMapFull %in% cvChannels]

          e <- flowCore::exprs(fcs) %>% `colnames<-`(descNameMapFull[colnames(.)]) %>%
            `[`(, descNameMap)
          if (min_zero)
            e <- plyr::aaply(e, 1, function(k) k - min(k, na.rm = TRUE))

          if (!is.null(stain_groups))
            names(stain_groups) <- descNameMapFull[names(stain_groups)]

          ## To calculate Voltration index (VI) or anything involving stained-vs.-unstained comparisons,
          ##   we need to use the distinct populations in 'stain_groups'.

          quantityMap <- structure(pData$quantity, .Names = pData$desc)
          # ee <- e %>% plyr::aaply(2, function(k) { cv_fun(data.matrix(k)) * cv_multiplier }) %>% t %>% plinth::dataframe()
          ee <- sapply(colnames(e),
            function(a) { cv_fun(unclass(e[, a]), stain_groups = stain_groups[[a]], quantity = quantityMap[a][1]) * cv_multiplier },
            simplify = FALSE) %>% plinth::dataframe()
          attr(ee, "channel_quantities") <- structure(pData$quantity, .Names = pData$desc)[colnames(e)]
          attr(ee, "machine_name") <- keywords$`$CYT`
          attr(ee, "experiment_name") <- keywords$FILENAME %>% dirname %>% basename

          ee
        }, simplify = FALSE)

      ## It probably won't happen, but this can separate channels that aren't on the same voltage for some reason.
      rr <- sapply(r,
        function(j)
        {
          channel_quantities <- attr(j, "channel_quantities")
          v <- Reduce(plyr::rbind.fill, sapply(unique(channel_quantities),
            function(k)
            {
              cbind(quantity = k, j[, channel_quantities %in% k])
            }, simplify = FALSE)) %>% dplyr::mutate(quantity = as.numeric(plinth::unfactor(quantity)))

          v
        }, simplify = FALSE)

      raw_cv <- Reduce(plyr::rbind.fill, rr)
      rrr <- raw_cv %>%
        naniar::replace_with_na_all(condition = replace_with_na_condition) %>%
        dplyr::arrange(quantity)

      d <- data.table::data.table(rrr)
      ## Remove voltage duplicates by averaging.
      #d <- d[, lapply(.SD, mean, na.rm = TRUE), by = .(quantity), .SDcols = tail(colnames(d), -1)]
      d <- d[, lapply(.SD, function(j) { if (all(is.na(j))) NA_real_ else mean(j, na.rm = TRUE) }), by = .(quantity), .SDcols = tail(colnames(d), -1)]
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
  report_dir = NULL,
  image_dir = report_dir,
  trans_fun = log10, # Also possibly 'identity'
  replace_with_na_condition = ~ is.infinite(.x) | is.nan(.x),
  max_cv_threshold = Inf, # Make anomalously high-valued CVs into missings
  keep_longest_continuous = FALSE,
  remove_empty_cols = TRUE,
  save_png = FALSE, png... = list(),
  default_color_fun = colorspace::rainbow_hcl,
  x_var_lab = c(PMT_voltage = "PMT Voltage"), y_var_lab = c(log10_CV = expression(paste(log[10], " CV"))),
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

  if (!is.null(report_dir) && !dir.exists(report_dir))
    dir.create(report_dir, recursive = TRUE)

  pngArgs <- list(
    #width = 12.5,
    width = 9.375,
    height = 7.3,
    units = "in",
    res = 600
  )
  pngArgs <- utils::modifyList(pngArgs, png..., keep.null = TRUE)

  x_var_lab <- head(x_var_lab, 1)
  if (is_invalid(names(x_var_lab)) || trimws(names(x_var_lab)) == "") names(x_var_lab) <- x_var_lab
  y_var_lab <- head(y_var_lab, 1)
  if (is_invalid(names(y_var_lab)) || trimws(names(y_var_lab)) == "") names(y_var_lab) <- y_var_lab

  ## Create plotting data set
  r <- sapply(x,
    function(a)
    {
      color <- plinth::wavelength2col(attr(a, "wavelength"))
      defaultColor <- default_color_fun(length(color))
      color[is.na(color)] <- defaultColor[is.na(color)]
      experiment_name <- attr(a, "experiment_name"); if (is_invalid(experiment_name)) experiment_name <- "[no experiment_name]"

      cat(sprintf("Processing object: %s\n", experiment_name)); flush.console()

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

          warning("Channel(s) ", paste(names(notAllNas)[!notAllNas], sep = ","), " has invalid CVs for all quantities.",
            immediate. = TRUE)
        }
      }

      series <- names(y)[2:NCOL(y)]

      Error <- function(e) {
        cat("Fatal error: ", e$message, ". Aborting & continuing next analysis.\n", sep = ""); flush.console()
        if (dev.cur() >= 0L) dev.off()
      }

      ### Do time-series plotting.

      #tryCatch({
        plot_seriesArgs <- list(
          x = y,
          series = series,
          x_var = "quantity",
          log = "",
          xlab = x_var_lab, ylab = y_var_lab,
          main = experiment_name,
          dev.new... = list(width = 9.375, height = 7.3),
          col = color, lwd = 1,
          trend = FALSE,
          segmented = FALSE, segmented... = list(breakpoints... = list(h = 3)),
          legend... = list(x = "topright")
        )

        ## Find smoothed series & derivative(s).
        deriv <- 1:2
        create_smooth_variablesArgs <- list(
          x = y,
          series = NULL,
          x_var = "quantity",
          pad_by = 3,
          keep_interpolated = TRUE,
          deriv = deriv,
          frfast... = list(smooth = "kernel"),
          interpolated_derivative = FALSE,
          loess... = list(span = 1.0),
          verbose = TRUE
        )
        create_smooth_variablesArgs <-
          utils::modifyList(create_smooth_variablesArgs, create_smooth_variables..., keep.null = TRUE)
        tictoc::tic("Create smoothed variables")
        z <- do.call(plinth::create_smooth_variables, create_smooth_variablesArgs)
        tictoc::toc()

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
            x1 <- zz[[i]]$quantity; y1 <- zz[[i]][[i]]
            if (any(is.na(y1))) { # This should mostly be false, but isn't always
              warning("Smoothed variable 'y1' contains some missings.", immediate. = TRUE)

              y1 <- drop(plinth::interpNA(y1, method = "fmm", unwrap = FALSE))
            }
            checkCurve <- inflection::check_curve(x1, y1)
            if (checkCurve$index == 1) { # This works! I should figure out why....
              bese <- inflection::bese(x1, y1, index = checkCurve$index, doparallel = FALSE)
              j <- x1 >= bese$iplast
              x1 <- x1[j]; y1 <- y1[j]
            }
            knee <- inflection::uik(x1, y1)
            cp <- plinth::nearest(x1, knee[1]) # Could just be 'cp <- knee[1]'.

            r <- list(
              x = zz[[i]]$quantity[cp] %>% plyr::round_any(round_to_nearest_volts),
              y = approx(x = y$quantity, y = y[[i]], xout = zz[[i]]$quantity)$y[cp]
            )

            r
          }, simplify = FALSE)
      #}, error = Error) # Could add 'warning = Error'.

      rv <- purrr::map_dfr(changepoints_cv, dataframe)
      rownames(rv) <- names(changepoints_cv)
      colnames(rv) <- c(names(x_var_lab), names(y_var_lab))

      if (plot_derivative) {} # See GitHub history for reimplementation if needed.

      rv <- rv %>%
        tibble::rownames_to_column(var = "channel") %>%
        dplyr::mutate(
          !!names(x_var_lab) := .data[[names(x_var_lab)]],
          !!names(y_var_lab) := .data[[names(y_var_lab)]] %>% round(2)
        )

      structure(
        list(plot_args = plot_seriesArgs, changepoints_cv = changepoints_cv, table = rv, experiment_name = experiment_name),
        derivatives = tibble::as_tibble(zz)
      )
    }, simplify = FALSE)

  ## Create plots
  plyr::l_ply(seq_along(r),
    function(a)
    {
      if (save_png) {
        do.call(grDevices::png,
          utils::modifyList(pngArgs,
            ## N.B. 'sprintf(0)' returns 0-length string for any NULL values; use 'format(NULL)' to output "NULL".
            list(filename = sprintf("%s/%03d - %s.png", format(image_dir), a, basename(r[[a]]$experiment_name))),
            keep.null = TRUE)
        )
      }

      do.call(plinth::plot_series, r[[a]]$plot_args)

      pointsArgs <- list(
        col = "black",
        pch = 4, cex = 1
      )
      pointsArgs <- utils::modifyList(pointsArgs, points..., keep.null = TRUE)
      plyr::l_ply(r[[a]]$changepoints_cv,
        function(b) { pointsArgs$x = b; do.call(points, pointsArgs) })

      if (save_png) dev.off()
    })

  ## Make Peak 2 report.
  if (!is.null(report_dir)) {
    rr <- sapply(r, function(a) { a$table }, simplify = FALSE)
    names(rr) <- stringr::str_trunc(basename(names(r)), 29, "center")

    duplicateNames <- names(rr) %>% intersect(.[duplicated(.)])
    for(i in duplicateNames) {
      dupIndex <- which(names(rr) == i)
      # Replace w/ sequential numbers:
      names(rr)[dupIndex] <-
        sapply(seq_along(dupIndex),
          function(j) sprintf("%s_%01d", names(rr)[dupIndex[j]], j))
    }

    fileName <- paste(report_dir, basename(report_dir) %_% ".xlsx", sep = "/")
    rio::export(rr, fileName, rowNames = FALSE)

    wb <- xlsx::loadWorkbook(fileName)
    plinth::poly_eval(xlsx_expression)

    ## Add plots to report.
    if (save_png) {
      ss <- xlsx::getSheets(wb)
      imageFiles <- list.files(report_dir, "^\\d{3} - .*?\\.png", full.names = TRUE, ignore.case = TRUE)

      plyr::l_ply(seq_along(ss),
        function(i) { xlsx::addPicture(imageFiles[i], ss[[i]], scale = 1, startRow = 1, startColumn = 4) })

      xlsx::saveWorkbook(wb, fileName)
    }
  }

  r
}
