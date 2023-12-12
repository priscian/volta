#' @import data.table
#' @importFrom magrittr %>% %<>%
#' @importFrom keystone %_% cordon is_invalid dataframe
#' @importFrom flowpipe `colnames<-`

## FCS file format: https://en.wikipedia.org/wiki/Flow_Cytometry_Standard

expand_parameters_data <- function(
  x,
  analysis_channels_sep = ",", quantity_sep = analysis_channels_sep,
  same_word = "SAME"
)
{
  pp <- dplyr::group_by(x, machine) %>%
    { structure(dplyr::group_split(., .keep = TRUE), .Names = dplyr::group_keys(.) %>%
      dplyr::pull(machine)) }

  ep <- sapply(pp,
    function(a)
    {
      ## Convert "serialized" channels strings to list of character vectors
      channels_ex <- a$channels
      for (i in seq(length(channels_ex))) {
        if (channels_ex[i] == same_word)
          channels_ex[i] <- channels_ex[i - 1]
      }
      channels_ex %<>%
        sapply(stringr::str_split_1, pattern = analysis_channels_sep,
          USE.NAMES = FALSE, simplify = FALSE)

      a %<>% dplyr::mutate(channels_ex = channels_ex, .after = "channels")

      ## Convert "serialized" quantity strings to list of character vectors
      quantity_ex <- mapply(a$quantity, a$channels_ex,
        FUN = function(d, e)
        {
          qex <- stringr::str_split_1(d, pattern = quantity_sep) %>%
            rep(length.out = length(e))

          qex
        }, USE.NAMES = FALSE, SIMPLIFY = FALSE)

      a %<>% dplyr::mutate(quantity_ex = quantity_ex, .after = "quantity")

      a
    }, simplify = FALSE)

  ep
}


#' @export
get_voltration_data <- function(
  x, # volta analysis parameters as a data frame or spreadsheet path
  create_si_data = FALSE,
  read.FCS... = list(),
  min_zero = FALSE,
  ## Default CV is BD Biosciences' robust coefficient of variation (BD-rCV):
  cv_fun = keystone::bd_rcv,
  cv_multiplier = 100,
  replace_with_na_condition = ~ is.infinite(.x) | is.nan(.x) | .x < 0.0,
  copy_results = TRUE,
  expand_parameters_data... = list(),
  ... # Arguments for function 'find_inflection_point()'
)
{
  if (missing(x) || is_invalid(x)) {
    params <- .volta$params
  } else if (is.character(x)) {
    params <- rio::import(x)
  }

  ## Expand parameters data to streamline voltration analysis
  # ep <- expand_parameters_data(params, ...)
  expand_parameters_dataArgs <- list(
    x = params
  )
  expand_parameters_dataArgs <-
    utils::modifyList(expand_parameters_dataArgs, expand_parameters_data..., keep.null = TRUE)

  ep <- do.call(expand_parameters_data, expand_parameters_dataArgs)

  p2 <- keystone::psapply(ep,
    function(i)
    {
      r <- plyr::alply(i, 1,
        function(j)
        {
          read.FCSArgs <- list(
            filename = j$path,
            truncate_max_range = FALSE,
            ## Retrieve only 1 event, but keep parameters, keywords, & value from FCS file's header:
            #which.lines = 1,
            transformation = FALSE
          )
          read.FCSArgs <- utils::modifyList(read.FCSArgs, read.FCS...)
          fcs <- do.call(flowCore::read.FCS, read.FCSArgs)
          exprs_tff <- flowpipe::get_fcs_expression_subset(list(fcs), sample = Inf)
          ## N.B. This is extra complicated for historical reasons:
          exprs_tff %<>% { `[`(., , -which(colnames(.) == "id"), drop = FALSE) } %>%
            `[`(, intersect(colnames(.), j$channels_ex[[1]]), drop = FALSE)

          stain_groups <- NULL
          if (create_si_data) {
            local({
              ## Separate events as stained/unstained for each channel
              stain_groups <<- sapply(colnames(exprs_tff),
                function(k)
                {
                  channel_data <- flowCore::exprs(fcs)[, k]
                  #channel_data <- flowpipe::rev_asinh(exprs_tff[, k], b = 1/150) # Should be same as previous
                  #channel_data <- exprs_tff[, k]
                  show_plots <- FALSE # For debugging
                  cutoffs <- multisect(channel_data, bins = 2, plot_cutoff = show_plots)
                  dens <- stats::density(channel_data)
                  fp_unstained <- keystone::find_peaks(dens$x[dens$x < cutoffs], dens$y[dens$x < cutoffs],
                    height_divisor = 2, truncate_overlaps = TRUE, force_results = TRUE, plot = show_plots) %>%
                    dplyr::arrange(desc(y)) %>% `[`(1, )
                  fp_stained <- keystone::find_peaks(dens$x[dens$x >= cutoffs], dens$y[dens$x >= cutoffs],
                    height_divisor = 2, truncate_overlaps = TRUE, force_results = TRUE, plot = show_plots) %>%
                    dplyr::arrange(desc(y)) %>% `[`(1, )
                  fp <- dplyr::bind_rows(fp_unstained, fp_stained) %>%
                    `rownames<-`(c("unstained", "stained"))
                  ## If 'channel_data' was transformed, check (partially) transformed 'fp' (should be close):
                  # fp %>% dplyr::mutate(across(all_of(c("x", "w_min", "w_max")), flowpipe::rev_asinh, b = 1/150))

                  fp
                }, simplify = FALSE) %>% purrr::compact()
            })
          }

          e <- flowCore::exprs(fcs) %>%
            `[`(, intersect(colnames(.), j$channels_ex[[1]]), drop = FALSE)
          if (min_zero)
            e <- plyr::aaply(e, 1, function(k) k - min(k, na.rm = TRUE))

          quantityMap <- structure(j$quantity_ex[[1]], .Names = j$channels_ex[[1]])
          ee <- sapply(colnames(e),
            function(a) { cv_fun(unclass(e[, a]), stain_groups = stain_groups[[a]],
              quantity = quantityMap[a][1]) * cv_multiplier },
            simplify = FALSE) %>% keystone::dataframe()
          attr(ee, "channel_quantities") <- quantityMap[colnames(e)]
          attr(ee, "experiment_name") <- j$machine

          ee
        })

      ## Low probability, but this separates channels not on the same voltage for some reason
      rr <- sapply(r,
        function(j)
        {
          channel_quantities <- attr(j, "channel_quantities")
          v <- Reduce(plyr::rbind.fill, sapply(unique(channel_quantities),
            function(k)
            {
              cbind(quantity = k, j[, channel_quantities %in% k])
            }, simplify = FALSE)) %>%
            dplyr::mutate(quantity = as.numeric(keystone::unfactor(quantity)))

          v
        }, simplify = FALSE)

      raw_cv <- Reduce(plyr::rbind.fill, rr)
      raw_cv_list <- list(raw_cv)

      rv <- sapply(raw_cv_list,
        function(j)
        {
          rrr <- j %>%
            naniar::replace_with_na_all(condition = replace_with_na_condition) %>%
            dplyr::arrange(quantity)

          d <- data.table::data.table(rrr)
          ## Remove voltage duplicates by averaging
          d <-
            d[, lapply(.SD, function(k) { if (all(is.na(k))) NA_real_ else mean(k, na.rm = TRUE) }),
            by = .(quantity), .SDcols = tail(colnames(d), -1)]
          rrr <- d %>% as.data.frame()

          ## Add some attributes.
          gcc <- guess_channel_colors(colnames(rrr)[2:NCOL(rrr)])

          attr(rrr, "experiment_name") <- attr(r[[1]], "experiment_name")
          attr(rrr, "color") <- gcc$color
          attr(rrr, "wavelength") <- gcc$wavelength
          attr(rrr, "raw_cv") <- j

          rrr
        }, simplify = FALSE)

      if (inherits(rv, "list") && length(rv) == 1L)
        rv <- rv[[1]]

      rv
    }, simplify = FALSE)#, .parallel = FALSE)

  fip <- find_inflection_point(p2, ...)

  p2 <- sapply(names(p2),
    function(i) { structure(p2[[i]], plot_data = fip[[i]]) %>% keystone::add_class("volta") })

  ## Copy 'p2' to the package global variable
  if (copy_results) {
    current <- get(".volta", envir = asNamespace("volta"))
    current$results <- p2
    assign(".volta", current, envir = asNamespace("volta"))
  }

  if (interactive()) {
    msg <- paste0(
r"---{
The volta results are ready to plot. To continue the analysis, type:

}---",
      "plot_voltration_data(save_png = TRUE)")
    message(msg); utils::flush.console()
  }

  invisible(p2)
}


find_inflection_point <- function(
  x, # 'p2' list of tables from 'get_voltration_data()'
  trans_fun = log10, # Also possibly 'identity'
  replace_with_na_condition = ~ is.infinite(.x) | is.nan(.x),
  max_cv_threshold = Inf, # Make anomalously high-valued CVs into missings
  ## If TRUE, keep only longest stretch of CVs without missing values:
  keep_longest_continuous = FALSE,
  remove_empty_cols = TRUE,
  create_smooth_variables... = list(),
  round_to_nearest_volts = 10, # Round to nearest multiple of this value
  round_any_fun = round,
  keep_derivatives = FALSE,
  ...
)
{
  ## Create tweaked voltration data set that includes curve inflection points
  r <- keystone::psapply(x,
    function(a)
    {
      color <- attr(a, "color")
      experiment_name <- attr(a, "experiment_name")

      cat(sprintf("Processing object: %s\n", experiment_name)); flush.console()

      y <- a %>%
        naniar::replace_with_na_all(condition =
          eval(substitute(~ .x > MCT, list(MCT = max_cv_threshold)))) %>%
        dplyr::mutate_at(dplyr::vars(2:NCOL(.)), trans_fun) %>%
        naniar::replace_with_na_all(condition = replace_with_na_condition) %>%
        dplyr::mutate_at(dplyr::vars(2:NCOL(.)),
          function(b)
          {
            if (keep_longest_continuous) {
              ## Keep all values in rightmost stretch
              # ii <- b[b %>% keystone::na_unwrap()] %>% is.na %>% `!` %>% rle %>%
              #   `$`("lengths") %>% keystone::cum_sum() %>% `c`(1, .)
              # i <- (ii[length(ii) - 1]):(length(b))
              # bb <- rep(NA_real_, length(b))
              # bb[i] <- b[i]

              ## Keep all values including + after longest stretch
              r <- b[b %>% keystone::na_unwrap()] %>% is.na %>% `!` %>% rle
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

          warning("Channel(s) ", paste(names(notAllNas)[!notAllNas], sep = ","),
            " has invalid CVs for all quantities.", immediate. = TRUE)
        }
      }

      series <- names(y)[2:NCOL(y)]

      Error <- function(e) {
        cat("Fatal error: ", e$message, ". Aborting & continuing next analysis.\n", sep = ""); flush.console()
        if (dev.cur() >= 0L) dev.off()
      }

      ### Find inflection point in voltration curve.

      #tryCatch({
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
        z <- do.call(keystone::create_smooth_variables, create_smooth_variablesArgs)
        tictoc::toc()

        create_smooth_variablesArgsAbo <- as.list(args(keystone::create_smooth_variables))
        o <- structure(vector("list", length = length(deriv)), .Names = as.character(deriv))
        for (i in as.character(deriv)) {
          o[[i]] <- list()

          o[[i]]$dv <- sprintf(create_smooth_variablesArgsAbo$deriv_suffix_template, as.numeric(i))
          o[[i]]$ll <- create_smooth_variablesArgsAbo$lower_ci_suffix
          o[[i]]$ul <- create_smooth_variablesArgsAbo$upper_ci_suffix
        }

        zz <- attr(z, "frfast")

        changepoints_cv <- sapply(series,
          function(i)
          {
            x1 <- zz[[i]][[create_smooth_variablesArgs$x_var]]; y1 <- zz[[i]][[i]]
            if (any(is.na(y1))) { # This should mostly be false, but isn't always
              warning("Smoothed variable 'y1' contains some missings.", immediate. = TRUE)

              y1 <- drop(keystone::interpNA(y1, method = "fmm", unwrap = FALSE))
            }
            checkCurve <- inflection::check_curve(x1, y1)
            # plot(x1, y1, main = sprintf("%s - %s", experiment_name, i))
            # mtext(sprintf("%s; index = %s; %s", checkCurve$ctype, checkCurve$index,
            #   inflection::bese(x1, y1, index = checkCurve$index, doparallel = FALSE)$iplast)) # For debugging
            if (checkCurve$index == 1) { # This works! I should figure out why....
              bese <- inflection::bese(x1, y1, index = checkCurve$index, doparallel = FALSE)
              if (!is_invalid(bese$iplast)) {
                j <- x1 >= bese$iplast
                x1 <- x1[j]; y1 <- y1[j]
              }
            }
            knee <- inflection::uik(x1, y1)
            cp <- keystone::nearest(x1, knee[1]) # Could just be 'cp <- knee[1]'.

            r <- list(
              x = zz[[i]][[create_smooth_variablesArgs$x_var]][cp] %>%
                plyr::round_any(round_to_nearest_volts, f = round_any_fun),
              y = approx(x = y[[create_smooth_variablesArgs$x_var]], y = y[[i]],
                xout = zz[[i]][[create_smooth_variablesArgs$x_var]])$y[cp]
            )

            r
          }, simplify = FALSE)
        #changepoints_cv <- rep(list(list(x = 0, y = 0)), length(series))
      #}, error = Error) # Could add 'warning = Error'.

      rv <- purrr::map_dfr(changepoints_cv, dataframe) %>%
        dplyr::mutate(channel = names(changepoints_cv), .before = 1) %>%
        dplyr::rename(inflection = "x")

      derivatives <- NULL
      if (keep_derivatives)
        derivatives <- sapply(zz, tibble::as_tibble, simplify = FALSE)

      structure(list(time_series = y, inflection_points = rv),
        experiment_name = experiment_name, derivatives = derivatives)
    }, simplify = FALSE)#, .parallel = FALSE)

  r
}


guess_channel_colors <- function(
  x, # Vector of channel names, possibly w/ some laser-color info,
  wavelength_re = formals(guess_parameters)$voltage_re,
  col_fun = colorspace::rainbow_hcl
)
{
  default_color_values <- col_fun(length(x))
  default_cols <- default_color_values %>% grDevices::col2rgb() %>% plyr::alply(2, identity) %>%
    `names<-`(NULL) %>% sapply(function(a) keystone::rgb2name(a["red"], a["green"], a["blue"]),
      USE.NAMES = FALSE) %>% names

  wavelengths <- stringr::str_extract(x, wavelength_re) %>% readr::parse_number()

  cols <- Vectorize(keystone:::nm_to_rgb, "wavelength", SIMPLIFY = FALSE)(wavelengths) %>%
    sapply(function(a) { if (is_invalid(a)) return (NA_character_); keystone::rgb2name(a$R, a$G, a$B) })
  names(cols)[is.na(cols)] <- default_cols[is.na(cols)]
  cols <- names(cols)

  list(wavelength = wavelengths, color = cols)
}


## N.B. This is best run interactively
#' @export
prepare_data <- function(
  x,
  ..., # Passed on to 'guess_parameters()'
  choose_excel = FALSE,
  copy_params = TRUE
)
{
  if (missing(x) || is_invalid(x)) {
    if (interactive()) {
      if (!choose_excel) {
        msg <-
r"---{
Choose a top-level directory containing subdirectories named after the machines
whose PMT voltages you wish to optimize; those subdirectories should in turn
contain FCS files collected over the voltration range for each machine.
}---"
        message(msg); utils::flush.console()
        #x <- keystone::choose_dir()
        x <- svDialogs::dlg_dir()$res

        if (is_invalid(x))
          stop("volta is unable to find a directory with FCS files")
      } else {
        msg <-
r"---{
Choose an Excel (XLS, XLSX) or CSV file that's been correctly prepared with
volta parameters (machine, quantity, channels, path) for analyzing a set of FCS
files collected over the voltration range for each machine.
}---"
        message(msg); utils::flush.console()
        #x <- keystone::choose_files("xlsx", "xls", "csv")
        x <- svDialogs::dlg_open(title = "Open volta parameters sheet",
          filters = svDialogs::dlg_filters[c("xls", "csv"), ])$res

        if (is_invalid(x))
          stop("volta is unable to find a spreadsheet file")
      }
    } else {
      stop("A valid 'x' value must be provided to run this function non-interactively")
    }
  }

  ## Is 'x[1]' a valid path to an XLSX file? If so, skip parameter guessing.
  is_spreadsheet <- FALSE
  if (is.character(x) &&
    stringr::str_detect(stringr::str_trim(tools::file_ext(x[1])),
      stringr::regex("^(xlsx|xls|csv)$", ignore_case = TRUE))) {
    is_spreadsheet <- TRUE
    x <- x[1]
  }

  if (!is_spreadsheet)
    params <- guess_parameters(path = x, ...)
  else
    params <- rio::import(x)

  ## Have interactive user review 'params':
  if (interactive()) {
    msg <-
r"---{
------------------------------------------------------
The volta parameters sheet for this analysis is ready.
------------------------------------------------------
You can:

  • REVIEW or edit the parameters sheet
  • SAVE sheet to file system for review or external editing
  • CONTINUE volta analysis using this parameters sheet
  • QUIT volta for now
}---"

    message(msg); utils::flush.console()
    repeat {
      choice <- keystone::ask_multiple_choice(c("Review", "Save", "Continue", "Quit"))

      ## React to each possible choice:
      switch(choice$lowercase,
        quit = {
          return (invisible(NULL))
        },

        review = {
          msg <-
r"---{
After making changes to the parameter sheet, click the "synchronize" button,
then click "Done". If no changes, click "Done" and return to the R console.
}---"
          message(msg); utils::flush.console()
          keystone::press_enter_to_continue()

          params <- DataEditR::data_edit(params)
        },

        save = {
          defaultFileName <- sprintf("volta-parameters_%s",
            keystone::make_current_timestamp(use_seconds = TRUE, seconds_sep = "+"))
          params_path <- svDialogs::dlg_save(default = defaultFileName, title = "Save volta parameters sheet",
            filters = svDialogs::dlg_filters[c("xls", "csv"), ])$res

          rio::export(params, params_path)

          msg <- paste0(
r"---{
To load this parameter sheet for a future volta analysis, use command:

}---",
            sprintf("prepare_data(\"%s\")\n", params_path))
          message(msg); utils::flush.console()
        },

        continue = {
          break
        }
      )
    }
  }

  ## Copy 'params' to the package global variable
  if (copy_params) {
    current <- get(".volta", envir = asNamespace("volta"))
    current$params <- params
    assign(".volta", current, envir = asNamespace("volta"))
  }

  msg <- paste0(
r"---{
The volta parameter sheet is ready to use. To continue the analysis, type:

}---",
    "get_voltration_data()")
  message(msg); utils::flush.console()

  ## Invisibly return 'params' as a function value, too
  invisible(params)
}

#' @export
guess_parameters <- function(
  path = ".", # If list, keep as is (check existence); vector gets processed
  ## 'analysis_channels_re' is one of: 1-element regex/character vector;
  ##   2+-element vector of channel names
  analysis_channels_re = stringr::regex("^(?!.*?(fsc|ssc|time|ir)).*(?<!-(h|w))$", ignore_case = TRUE),
  ## These calcs may rely on data stored in ad-hoc keywords
  ##   for parameter n, e.g. PnV (_V_oltage), PnG (_G_ain):
  cv_keyword_fmt = "$P%sV",
  parse_quantity = expression(pData %>% dplyr::mutate(quantity = quantity %>% readr::parse_number())),
  analysis_channels_sep = ",", quantity_sep = analysis_channels_sep,
  channels_callback = NULL,
  same_word = "SAME", spill_word = "SPILL",
  use_spillover_channels = FALSE,
  list.dirs... = list(),
  list.files... = list(),
  use_filename_voltage = TRUE,
  voltage_re = "(?<=^|[^\\dA-Za-z])[1-9][\\d]{2}(?=[^\\dA-Za-z]|$)", default_voltage = 300
)
{
  ## 'path' can be a list, or a vector mix of directory- & file paths; sort them out:
  if (!is.list(path)) { # A list from 'prepare_data()' should be ready to go
    pathAbo <- path
    path <- sapply(path, tools::file_path_as_absolute, USE.NAMES = FALSE)
    isDirectory <- sapply(path, function(a) utils::file_test("-d", a))
    singleFiles <- structure(path[!isDirectory], .Names = dirname(path[!isDirectory])) %>% as.list
    path <- path[isDirectory]

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
        filePaths <- do.call(keystone::list_files, list.filesArgs)

        filePaths
      }, simplify = FALSE) %>% purrr::compact() # Remove blank elements

    ## Make sure any individual files actually in the same directory are grouped together:
    ff <- purrr::reduce(list(ff, singleFiles), merge_named_list_elements)
  } else {
    ff <- path
    ## TODO: Check all files for existence?
  }

  if (is_invalid(ff))
    stop("No FCS files found")

  d <- plyr::ldply(names(ff),
    function(a)
    {
      machine <- basename(a) # 'machine' is required in the spreadsheet

      dd <- plyr::ldply(ff[[a]],
        function(fcs_path) # 'fcs_path' is required in the spreadsheet
        {
          fileName <- basename(fcs_path)
          ## Can we guess the voltage/gain from the file name?
          fileNameVG <- stringr::str_extract(fileName, voltage_re) # Use 'string_extract_all()' here?

          ## Now let's get info from the file itself
          fcs1 <- flowCore::read.FCS(
            fcs_path,
            truncate_max_range = FALSE,
            ## Retrieve only single event, but keep parameters, keywords, & value from header of FCS file:
            which.lines = 1,
            transformation = FALSE,
          )
          keywords <- fcs1 %>% flowCore::keyword()

          ## If we want to try & use keyword values of the quantity of interest:
          pData <- fcs1 %>% flowCore::parameters() %>% flowCore::pData() %>%
            tibble::rownames_to_column(var = "parameter") %>%
            dplyr::mutate(
              quantity_keyword = sprintf(cv_keyword_fmt, stringr::str_extract(parameter, "\\d+$")),
              desc = dplyr::case_when(is.na(desc) ~ name, TRUE ~ desc)
            ) %>%
            dplyr::left_join(
              keystone::dataframe(
                parameter = .$parameter,
                quantity =
                  sapply(.$quantity_keyword,
                    function(a) { r <- keywords[[a]]; if (is_invalid(r)) r <- NA_character_; r })
              ), by = "parameter"
            )

          pData <- keystone::poly_eval(parse_quantity)

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

          ## Now, identify all the channels to be used in the analysis
          if (is.numeric(analysis_channels_re)) {
            analysis_channels <- flowCore::colnames(fcs1)[analysis_channels_re]
          } else if (is.character(analysis_channels_re)) {
            if (length(analysis_channels_re) > 1) {
              analysis_channels <- sapply(analysis_channels_re,
                function(b) stringr::str_subset(zz, stringr::regex(b, ignore_case = TRUE)), simplify = FALSE) %>%
                unlist(use.names = FALSE) %>% unique
            } else {
              if (stringr::str_detect(analysis_channels_re, # Let's try the channels in the spillover matrix
                stringr::regex(sprintf("^\\s*%s\\s*$", spill_word), ignore_case = TRUE))) {
                ## Check for spillover matrix
                spilloverKey <- names(keywords) %>% stringr::str_subset(stringr::regex("spill", ignore_case = TRUE))
                hasSpilloverMatrix <- spilloverKey %>% `[[`(keywords, .) %>% is.matrix
                if (hasSpilloverMatrix)
                  analysis_channels <- spilloverKey %>% `[[`(keywords, .) %>% colnames
                else # If no spillover matrix, default to using all channels:
                  analysis_channels <- flowCore::colnames(fcs1)
              } else {
                analysis_channels <- stringr::str_subset(flowCore::colnames(fcs1), analysis_channels_re)
              }
            }
          } else {
            analysis_channels <- flowCore::colnames(fcs1)
          }

          analysis_channels <- intersect(analysis_channels, flowCore::colnames(fcs1))
          if (is_invalid(analysis_channels))
            analysis_channels <- flowCore::colnames(fcs1)

          ## At this point, allow an expression or function "callback" for trickier channel extractions,
          ##   using e.g. 'fileName', 'keywords', 'pData', 'analysis_channels'
          keystone::poly_eval(channels_callback)

          ## Now back to finding the voltage/gain quantities; start w/ info in 'pData':
          pDataVG <- pData %>% dplyr::filter(name %in% analysis_channels) %>% dplyr::pull(quantity)
          if (length(unique(pDataVG)) > 1)
            probable_vg <- paste(pDataVG, collapse = quantity_sep)
          else
            probable_vg <- unique(pDataVG)

          ## We should favor 'pDataVG', but if it doesn't work, here are some Hail Marys ...
          if (use_filename_voltage || is_invalid(probable_vg)) {
            ## If most of the channels are tested, then several keywords should contain the same voltage/gain:
            keywordsVG <- sapply(keywords,
              function(b) { if (!inherits(b, "character")) return (NA); stringr::str_extract(b, "^\\s*[1-9][\\d]0\\s*$") },
              simplify = TRUE)
            ## Return the most frequent voltage/gain-like number(s) among the keywords
            probable_vg <- keystone::stat_mode(c(fileNameVG, keywordsVG) %>% readr::parse_number(), na.rm = TRUE)

            ## Reluctantly, I'm going to assume that the file name is our best source:
            if (use_filename_voltage && !is_invalid(fileNameVG))
              probable_vg <- fileNameVG

            ## These are sort of last resorts:
            if (is_invalid(probable_vg)) {
              probable_vg <- default_voltage
              warning(sprintf("For machine %s, file %s: no voltage/gain value found; default of %s used",
                machine, fileName, default_voltage))
            }
            if (length(probable_vg) > 1)
              probable_vg <- sample(probable_vg, 1) # 'probable_vg' is required in the spreadsheet
          }

          analysis_channels <- # 'analysis_channels' is required in the spreadsheet
            paste(analysis_channels, collapse = analysis_channels_sep)
          ## To reverse: stringr::str_split_1(analysis_channels, analysis_channels_sep)

          keystone::dataframe(
            machine = machine,
            quantity = probable_vg %>% as.character,
            channels = analysis_channels,
            path = fcs_path) %>% dplyr::arrange(quantity)
        })

      ii <- plyr::alply(dd, 1, function(a) { stringr::str_split_1(a$channels, analysis_channels_sep) %>%
        sort %>% paste(collapse = analysis_channels_sep) }) %>% unlist(use.names = FALSE) %>%
        as.factor %>% as.numeric %>% vctrs::vec_duplicate_id() %>% unique
      same <- rep(same_word, NROW(dd)); same[ii] <- dd$channels[ii]
      dd$channels <- same

      dd
    }
  )

  ## Return data frame defining parameters for entire volta analysis:
  d
}
