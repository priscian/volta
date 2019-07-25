## Fit segmented linear models to selected data.
#' @export
fit_segmented_model <- function(x,
  series,
  col = colorspace::rainbow_hcl(length(series)),
  start = NULL, end = NULL,
  yearly = FALSE, year_var = NULL,
  breakpoints... = list(),
  segmented... = list(), seg.control... = list(seed = 100),
  ...
)
{
  r <- list(data = x, series = series)
  r$range <- list(start = start, end = end)
  r$col <- col
  length(r$col) <- length(r$series); names(r$col) <- r$series

  if (!yearly) {
    g <- r$data
  }
  else { # This won't work for non-climeseries data.
    g <- as.data.frame(make_yearly_data(r$data))
    if (!is.null(start)) start <- trunc(start)
    if (!is.null(end)) end <- trunc(end)
  }

  yearVar <- ifelse(yearly, "year", "yr_part")
  if (!is.null(year_var))
    yearVar <- year_var

  r$piecewise <- list()
  for (i in r$series) {
    r$piecewise[[i]] <- list()
    r$piecewise[[i]]$col <- r$piecewise$col[i]

    h <- oss(g, i, common_columns = year_var)[na_unwrap(g[[i]]), , drop=FALSE]
    h <- h[h[[yearVar]] >= ifelse(!is.null(start), start, -Inf) & h[[yearVar]] <= ifelse(!is.null(end), end, Inf), ]

    breakpointsArgs <- list(
      formula = eval(substitute(Y ~ X, list(X = as.name(yearVar), Y = as.name(i)))),
      data = h,
      breaks = NULL
    )
    breakpointsArgs <- modifyList(breakpointsArgs, breakpoints...)
    r$piecewise[[i]]$bp <- do.call("breakpoints", breakpointsArgs)

    r$piecewise[[i]]$breaks <- r$piecewise[[i]]$bp$X[, yearVar][r$piecewise[[i]]$bp$breakpoint]

    seg.controlArgs <- list(
      stop.if.error = TRUE,
      K = length(r$piecewise[[i]]$breaks),
      n.boot = 250,
      random = FALSE,
      h = 0.3
    )
    seg.controlArgs <- modifyList(seg.controlArgs, seg.control...)
    segControl <- do.call("seg.control", seg.controlArgs)

    r$piecewise[[i]]$lm <- lm(breakpointsArgs$formula, data = h, x = TRUE, y = TRUE)

    segmentedArgs <- list(
      obj = r$piecewise[[i]]$lm,
      seg.Z = as.formula(paste("~", yearVar)),
      psi = r$piecewise[[i]]$breaks,
      control = segControl
    )
    segmentedArgs <- modifyList(segmentedArgs, segmented...)

    run_segmented <- function()
    {
      mf <- model.frame(r$piecewise[[i]]$lm)

      while (TRUE) {
        withRestarts({
          sm <- do.call("segmented", segmentedArgs)
          break
        },
          restart = function() {
            ## Which breakpoint is closest to the start or end of the time series?
            if (length(segmentedArgs$psi) > 1L) {
              segmentedArgs$psi <<- segmentedArgs$psi[-which.min(pmin(segmentedArgs$psi, NROW(mf) - segmentedArgs$psi + 1))]
            }
          })
      }

      sm
    }

    tryCatch({
      withCallingHandlers({
          sm <- run_segmented()
        },
          error = function(e) {
            message("Error: ", e$message); flush.console()
            if (any(grepl("one coef is NA: breakpoint(s) at the boundary", e$message, fixed = TRUE)) && length(segmentedArgs$psi) > 1L)
              invokeRestart("restart")
          }
      )

      r$piecewise[[i]]$sm <- sm
    }, error = function(e) { message("Warning: No breakpoint(s) found"); flush.console() })
  }

  r
}
