#' @export
plot_series <- function(
  x, # This should be a data frame.
  series,
  year_var,
  start = NULL, end = NULL,
  ma = NULL, ma_sides = 1L,
  plot_type = c("single", "multiple"),
  type = "l",
  main = "", xlab = "", ylab = "",
  unit = NULL,
  col = NULL, col_fun = colorspace::rainbow_hcl, col_fun... = list(l = 65), alpha = 0.5, lwd = 2,
  legend... = list(),
  trend = FALSE, trend_lwd = lwd, trend_legend_inset = c(0.2, 0.2),
  add = FALSE,
  segmented = FALSE, segmented... = list(),
  loess = FALSE, loess... = list(),
  plot.segmented... = list(),
  mark_segments = FALSE, vline... = list(),
  start_callback = NULL, end_callback = NULL,
  save_png = FALSE, save_png_dir = ".",
  png... = list(),
  ...
)
{
  plot_type <- match.arg(plot_type)

  ## This is to avoid an roxygen error described here: https://github.com/klutometis/roxygen/issues/592
  if (is.null(unit))
    unit <- "\u00b0C"

  y <- zoo(x[, c(year_var, series)], order.by = x[[year_var]])
  y <- subset(y, na_unwrap(as.matrix(y[, series]))) # Remove trailing NAs.
  w <- interpNA(y, "linear", unwrap = TRUE)

  ## Create moving-average variables if requested.
  w[, series] <- MA(w[, series], ma, sides = ma_sides); w <- zoo(w, order.by = index(y))
  maText <- ""
  if (!is.null(ma))
    maText <- "(" %_% ma %_% "-month moving average)"

  y <- window(y, start = start, end = end, extend = TRUE) # N.B. This can't be 'extend()'ed.
  w <- window(w, start = start, end = end, extend = TRUE)
  #wz <- as.zoo(w) # Unnecessary, but for legacy's sake.
  wz <- w

  ## Series colors.
  if (is.null(col)) {
    col <- seq_along(series)
    col_funArgs <- list(
      n = length(col)
    )
    col_funArgs <- modifyList(col_funArgs, col_fun...)
    col <- suppressWarnings(do.call(col_fun, col_funArgs))
    ## Some other possible function calls:
    #col <- rainbow(length(col))
    #col <- terrain.colors(length(col))
    #col <- topo.colors(length(col))
    #col <- suppressWarnings(brewer.pal(length(col), "Spectral")) # Or "Paired".
    #col <- matlab.like2(length(col)) # From package "colorRamps".
  }
  col <- rep(col, length.out = length(series))
  col <- alpha(col, alpha)
  names(col) <- series

  ## Arguments for saving the plot.
  imageDir <- save_png_dir

  if (save_png) {
    pngArgs <- list(
      filename = paste(imageDir, filename, sep = "/"),
      width = 12.5,
      height = 7.3,
      units = "in",
      res = 600
    )
    pngArgs <- modifyList(pngArgs, png...)
    do.call("png", pngArgs)
  }

  ## Arguments for plotting.
  xaxt <- "n"
  if (dev.cur() == 1L) # If a graphics device is active, plot there instead of opening a new device.
    dev.new(width = 12.5, height = 7.3) # New default device of 1200 × 700 px at 96 DPI.
    #dev.new(width = 7.3, height = 7.3) # New default device of 700 × 700 px at 96 DPI.

  if (!add)
    plot(w[, series], plot.type = plot_type, type = "n", xaxs = "r", xaxt = "s", xlab = xlab, ylab = ylab, main = main, ...)
  if (maText != "") mtext(maText, 3L)

  if (!add) {
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"))
  }

  ## Evaluate expression after creating graphics device but before plotting series.
  if (!is.null(start_callback))
    eval(start_callback)

  par(new = TRUE)
  if (add) {
    plotSeries <- series
    for (i in seq_along(plotSeries))
      lines(wz[, plotSeries[i]], type = type, col = col[i], lwd = lwd, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ...) # I.e. 'plot.zoo()'.
  }
  else
    plot(wz[, series], screens = 1L, plot.type = plot_type, type = type, col = col, lwd = lwd, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ...) # I.e. 'plot.zoo()'.

  legendArgs <- list(
    x = "topleft",
    legend = series %_% ifelse(loess, " (+ LOESS)", ""),
    col = col,
    lwd = lwd,
    bty = "n",
    cex = 0.8
  )
  legendArgs <- modifyList(legendArgs, legend...)
  do.call("legend", legendArgs)

  ## LOESS smooth.
  if (loess) {
    for (s in series) {
      loessArgs = list(
        formula = eval(substitute(s ~ x, list(s = as.name(s), x = as.name(year_var)))),
        data = y[, c(year_var, series)],
        span = 0.2
      )
      loessArgs <- modifyList(loessArgs, loess...)

      l <- do.call("loess", loessArgs)

      lwd <- 2
      if (!is.null(loessArgs$lwd))
        lwd <- loessArgs$lwd
      lines(l$x, l$fit, col = col[s], lwd = lwd)
    }
  }

  r <- list()

  ## Linear trends.
  if (trend) {
    m <- list()
    m$series <- series
    m$range <- list(start = start(wz), end = end(wz))
    m$col <- col
    m$data <- y
    for (s in m$series) {
      m[[s]]$lm <- lm(eval(substitute(b ~ x, list(b = as.name(s), x = as.name(year_var)))), data = m$data, x = TRUE)
      m[[s]]$warming <- coef(m[[s]]$lm)[2] * diff(range(m[[s]]$lm$model[, 2]))
      m[[s]]$rate <- coef(m[[s]]$lm)[2] * 10
      m[[s]]$rateText <- eval(substitute(expression(paste(Delta, " = ", r, phantom(l), unit, "/dec.", sep = "")), list(r = sprintf(m[[s]]$rate, fmt = "%+1.3f"), unit = unit)))
      m[[s]]$col <- m$col[s]
    }

    legendText <- NULL
    for (s in m$series) {
      ## Set clipping region for 'abline()'.
      xRange <- range(wz[!is.na(wz[, s]), year_var], na.rm = TRUE)
      yRange <- range(wz[, s], na.rm = TRUE)
      usr <- par("usr")
      clip(xRange[1], xRange[2], yRange[1], yRange[2])

      abline(m[[s]]$lm, col = m[[s]]$col, lwd = trend_lwd, untf = TRUE)

      ## Reset clipping to plot region.
      do.call("clip", as.list(usr))

      legendText <- c(legendText, m[[s]]$rateText)
    }

    if (!is.null(trend_legend_inset))
      #legend("bottomright", inset = trend_legend_inset, legend = legendText, col = m$col, lwd = trend_lwd, bty = "n", cex = 0.8)
      legend("topright", legend = legendText, col = m$col, lwd = trend_lwd, bty = "n", cex = 0.8)

    r$trend <- m
  }

  ## Segmented linear regression.
  if (segmented) {
    segmentedArgs <- list(
      x = x,
      series = series,
      col = col,
      start = start,
      end = end,
      year_var = year_var
    )
    segmentedArgs <- modifyList(segmentedArgs, segmented...)
    sm <- do.call("fit_segmented_model", segmentedArgs)

    for (i in names(sm$piecewise)) {
      ## Set clipping region for 'plot.segmented()' and 'abline()'.
      xRange <- range(wz[!is.na(wz[, i]), year_var], na.rm = TRUE)
      yRange <- range(wz[, i], na.rm = TRUE)
      usr <- par("usr")
      clip(xRange[1], xRange[2], yRange[1], yRange[2])

      x <- sm$piecewise[[i]]$sm

      if (!is.null(x)) {
        plot.segmentedArgs <- list(
          x = x,
          add = TRUE,
          rug = FALSE,
          lwd = 2,
          #lty = "longdash",
          col = col[i],
          alpha = alpha
        )
        plot.segmentedArgs <- modifyList(plot.segmentedArgs, plot.segmented...)
        dev_null <- do.call("plot", plot.segmentedArgs)

        if (mark_segments) {
          vlineArgs <- list(
            mark_years = sprintf(sm$piecewise[[i]]$sm$psi[, 2], fmt = "%1.1f")
          )
          vlineArgs <- modifyList(vlineArgs, vline...)
          do.call("vline", vlineArgs)
        }
      } else {
        lwd <- ifelse(is.null(plot.segmented...$lwd), 2, plot.segmented...$lwd)
        abline(sm$piecewise[[i]]$lm, col = col[i], lwd = lwd, untf = TRUE)
      }

      ## Reset clipping to plot region.
      do.call("clip", as.list(usr))
    }

    r$segmented <- sm
  }

  if (!is.null(end_callback))
    eval(end_callback)

  if (length(r) > 0L)
    return (invisible(r))

  return (nop())
}


color_nm_map <- c(
  red = 700,
  orange = 620,
  yellow = 580,
  green = 530,
  cyan = 500,
  blue = 470,
  indigo = 450,
  violet = 420
)

## V. JavaScript source for https://academo.org/demos/wavelength-to-colour-relationship/
nm_to_rgb <- function(wavelength, Gamma = 0.8, IntensityMax = 255)
{
  if ((wavelength >= 380) && (wavelength < 440)) {
    red <- -(wavelength - 440) / (440 - 380)
    green <- 0.0
    blue <- 1.0
  } else if ((wavelength >= 440) && (wavelength < 490)) {
    red <- 0.0
    green <- (wavelength - 440) / (490 - 440)
    blue <- 1.0
  } else if ((wavelength >= 490) && (wavelength < 510)) {
    red <- 0.0
    green <- 1.0
    blue <- -(wavelength - 510) / (510 - 490)
  } else if ((wavelength >= 510) && (wavelength < 580)) {
    red <- (wavelength - 510) / (580 - 510)
    green <- 1.0
    blue <- 0.0
  } else if ((wavelength >= 580) && (wavelength < 645)) {
    red <- 1.0
    green <- -(wavelength - 645) / (645 - 580)
    blue <- 0.0
  } else if ((wavelength >= 645) && (wavelength < 781)) {
    red <- 1.0
    green <- 0.0
    blue <- 0.0
  } else {
    red <- 0.0
    green <- 0.0
    blue <- 0.0
  }

  ## Let intensity fall off near vision limits.
  if((wavelength >= 380) && (wavelength < 420)) {
    factor <- 0.3 + 0.7 * (wavelength - 380) / (420 - 380)
  } else if((wavelength >= 420) && (wavelength < 701)) {
    factor <- 1.0
  } else if((wavelength >= 701) && (wavelength < 781)) {
    factor <- 0.3 + 0.7 * (780 - wavelength) / (780 - 700)
  } else {
    factor <- 0.0
  }

  if (red != 0) {
    red <- round(IntensityMax * (red * factor)^Gamma)
  }
  if (green != 0) {
    green <- round(IntensityMax * (green * factor)^Gamma)
  }
  if (blue != 0) {
    blue <- round(IntensityMax * (blue * factor)^Gamma)
  }

  list(R = red, G = green, B = blue)
}


#' @export
wavelength2col <- Vectorize(function(wavelength, Gamma = 0.8, IntensityMax = 255, ...)
{
  #RGB <- colorscience::heuristic.wlnm2RGB(wavelength, Gamma, IntensityMax)
  RGB <- nm_to_rgb(wavelength, Gamma, IntensityMax)

  grDevices::rgb(RGB$R, RGB$G, RGB$B, maxColorValue = IntensityMax, ...) # Can add e.g. 'alpha' here.
})
