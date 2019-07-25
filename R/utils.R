#' @export
char_sort <- function(x, s)
{
  x[which(x %in% s)] <- s[which(s %in% x)]

  x
}


#' @export
only_selected_series <- function(x, series, common_columns, sort = FALSE, range = NULL, year_var = NULL, ...)
{
  if (missing(series))
    series <- NULL

  colNames <- c(intersect(colnames(x), c(common_columns, series)))
  if (!sort)
    colNames <- char_sort(colNames, series)
  r <- x[, colNames, ...]

  if (!is.null(range)) {
    if (!is.null(year_var))
      yearVar <- year_var
    else if ("yr_part" %in% colNames)
      yearVar <- "yr_part"
    else if ("year" %in% colNames)
      yearVar <- "year"

    r <- r[r[[yearVar]] >= ifelse(is.na(range[1]), min(r[[yearVar]], na.rm = TRUE), range[1]) & r[[yearVar]] <= ifelse(is.na(range[2]), max(r[[yearVar]]), range[2]), ]
  }

  r
}

#' @export
oss <- only_selected_series


#' @export
view_only_selected_series <- function(..., fun = View)
{
  fun(only_selected_series(...))
}

#' @export
vss <- view_only_selected_series
