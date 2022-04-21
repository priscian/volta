## V. https://www.thermofisher.com/us/en/home/references/newsletters-and-journals/bioprobes-journal-of-cell-biology-applications/bioprobes-78/photomultiplier-tube-pmt-optimization-attune-nxt-flow-cytometer.html

#' @export
staining_index <- function(
  x,
  stain_groups,
  na.rm = FALSE
)
{
  negative <- stain_groups$unstained %>% unclass
  positive <- stain_groups$stained %>% unclass

  (stats::median(positive, na.rm = na.rm) - stats::median(negative, na.rm = na.rm)) /
    (2 * stats::sd(negative, na.rm = na.rm))
}


#' @export
alt_staining_index <- function(
  x,
  stain_groups,
  na.rm = FALSE
)
{
  negative <- stain_groups$unstained %>% unclass
  positive <- stain_groups$stained %>% unclass

  stats::median(positive, na.rm = na.rm) / stats::sd(negative, na.rm = na.rm)
}


#' @export
voltration_index <- function(
  x,
  stain_groups,
  quantity,
  na.rm = FALSE
)
{

  ## Staining Index (SI)
  #si <- staining_index(NULL, stain_groups)

  ## Alternative Staining Index (Alt SI)
  alt_si <- staining_index(NULL, stain_groups)

  vi <- alt_si / sqrt(quantity)

  vi
}
