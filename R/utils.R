## https://stackoverflow.com/questions/9519543/merge-two-lists-in-r/9519964#9519964
merge_named_list_elements <- function(x, val)
{
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
      merge_named_list_elements(x[[v]], val[[v]])
    else c(x[[v]], val[[v]])
  }

  x
}
