.onLoad <- function(...)
{
  #assign("%_%", keystone::`%_%`, envir = .GlobalEnv)
  #.reload_all("keystone", redocument = FALSE)
}


.onAttach <- function(...)
{
  keystone:::.onAttach()
}


.onDetach <- function(...)
{
  keystone:::.onDetach()
}
