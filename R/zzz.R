.onLoad <- function(...)
{
  #assign("%_%", keystone::`%_%`, envir = .GlobalEnv)
  #.reload_all("keystone", redocument = FALSE)
}


.onAttach <- function(...)
{
  keystone:::.onAttach()

  ## Unlock '.[package]' variable to allow its modification:
  unlockBinding(".volta", asNamespace("volta"))

  ## Startup message
  msg <- voltaStartupMessage()
  if (!interactive())
    msg[1] <- paste("Package 'volta' version", packageVersion("volta"))

  packageStartupMessage(msg)

  invisible()
}


.onDetach <- function(...)
{
  keystone:::.onDetach()
}


voltaStartupMessage <- function()
{
  msg <- c(paste0(
r"---{__       ______________________      ___
\ \     / / __  / /___  ___/ _ \   _/  /
 \ \   / / / / / /   / /  / /_\ \ /  _/
  \ \_/ / /_/ / /___/ /  / _____ \  /
   \___/_____/_____/_/  /_/    /\_\/
                             _/ /}---",
sprintf("\nVersion %-12s        /__/\n", packageVersion("volta")),
r"---{https://CART.urmc.edu      //
                          /'
Type 'citation("volta")' to acknowledge this R package in publications.

Use 'prepare_data()' or 'prepare_data(choose_excel = TRUE)' to get started.
}---")
  )

  return(msg)
}
