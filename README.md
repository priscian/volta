# volta
Flow cytometers: calculate minimum per-channel voltage for optimal resolution sensitivity to dim vs. unstained events.

## Preliminaries
The *volta* R package is fairly easy to set up. In an R session:
```
install.packages("devtools") # If necessary.
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true") # https://github.com/r-lib/remotes#environment-variables
devtools::install_github("priscian/volta")
library(volta)

## Once the package has been installed as described above, all you need to use it is:
library(volta)
```

## Using *volta*
```
## Analyze some test data:
test_volta()
```

This will produce a plot and a data frame containing the minimum voltage for each channel.

![Some major monthly global average temperature time series.](<inst/images/001 - Maecker & Trotter 2006.png>)

Here's the test code in detail:
```
### Maecker & Trotter 2006

e <- new.env()
load(system.file("extdata/peak2-maecker+trotter2006_2019-09-19+54111.RData",
  package = "volta"), envir = e)
attr(e$d, "wavelength") <- c(520, 785, 578, 695) # This attribute is necessary!
attr(e$d, "experiment_name") <- "Maecker & Trotter 2006" # Not necessary to set this attribute.
v2006 <- plot_peak2_data(x = list(e$d), image_dir = system.file("images", package = "volta"),
  save_png = FALSE)
print(v2006)
```

### More information
*volta* is presented here as a working beta. For more information on what the package offers, check out
```
library(help = volta)
```
from the R command line.
