# volta
Flow cytometers: calculate minimum per-channel voltage for optimal resolution sensitivity to dim vs. unstained events.

## Preliminaries
The *volta* R package is fairly easy to set up. In an R session:
```r
install.packages("devtools") # If necessary.
# https://github.com/r-lib/remotes#environment-variables
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
devtools::install_github("priscian/volta")
library(volta)

## Once the package has been installed as described above, all you need to use it is:
library(volta)
```

## Trying out *volta*
```
## Analyze some test data:
test_volta()
```

This will produce a plot and a data frame containing the minimum voltage for each channel.

![Some major monthly global average temperature time series.](<inst/images/001 - Maecker & Trotter 2006.png>)

Here's the test code in detail:
```r
### Maecker & Trotter 2006

library(volta)
e <- new.env()
load(system.file("extdata/peak2-maecker+trotter2006_2019-09-19+54111.RData",
  package = "volta"), envir = e)
attr(e$d, "wavelength") <- c(520, 785, 578, 695) # This attribute is necessary!
## N.B. Not necessary to set this attribute:
attr(e$d, "experiment_name") <- "Maecker & Trotter 2006"
v2006 <-
  plot_peak2_data(
    x = list(e$d),
    image_dir = system.file("images", package = "volta"),
    plot_series... = list(x_var = "voltage"),
    create_smooth_variables... = list(x_var = "voltage"),
    save_png = FALSE
  )

print(v2006[[1]]$table)
```

### More information
*volta* is presented here as a working beta. For more information on what the package offers, check out
```r
library(help = volta)
```
from the R command line.

## Example from scratch

*volta*'s installation directory includes Peak 2 FCS files from the URMC (Rochester, New York, USA) Flow Core's suite of instruments; featured are: Animal (BD LSR II), Dr. Teeth (BD LSRFortessa™), and Statler (BD FACSAria™ II)—guess the naming theme! Runs of Peak 2 beads at various voltages are stored in directories named after each machine. Let's assume that we have enough data to determine the minimum per-channel voltage for optimal resolution sensitivity for each instrument, and let's furthermore put *volta* into action finding those voltages for us.

```r
library(volta)

```
