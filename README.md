# volta
Flow cytometers: calculate minimum per-channel voltage for optimal resolution sensitivity to dim vs. unstained events.

## Preliminaries
The *volta* R package is fairly easy to set up. In an R session:
```r
install.packages("remotes") # If necessary.
# https://github.com/r-lib/remotes#environment-variables
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
remotes::install_github("priscian/volta", INSTALL_opts = c("--no-multiarch"))

## Once the package has been installed as described above, all you need to use it is:
options(keystone_parallel = TRUE); library(volta)
```

## Trying out *volta*
```
## Analyze some test data:
test_volta()
```

This will produce a plot and a data frame containing the minimum voltage for each channel.

![Reproduction of Figure 1B in Maecker & Trotter 2006 with volta-calculated inflection points.](<inst/images/001 - Maecker & Trotter 2006.png>)

Here's the test code in detail:
```r
### Maecker & Trotter 2006

options(keystone_parallel = TRUE)
library(volta); library(magrittr)
e <- new.env()
load(system.file("extdata/peak2-maecker+trotter2006_2019-09-19+54111.RData",
  package = "volta"), envir = e)
p2 <- list(`Maecker & Trotter 2006` = structure(e$d %>% dplyr::rename(quantity = 1),
  color = keystone::wavelength2col(c(520, 785, 578, 695)), # This attribute is necessary!
  ## N.B. Not necessary to set this attribute:
  experiment_name = "Maecker & Trotter 2006"))
fip <- p2 %>% volta:::find_inflection_point()
p2 <- sapply(names(p2),
  function(i) { (structure(p2[[i]], plot_data = fip[[i]]) %>%
    keystone::add_class("volta")) },
  simplify = FALSE)

plot(p2[[1]])

attr(p2[[1]], "plot_data")$inflection_points %>%
  dplyr::rename(
    !!names(formals(volta:::plot.volta)$x_var_lab %>% keystone::poly_eval()) := "inflection",
    !!names(formals(volta:::plot.volta)$y_var_lab %>% keystone::poly_eval()) := "y"
  ) %>% print
```

### More information
*volta* is presented here as an early-version release. For more information on what the package offers, check out
```r
library(help = volta)
```
from the R command line.

## Example from scratch

*volta*'s installation directory includes Peak 2 FCS files from the URMC (Rochester, New York, USA) Flow Core's suite of instruments; featured are: Animal (BD LSR II), Dr. Teeth (BD LSRFortessa™), and Statler (BD FACSAria™ II)—guess the naming theme! Runs of Peak 2 beads at various voltages are stored in directories named after each machine. Let's assume that we have enough data to determine the minimum per-channel voltage for optimal resolution sensitivity for each instrument, and let's furthermore put *volta* into action finding those voltages for us.

```r
options(keystone_parallel = TRUE) # Allow parallel processing via the "future" package
options(volta_interactive_off = TRUE) # Turn off the interactive R console "wizard"
library(volta)
## Create parameters sheet for voltration analysis from FCS files
p <- prepare_data(system.file("extdata/peak2-fcs", package = "volta"))
## Calculate voltration data set with inflection points
v <- get_voltration_data(p)
## Plot CV vs. PMT voltage curves calculated from Peak 2 data
plot_voltration_data(v)
```

### Graphical output
<!--<br/>-->
[//]: <br/>

![Results of volta run on URMC Peak 2 data: Animal.](<inst/images/002 - Animal.png>)

![Results of volta run on URMC Peak 2 data: Dr. Teeth.](<inst/images/003 - Dr Teeth.png>)

![Results of volta run on URMC Peak 2 data: Statler.](<inst/images/004 - Statler.png>)

### Detail of return value from function `get_voltration_data()`
<!--<br/>-->
[//]: <br/>

```r
## Show table of CV-vs-V inflection points for each channel & each instrument
> sapply(v, function(a) attr(a, "plot_data")$inflection_points %>%
+   dplyr::rename(PMT_voltage = 2, log10_CV = 3), simplify = FALSE)
$`Animal Pk2 122221`
             channel PMT_voltage  log10_CV
1    Blue B 515/20-A         370 1.7023795
2    Blue A 710/50-A         520 1.7550937
3  Violet H 450/50-A         330 0.7260804
4  Violet G 550/40-A         450 1.3776819
5  Violet F 560/40-A         470 1.3634206
6  Violet E 585/42-A         440 1.3478023
7  Violet D 605/40-A         440 1.3997432
8  Violet C 660/40-A         490 1.3815358
9  Violet B 705/70-A         460 1.6562864
10 Violet A 780/60-A         540 2.0365942
11    Red C 660/20-A         430 1.2968331
12    Red B 710/50-A         420 1.5210990
13    Red A 780/60-A         460 1.7625997
14  Green E 575/25-A         380 1.0859096
15  Green D 610/20-A         420 1.1867483
16  Green C 660/40-A         460 1.4167239
17  Green B 710/50-A         500 1.4764845
18  Green A 780/40-A         540 1.8206604

$`Dr Teeth Pk 2 122221`
             channel PMT_voltage  log10_CV
1    Blue B 530/30-A         450 1.3669519
2    Blue A 710/50-A         520 1.7164870
3     Red C 670/14-A         450 1.2678660
4     Red B 730/45-A         450 1.3876899
5     Red A 780/60-A         470 1.5906573
6  Violet F 431/28-A         420 0.6520710
7  Violet E 525/50-A         430 0.8359811
8  Violet D 610/20-A         500 1.5111939
9  Violet C 660/20-A         510 1.5843604
10 Violet B 710/50-A         560 1.6994997
11 Violet A 780/60-A         570 2.0850087
12     YG E 586/15-A         430 1.3781064
13     YG D 610/20-A         440 1.4030870
14     YG C 670/30-A         430 1.5033048
15     YG B 710/50-A         470 1.7330610
16     YG A 780/60-A         540 1.7998980
17     UV B 379/28-A         460 1.3941904
18     UV A 740/35-A         600 2.3633671

$`Peak 2 statler LN Z01 122221`
             channel PMT_voltage  log10_CV
1     Red C 670/14-A         420 1.3536483
2     Red B 730/45-A         430 1.4581357
3     Red A 780/60-A         480 1.7466976
4    Blue B 510/21-A         510 1.7370133
5    Blue A 710/50-A         390 2.1336790
6  Violet G 450/50-A         490 0.9208474
7  Violet F 550/50-A         420 1.4459166
8  Violet E 585/42-A         410 1.5869469
9  Violet D 605/40-A         430 1.6352509
10 Violet C 660/20-A         510 2.0387778
11 Violet B 705/70-A         440 2.0057206
12 Violet A 780/60-A         520 2.4778029
13  Green E 575/25-A         460 1.3724170
14  Green D 610/20-A         490 1.4257077
15  Green C 660/20-A         510 1.7861783
16  Green B 710/50-A         520 1.8279439
17  Green A 780/60-A         590 2.2966084
```

## Interactive session
<!--<br/>-->
[//]: <br/>

If `options(volta_interactive_off = FALSE)` or it isn't set, *volta* runs in interactive mode on the R console. (Use `options(volta_interactive_off = TRUE)` before starting a session to turn off *volta*'s interactive mode.) Interactive mode will try to guide the user through a complete analysis in as simple a manner as possible.

![Complete interactive session of volta on R console.](<inst/images/volta-1.0.0-interactive-session-cropped.png>)

## Function signatures
<!--<br/>-->
[//]: <br/>

This is the function signature of `prepare_data()` with annotations:

```r
prepare_data <- function(
  x,
  ..., # Passed on to 'guess_parameters()'
  choose_excel = FALSE,
  copy_params = TRUE
)
{...}
```

This is the function signature of `get_voltration_data()` with annotations:

```r
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
{...}
```

This is the function signature of `plot_voltration_data()` with annotations:

```r
plot_voltration_data <- function(
  x, # Result from 'get_voltration_data()'
  report_path = NULL,
  image_dir = NULL,
  save_png = FALSE, # If TRUE, save PNG plots to report directory
  png... = list(),
  x_var_lab = c(PMT_voltage = "PMT Voltage"),
  ## N.B. 'names(y_var_lab)' will correctly extract "log10_CV":
  y_var_lab = c(log10_CV = expression(paste(log[10], " CV"))),
  plot_series... = list(),
  points... = list(),
  xlsx_expression = NULL,
  ## Noli me tangere; when TRUE, used for IRR checks:
  plot_individual_channels = FALSE
)
{...}
```
