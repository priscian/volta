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

![Reproduction of Figure 1B in Maecker & Trotter 2006 with volta-calculated inflection points.](<inst/images/001 - Maecker & Trotter 2006.png>)

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
## Retrieve Peak 2 data
x <- get_peak2_data(
  system.file("extdata/peak2-fcs", package = "volta"),
  list.dirs... = list(recursive = TRUE), use_spillover_channels = TRUE)
## Plot CV vs. PMT voltage curves calculated from Peak 2 data
p <- plot_peak2_data(x, max_cv_threshold = 2000, keep_longest_continuous = TRUE,
  round_to_nearest_volts = 10)
```

### Graphical output
<!--<br/>-->
[//]: <br/>

![Results of volta run on URMC Peak 2 data: Animal.](<inst/images/002 - Animal.png>)

![Results of volta run on URMC Peak 2 data: Dr. Teeth.](<inst/images/003 - Dr Teeth.png>)

![Results of volta run on URMC Peak 2 data: Statler.](<inst/images/004 - Statler.png>)

### Function return value
<!--<br/>-->
[//]: <br/>

```r
## Show table of CV-vs-V inflection points for each channel & each instrument
> sapply(p, `[[`, i = "table", simplify = FALSE) %>% print
$`[...]/Animal Pk2 122221`
             channel PMT_voltage log10_CV
1    Blue B 515/20-A         370     1.70
2    Blue A 710/50-A         540     1.67
3  Violet H 450/50-A         330     0.73
4  Violet G 550/40-A         450     1.38
5  Violet F 560/40-A         470     1.36
6  Violet E 585/42-A         440     1.35
7  Violet D 605/40-A         440     1.40
8  Violet C 660/40-A         490     1.38
9  Violet B 705/70-A         460     1.66
10 Violet A 780/60-A         540     2.04
11    Red C 660/20-A         430     1.30
12    Red B 710/50-A         420     1.52
13    Red A 780/60-A         460     1.76
14  Green E 575/25-A         380     1.09
15  Green D 610/20-A         420     1.19
16  Green C 660/40-A         460     1.42
17  Green B 710/50-A         500     1.48
18  Green A 780/40-A         540     1.82

$`[...]/Dr Teeth Pk 2 122221`
             channel PMT_voltage log10_CV
1    Blue B 530/30-A         450     1.37
2    Blue A 710/50-A         520     1.72
3     Red C 670/14-A         450     1.27
4     Red B 730/45-A         450     1.39
5     Red A 780/60-A         470     1.59
6  Violet F 431/28-A         420     0.65
7  Violet E 525/50-A         430     0.84
8  Violet D 610/20-A         500     1.51
9  Violet C 660/20-A         510     1.58
10 Violet B 710/50-A         560     1.70
11 Violet A 780/60-A         570     2.09
12     YG E 586/15-A         430     1.38
13     YG D 610/20-A         440     1.40
14     YG C 670/30-A         430     1.50
15     YG B 710/50-A         470     1.73
16     YG A 780/60-A         540     1.80
17     UV B 379/28-A         480     1.33
18     UV A 740/35-A         630     2.32

$`[...]/Peak 2 statler LN Z01 122221`
             channel PMT_voltage log10_CV
1     Red C 670/14-A         420     1.35
2     Red B 730/45-A         430     1.46
3     Red A 780/60-A         480     1.75
4    Blue B 510/21-A         510     1.74
5    Blue A 710/50-A         390     2.13
6  Violet G 450/50-A         490     0.92
7  Violet F 550/50-A         420     1.45
8  Violet E 585/42-A         410     1.59
9  Violet D 605/40-A         430     1.64
10 Violet C 660/20-A         510     2.04
11 Violet B 705/70-A         440     2.01
12 Violet A 780/60-A         520     2.48
13  Green E 575/25-A         460     1.37
14  Green D 610/20-A         490     1.43
15  Green C 660/20-A         510     1.79
16  Green B 710/50-A         520     1.83
17  Green A 780/60-A         590     2.30
```

### Function signatures
<!--<br/>-->
[//]: <br/>

This is the function signature of `get_peak2_data()` with annotations:

```r
get_peak2_data <- function(
  path = ".",
  ## These calcs typically rely on data stored in ad-hoc keywords
  ##   for parameter n, e.g. PnV (_V_oltage), PnG (_G_ain):
  cv_keyword_fmt = "$P%sV",
  ## 'analysis_channels_re' is one of: 1-element regex/character vector;
  ##   2+-element vector of channel names
  analysis_channels_re = stringr::regex("^(?!(fsc|ssc|time).*)",
    ignore_case = TRUE),
  ## 'use_spillover_channels' overrides 'analysis_channels_re' if TRUE;
  ##   if TRUE, use channel names from spillover matrix
  use_spillover_channels = FALSE,
  create_si_data = FALSE,
  prepare_augmented_fcs_data... = list(),
  list.dirs... = list(),
  list.files... = list(),
  read.FCS... = list(),
  min_zero = FALSE,
  ## Default CV is BD Biosciences' robust coefficient of variation (BD-rCV):
  cv_fun = keystone::bd_rcv,
  cv_multiplier = 100,
  replace_with_na_condition = ~ is.infinite(.x) | is.nan(.x) | .x < 0.0,
  color_nm_map = keystone::color_nm_map
)
{...}
```

This is the function signature of `plot_peak2_data()` with annotations:

```r
plot_peak2_data <- function(
  x,
  report_dir = NULL, # Optional path to report's output directory
  image_dir = report_dir,
  trans_fun = log10, # Also possibly 'identity'
  replace_with_na_condition = ~ is.infinite(.x) | is.nan(.x),
  max_cv_threshold = Inf, # Make anomalously high-valued CVs into missings
  ## If TRUE, keep only longest stretch of CVs without missing values:
  keep_longest_continuous = FALSE,
  remove_empty_cols = TRUE,
  save_png = FALSE, # If TRUE, save PNG plots to report directory
  png... = list(),
  default_color_fun = colorspace::rainbow_hcl,
  x_var_lab = c(PMT_voltage = "PMT Voltage"),
  y_var_lab = c(log10_CV = expression(paste(log[10], " CV"))),
  plot_series... = list(),
  create_smooth_variables... = list(),
  plot_derivative = FALSE,
  points... = list(),
  debug = FALSE,
  round_to_nearest_volts = 1, # Round to nearest multiple of this value
  round_any_fun = round,
  xlsx_expression = NULL,
  ## Noli me tangere; when TRUE, used for IRR checks:
  plot_individual_channels = FALSE
)
{...}
```
