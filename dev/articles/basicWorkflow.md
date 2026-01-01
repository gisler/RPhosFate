# Basic modelling workflow

This vignette provides an overview of the basic modelling workflow with
the `RPhosFate` package.

------------------------------------------------------------------------

## Preparations

Load the package and obtain a copy of the demonstration project:

``` r
library(RPhosFate)

cv_dir <- demoProject()
```

## Project initialisation

Use
[`RPhosFate()`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md)
or
[`catchment()`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md)
to initialise the project:

``` r
x <- RPhosFate(
  cv_dir = cv_dir,
  ns_dep_ovl = 25e-4,
  ns_dep_cha = 0.0,
  nv_tfc_inl = c(SS = 0.6, PP = 0.6),
  nv_enr_rto = c(PP = 2.0),
  nm_olc = matrix(c(4704255, 2795195), ncol = 2L),
  df_cdt = read.table(
    file.path(cv_dir, "cdt.txt"),
    header = TRUE,
    stringsAsFactors = FALSE
  )
)
```

## First model run for suspended solids (SS)

[`firstRun()`](https://gisler.github.io/RPhosFate/dev/reference/firstRun-RPhosFate-method.md)
calls all low-level interface methods for the specified substance in the
required order:

``` r
x <- firstRun(x, substance = "SS")
```

## Calibration quality of SS

Snap coordinates of provided calibration gauges to the respective
midpoint of the nearest channel cell if necessary and check calibration
quality:

``` r
x <- snapGauges(x)

metrics <- calibrationQuality(x, substance = "SS", col = "SS_load")
```

## Calibrate SS

SS is calibrated by iteratively specifying better parameter values for
`ns_dep_ovl` (overland deposition rate) and/or `ns_dep_cha` (channel
deposition rate) as well as calling
[`subsequentRun()`](https://gisler.github.io/RPhosFate/dev/reference/subsequentRun-RPhosFate-method.md)
for SS afterwards until pleased with the metrics. By default,
[`subsequentRun()`](https://gisler.github.io/RPhosFate/dev/reference/subsequentRun-RPhosFate-method.md)
only calls the
[`transport()`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md)
low-level interface method for the specified substance:

``` r
x <- setParameter(x, ns_dep_ovl = 15e-4)

x <- subsequentRun(x, substance = "SS")

metrics <- calibrationQuality(x, substance = "SS", col = "SS_load")
```

The
[`autoCalibrate()`](https://gisler.github.io/RPhosFate/dev/reference/autoCalibrate-RPhosFate-method.md)
and
[`autoCalibrate2()`](https://gisler.github.io/RPhosFate/dev/reference/autoCalibrate2-RPhosFate-method.md)
methods may provide more comfortable alternatives to this process.

## Calibration quality of particulate phosphorus (PP)

First, a further call to
[`subsequentRun()`](https://gisler.github.io/RPhosFate/dev/reference/subsequentRun-RPhosFate-method.md)
for PP is necessary:

``` r
x <- subsequentRun(x, substance = "PP")

metrics <- calibrationQuality(x, substance = "PP", col = "PP_load")
```

## Calibrate PP

Same procedure as with SS apart from iteratively specifying better
parameter values for the enrichment ratio:

``` r
x <- setParameter(x, nv_enr_rto = c(PP = 1.4))

x <- subsequentRun(x, substance = "PP")

metrics <- calibrationQuality(x, substance = "PP", col = "PP_load")
```

In case the only substance of interest is PP, it is possible to set its
enrichment ratio to one and directly calibrate it via `ns_dep_ovl`
and/or `ns_dep_cha`.

## Save state

Write parameters to disk:

``` r
saveState(x)
```

------------------------------------------------------------------------

## Return to the calibrated project at a later time

Simply load the previously saved state of the project via the `ls_ini`
argument of
[`RPhosFate()`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md)
or
[`catchment()`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md):

``` r
x <- RPhosFate(
  cv_dir = cv_dir,
  ls_ini = TRUE
)
```
