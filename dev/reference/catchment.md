# Initialise project

Initialises a project from scratch or loads the state of an existing one
utilising *GeoTIFF* (\*.tif) raster files from, by convention, the
following three project root subdirectories:

- *Input*

- *Intermediate*

- *Result*

See subdirectory sections for further information.

`catchment` is an alias for `RPhosFate`.

## Usage

``` r
RPhosFate(...)

catchment(...)
```

## Arguments

- ...:

  Arguments used to initialise the project. See argument sections for
  further information.

## Value

An S4
[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
river catchment object.

## *Input* subdirectory

This directory holds all possible user input raster data (flow obstacles
like roads must be considered during the generation of the flow
accumulation layer and must also be cut out from it in order to be
properly respected):

- *acc_inf:* D-infinity flow accumulations in number of upslope grid
  cells required for everything.

- *CFa:* (R)USLE C-factors required for
  [`erosion`](https://gisler.github.io/RPhosFate/dev/reference/erosion-RPhosFate-method.md).

- *cha:* Channel cells required for everything (`1`: channel cell, `NA`:
  no channel cell).

- *clc:* Clay contents of top soils in % required for substance
  [`emission`](https://gisler.github.io/RPhosFate/dev/reference/emission-RPhosFate-method.md)s.

- *dem:* Digital elevation model in m a.s.l. (optional).

- *dir_inf:* D-infinity flow directions in azimuth degrees measured from
  north (0 to 360 clockwise) required for
  [`transportPrerequisites`](https://gisler.github.io/RPhosFate/dev/reference/transportPrerequisites-RPhosFate-method.md)
  and substance
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md).

- *fid:* Field IDs (optional).

- *KFa:* (R)USLE K-factors required for
  [`erosion`](https://gisler.github.io/RPhosFate/dev/reference/erosion-RPhosFate-method.md).

- *lue:* Land use classes (optional).

- *man:* Manning's roughness coefficients required for substance
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md).

- *xxc:* Substance contents of top soils in mg/kg required for substance
  [`emission`](https://gisler.github.io/RPhosFate/dev/reference/emission-RPhosFate-method.md)s,
  for example, *ppc* for PP top soil contents.

- *rds:* Road cells required for
  [`transportPrerequisites`](https://gisler.github.io/RPhosFate/dev/reference/transportPrerequisites-RPhosFate-method.md)
  (`0`: road cell without subsurface drainage, `1`: road cell with
  subsurface drainage, `NA`: no road cell).

- *RFa:* (R)USLE R-factors required for
  [`erosion`](https://gisler.github.io/RPhosFate/dev/reference/erosion-RPhosFate-method.md).

- *slp_inf:* D-infinity slopes in % required for everything.

- *wsh:* Watershed (optional).

## *Intermediate* subdirectory

This directory holds intermediate calculations:

- *inl:* Cells representing inlets at roads (storm drains).

- *LFa:* L-factors.

- *rip:* Cells representing the riparian zones within channel cells.

- *SFa:* RUSLE S-factors.

- *slp_cap:* Capped slopes in %.

## *Result* subdirectory

This directory holds the model results:

- *ero:* Erosion in t/cell/yr.

- *xxe:* Substance emissions in kg/cell/yr, for example, *ppe* for PP
  emissions.

- *xxr:* Substance retentions in t/cell/yr (SS) or kg/cell/yr, for
  example, *ppr* for PP retentions.

- *xxt:* Substance transports in t/cell/yr (SS) or kg/cell/yr, for
  example, *ppt* for PP transports.

- *xxt_cld:* Substance cell loads in t/cell/yr (SS) or kg/cell/yr, for
  example, *ppt_cld* for PP cell loads.

- *xxt_ctf:* Substance cell transfers in t/cell/yr (SS) or kg/cell/yr,
  for example, *ppt_ctf* for PP transfers.

- *xxt_inp:* Substance inputs into surface waters in t/cell/yr (SS) or
  kg/cell/yr, for example, *ppt_inp* for PP inputs into surface waters.

- *xxt_out:* Substance outlet loads of subsurface drainages in t/cell/yr
  (SS) or kg/cell/yr, for example, *ppt_out* for PP outlet loads.

## Data management and processing arguments

- `cv_dir`: A character vector specifying the project root (first
  position) and optionally the Monte Carlo input data directory (second
  position).

- `ls_ini`: A logical scalar specifying if the state of an existing
  project shall be loaded from disk (defaults to `FALSE`). Parameters or
  substance parameter values specified via the `...` argument take
  precedence over loaded ones.

- `is_ths`: An integer scalar holding the number of threads to use for
  processing (defaults to 1).

- `is_MCi`: An integer scalar specifying the current Monte Carlo
  iteration if applicable (defaults to
  [`integer()`](https://rdrr.io/r/base/integer.html), which means Monte
  Carlo simulation mode is disabled).

- `cv_MCl`: A character vector specifying the names of the layers, which
  shall be written to disk with the associated Monte Carlo iteration in
  their filenames upon calling the appropriate methods (defaults to
  `"xxt"`; no effect in case Monte Carlo simulation mode is disabled).

## Model parameter arguments

- `ns_slp_min`: A numeric scalar specifying the minimum bounding slope
  in % (defaults to `1.0`).

- `ns_slp_max`: A numeric scalar specifying the maximum bounding slope
  in % (defaults to `999.0`).

- `ns_rhy_a`: A numeric scalar specifying a network constant depending
  on the discharge frequency needed for the calculation of the hydraulic
  radius, which in turn is a prerequisite for substance transport
  (defaults to `0.09` representing a discharge frequency of
  approximately six years).

- `ns_rhy_b`: A numeric scalar specifying a geometry scaling exponent
  depending on the discharge frequency needed for the calculation of the
  hydraulic radius, which in turn is a prerequisite for substance
  transport (defaults to `0.50` representing a discharge frequency of
  approximately six years).

- `ns_cha_rto`: A numeric scalar specifying the ratio of the channel to
  the cell width determining the widths of the riparian zones required
  for substance
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md)
  (defaults to `0.5`).

- `ns_man_rip`: A numeric scalar specifying Manning's roughness
  coefficient of the riparian zones within channel cells required for
  substance
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md)
  (defaults to `0.32`).

- `ns_man_cha`: A numeric scalar specifying Manning's roughness
  coefficient of the channel within channel cells required for substance
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md)
  (defaults to `0.04`).

- `ns_dep_ovl`: A numeric scalar specifying the overland deposition rate
  per second required for substance
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md)
  (calibration parameter; no default).

- `ns_dep_cha`: A numeric scalar specifying the channel deposition rate
  per second required for substance
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md)
  (calibration parameter; no default).

- `nv_tfc_inl`: A named numeric vector specifying the inlet transfer
  coefficients required for substance
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md),
  for example, `c(SS = 0.6, PP = 0.6)` (no default).

- `nv_enr_rto` A named numeric vector specifying the substance
  enrichment ratios required for substance except SS
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md),
  for example, `c(PP = 2.0)` (calibration parameter; no default).

- `nm_olc`: A two column numeric
  [`matrix`](https://rdrr.io/r/base/matrix.html) specifying one or more
  catchment outlet coordinates required for the in-channel retention
  ratio of
  [`calibrationQuality`](https://gisler.github.io/RPhosFate/dev/reference/calibrationQuality-RPhosFate-method.md)
  (no default).

- `df_cdt`: A [`data.frame`](https://rdrr.io/r/base/data.frame.html)
  with calibration data, which must have at least the following three
  columns and one or more columns with substance river loads in t/yr (no
  default):

  - *ID:* ID(s) of the gauge(s)

  - *x:* x-coordinate(s) of the gauge(s)

  - *y:* y-coordinate(s) of the gauge(s)

## Monte Carlo simulation mode

This mode can make use of repeated random samples, i.e. raster data, of
distributions of about all input data. The filenames of the Monte Carlo
input raster data must contain the specified iteration, for example,
*CFa12.tif* for the twelfth iteration of the C-factors input data, and
can reside in a separate directory. In case no Monte Carlo raster file
is found for a certain layer in the designated directory, the respective
project root subdirectory is searched for one and finally the “normal”
project input raster data is utilised.

## See also

[`saveState`](https://gisler.github.io/RPhosFate/dev/reference/saveState-RPhosFate-method.md),
[`demoProject`](https://gisler.github.io/RPhosFate/dev/reference/demoProject.md)

## Examples

``` r
# \donttest{
# temporary demonstration project copy
cv_dir <- demoProject()
#> Warning: A folder called "demoProject" already exists and is left as is.

# initialise project from scratch
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

# load state of existing project in Monte Carlo simulation mode
x <- RPhosFate(
  cv_dir = c(
    cv_dir,
    system.file("tinytest", "testProject", package = "RPhosFate")
  ),
  ls_ini = TRUE,
  is_MCi = 1L,
  cv_MCl = c("xxt", "xxt_cld")
)# }
```
