# Set parameter(s)

Sets one or more model parameters or substance parameter values.

## Usage

``` r
# S4 method for class 'RPhosFate'
setParameter(x, ...)
```

## Arguments

- x:

  An S4
  [`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
  river catchment object.

- ...:

  Names and values of the parameters to set. See model parameter
  arguments section for further information.

## Value

An S4
[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
river catchment object.

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

## See also

[`getParameter`](https://gisler.github.io/RPhosFate/dev/reference/getParameter-RPhosFate-method.md)

## Examples

``` r
# temporary demonstration project copy
cv_dir <- demoProject()
#> Warning: A folder called "demoProject" already exists and is left as is.
# load temporary demonstration project
x <- RPhosFate(
  cv_dir = cv_dir,
  ls_ini = TRUE
)

x <- setParameter(x, ns_dep_ovl = 15e-4)
x <- setParameter(
  x,
  nv_tfc_inl = c(SS = 0.6, PP = 0.6),
  nv_enr_rto = c(PP = 1.4)
)
```
