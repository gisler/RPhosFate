# Get parameter(s)

Obtains a single model parameter or all model parameters at once.

## Usage

``` r
# S4 method for class 'RPhosFate'
getParameter(x, parameter = NULL)
```

## Arguments

- x:

  An S4
  [`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
  river catchment object.

- parameter:

  A character string specifying a parameter name or `NULL` for a
  [`list`](https://rdrr.io/r/base/list.html) of all parameters. See
  model parameter arguments section for further information.

## Value

Depends on the queried parameter or a
[`list`](https://rdrr.io/r/base/list.html) in case of all parameters.
See model parameter arguments section for further information.

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

[`setParameter`](https://gisler.github.io/RPhosFate/dev/reference/setParameter-RPhosFate-method.md)

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

getParameter(x)
#> $ns_slp_min
#> [1] 1
#> 
#> $ns_slp_max
#> [1] 999
#> 
#> $ns_rhy_a
#> [1] 0.09
#> 
#> $ns_rhy_b
#> [1] 0.5
#> 
#> $ns_cha_rto
#> [1] 0.5
#> 
#> $ns_man_rip
#> [1] 0.32
#> 
#> $ns_man_cha
#> [1] 0.04
#> 
#> $ns_dep_ovl
#> [1] 0.002
#> 
#> $ns_dep_cha
#> [1] 0
#> 
#> $nv_tfc_inl
#>  SS  PP 
#> 0.6 0.6 
#> 
#> $nv_enr_rto
#> PP 
#>  2 
#> 
#> $nm_olc
#>         [,1]    [,2]
#> [1,] 4704255 2795195
#> 
#> $df_cdt
#>   ID       x       y SS_load  PP_load
#> 1 G3 4704246 2795192  28.082 0.027168
#> 2 G2 4704193 2795374  19.425 0.020549
#> 3 G1 4704054 2795604   2.387 0.002127
#> 
getParameter(x, "ns_dep_ovl")
#> [1] 0.002
```
