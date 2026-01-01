# Transport

Calculates and writes substance retentions, transports and cell loads as
well as transfers to disk.

## Usage

``` r
# S4 method for class 'RPhosFate'
transport(x, substance = "PP")
```

## Arguments

- x:

  An S4
  [`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
  river catchment object.

- substance:

  A character string specifying the substance to calculate.

## Value

An S4
[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
river catchment object and side effects in the form of raster files.

## References

Engman, E.T., 1986. Roughness coefficients for routing surface runoff.
Journal of Irrigation and Drainage Engineering 112, 39–53.

Molnár, P., Ramírez, J.A., 1998. Energy dissipation theories and optimal
channel characteristics of river networks. Water Resources Research 34,
1809–1818. https://doi.org/10.1029/98WR00983

## See also

[`firstRun`](https://gisler.github.io/RPhosFate/dev/reference/firstRun-RPhosFate-method.md),
[`subsequentRun`](https://gisler.github.io/RPhosFate/dev/reference/subsequentRun-RPhosFate-method.md)

## Examples

``` r
# \donttest{
# temporary demonstration project copy
cv_dir <- demoProject()
#> Warning: A folder called "demoProject" already exists and is left as is.
# load temporary demonstration project
x <- RPhosFate(
  cv_dir = cv_dir,
  ls_ini = TRUE
)
# presupposed method calls
x <- erosionPrerequisites(x)
x <- erosion(x)
x <- emission(x, "PP")
x <- transportPrerequisites(x)

x <- transport(x, "PP")# }
```
