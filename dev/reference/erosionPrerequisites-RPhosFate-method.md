# Erosion prerequisites

Calculates and writes capped slopes, L- and RUSLE S-factors (equations
for summer conditions and slopes \\\geq\\ 15 ft) to disk.

## Usage

``` r
# S4 method for class 'RPhosFate'
erosionPrerequisites(x)
```

## Arguments

- x:

  An S4
  [`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
  river catchment object.

## Value

An S4
[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
river catchment object and side effects in the form of raster files.

## References

Desmet, P.J.J., Govers, G., 1996. A GIS procedure for automatically
calculating the USLE LS factor on topographically complex landscape
units. Journal of Soil and Water Conservation 51, 427–433.

Renard, K.G., Foster, G.R., Weesies, G.A., McCool, D.K., Yoder, D.C.,
1997. Predicting soil erosion by water: a guide to conservation planning
with the Revised Universal Soil Loss Equation (RUSLE), Agriculture
Handbook. U.S. Government Printing Office, Washington, DC.

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

x <- erosionPrerequisites(x)# }
```
