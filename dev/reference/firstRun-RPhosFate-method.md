# First run

Calls
[`erosionPrerequisites`](https://gisler.github.io/RPhosFate/dev/reference/erosionPrerequisites-RPhosFate-method.md),
[`erosion`](https://gisler.github.io/RPhosFate/dev/reference/erosion-RPhosFate-method.md),
[`emission`](https://gisler.github.io/RPhosFate/dev/reference/emission-RPhosFate-method.md),
[`transportPrerequisites`](https://gisler.github.io/RPhosFate/dev/reference/transportPrerequisites-RPhosFate-method.md)
and
[`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md)
in the mentioned order. While
[`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md)
is called for the specified substance only,
[`emission`](https://gisler.github.io/RPhosFate/dev/reference/emission-RPhosFate-method.md)
is called for all substances whose top soil concentrations have been
provided.

## Usage

``` r
# S4 method for class 'RPhosFate'
firstRun(x, substance = "PP")
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

## See also

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

x <- firstRun(x, "SS")# }
```
