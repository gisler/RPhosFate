# Subsequent run

Calls
[`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md)
for the specified substance and optionally
[`erosionPrerequisites`](https://gisler.github.io/RPhosFate/dev/reference/erosionPrerequisites-RPhosFate-method.md),
[`erosion`](https://gisler.github.io/RPhosFate/dev/reference/erosion-RPhosFate-method.md),
[`emission`](https://gisler.github.io/RPhosFate/dev/reference/emission-RPhosFate-method.md)
and
[`transportPrerequisites`](https://gisler.github.io/RPhosFate/dev/reference/transportPrerequisites-RPhosFate-method.md)
beforehand.

## Usage

``` r
# S4 method for class 'RPhosFate'
subsequentRun(
  x,
  substance = "PP",
  erosionPrerequisites = FALSE,
  erosion = FALSE,
  emission = FALSE,
  transportPrerequisites = FALSE
)
```

## Arguments

- x:

  An S4
  [`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
  river catchment object.

- substance:

  A character string specifying the substance to calculate.

- erosionPrerequisites:

  A logical scalar specifying if
  [`erosionPrerequisites`](https://gisler.github.io/RPhosFate/dev/reference/erosionPrerequisites-RPhosFate-method.md)
  is called.

- erosion:

  A logical scalar specifying if
  [`erosion`](https://gisler.github.io/RPhosFate/dev/reference/erosion-RPhosFate-method.md)
  is called.

- emission:

  A logical scalar specifying if
  [`emission`](https://gisler.github.io/RPhosFate/dev/reference/emission-RPhosFate-method.md)
  is called. It is never called with `substance = "SS"` though.

- transportPrerequisites:

  A logical scalar specifying if
  [`transportPrerequisites`](https://gisler.github.io/RPhosFate/dev/reference/transportPrerequisites-RPhosFate-method.md)
  is called.

## Value

An S4
[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
river catchment object and side effects in the form of raster files.

## See also

[`firstRun`](https://gisler.github.io/RPhosFate/dev/reference/firstRun-RPhosFate-method.md)

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
# presupposed method call
x <- firstRun(x, "SS")

x <- subsequentRun(x, "PP")# }
```
