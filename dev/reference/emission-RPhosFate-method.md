# Emission

Calculates and writes substance emissions to disk.

## Usage

``` r
# S4 method for class 'RPhosFate'
emission(x, substance = "PP")
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

x <- emission(x, "PP")# }
```
