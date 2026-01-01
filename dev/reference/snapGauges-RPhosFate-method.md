# Snap gauge(s)

Snaps the coordinates of the provided calibration gauges to the
respective midpoint of the nearest channel cell.

## Usage

``` r
# S4 method for class 'RPhosFate'
snapGauges(x)
```

## Arguments

- x:

  An S4
  [`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
  river catchment object.

## Value

An S4
[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
river catchment object.

## See also

[`calibrationQuality`](https://gisler.github.io/RPhosFate/dev/reference/calibrationQuality-RPhosFate-method.md),
[`autoCalibrate`](https://gisler.github.io/RPhosFate/dev/reference/autoCalibrate-RPhosFate-method.md),
[`autoCalibrate2`](https://gisler.github.io/RPhosFate/dev/reference/autoCalibrate2-RPhosFate-method.md)

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

x <- snapGauges(x)
```
