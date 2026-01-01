# Save state

Saves parameters *(parameters.yaml)* to disk.

## Usage

``` r
# S4 method for class 'RPhosFate'
saveState(x)
```

## Arguments

- x:

  An S4
  [`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
  river catchment object.

## Value

`NULL` invisibly and side effects in the form of files.

## See also

[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md),
[`catchment`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md)

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

saveState(x)
```
