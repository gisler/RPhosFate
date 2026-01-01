# Convert *ERDAS IMAGINE* to *GeoTIFF* raster files

Converts all *ERDAS IMAGINE* raster files in a directory and its
subdirectories into *GeoTIFF* raster files.

## Usage

``` r
img2tif(cs_dir, cs_crs = NULL)
```

## Arguments

- cs_dir:

  A character string specifying an existing directory.

- cs_crs:

  An optional character string used to set the coordinate reference
  system of all output raster files. See
  [`terra::crs`](https://rspatial.github.io/terra/reference/crs.html)
  for further information.

## Value

A character vector containing the paths to the processed *ERDAS IMAGINE*
raster files.
