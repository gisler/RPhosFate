# DEM related input

Clips, pre-processes and calculates or determines all input data related
to the digital elevation model (DEM) in the broader sense: *acc_inf,
cha, dem, dir_inf, rds, slp_inf,* and *wsh.*

Requires the
*[WhiteboxTools](https://www.whiteboxgeo.com/download-whiteboxtools/)*
binary
([`whitebox::install_whitebox`](https://whiteboxr.gishub.org/reference/install_whitebox.html))
to be installed on your computer.

## Usage

``` r
DEMrelatedInput(
  cv_dir,
  cs_dem,
  cs_cha,
  sp_msk,
  sp_olp,
  sp_sds,
  cs_rds = NULL,
  ns_cha = NULL,
  ns_brn = 50,
  is_adj = 1L,
  is_ths = 1L,
  ls_mD8 = FALSE,
  ls_tmp = FALSE
)
```

## Arguments

- cv_dir:

  A character vector specifying the desired project root directory
  (first position).

- cs_dem:

  A character string specifying a path to a potentially large raster
  digital elevation model.

- cs_cha:

  A character string specifying a path to a potentially large raster
  providing channels.

- sp_msk:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  providing a somewhat oversized catchment polygon mask used to clip the
  potentially large input rasters for further processing.

- sp_olp:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  providing the desired catchment outlet point(s).

- sp_sds:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  providing channel source points.

- cs_rds:

  An optional character string specifying a path to a potentially large
  raster providing roads.

- ns_cha:

  An optional numeric scalar specifying the minimum D8 flow accumulation
  in number of upslope grid cells determining a channel.

- ns_brn:

  A numeric scalar specifying the stream burning step size in m.

- is_adj:

  A numeric scalar specifying how many cells adjacent to channels shall
  be burnt.

- is_ths:

  An integer scalar specifying the number of threads to use for
  processing, where applicable.

- ls_mD8:

  A logical scalar specifying if D8 flow directions shall be mimicked,
  i.e. the D-infinity flow directions are rounded to the nearest
  multiple of 45 degrees. Please note that this treatment is always
  applied to channel cells independently of this argument.

- ls_tmp:

  A logical scalar specifying if the temporary files created during
  computation shall be kept.

## Value

A two column numeric [`matrix`](https://rdrr.io/r/base/matrix.html)
specifying one or more catchment outlet coordinates and side effects in
the form of raster files.

## Details

This function applies the following (pre-processing) steps to ensure
hydrologic consistency of the generated input data:

- Stream burning and orientation of cells adjacent to channel cells
  approximately into the direction of channel cells (no effect with
  `ns_brn = 0`).

- Depression breaching.

- Tracing of downslope flowpaths from the provided channel sources.

When roads are provided, they are considered as flow obstacles breaking
the continuity of the calculated flow accumulations.

`ns_cha` can be used to enhance the channel network obtained by the
tracing of downslope flowpaths from the provided channel sources.

*dem* represents the breached DEM with reversed stream burning if
applicable. The basis for the calculation of the D-infinity slopes
provided by *slp_inf,* however, is the original DEM.

## References

Lindsay, J.B., 2016. Efficient hybrid breaching-filling sink removal
methods for flow path enforcement in digital elevation models.
Hydrological Processes 30, 846–857. https://doi.org/10.1002/hyp.10648

Tarboton, D.G., 1997. A new method for the determination of flow
directions and upslope areas in grid digital elevation models. Water
Resour. Res. 33, 309–319. https://doi.org/10.1029/96WR03137

## See also

[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md),
[`catchment`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# obtain temporary project root directory
cv_dir <- normalizePath(
  tempfile("cmt"),
  winslash = .Platform$file.sep,
  mustWork = FALSE
)
# obtain directory holding "large" rasters and other required data sets
cs_dir_lrg <- system.file("tinytest", "largeData", package = "RPhosFate")

nm_olc <- DEMrelatedInput(
  cv_dir = cv_dir,
  cs_dem = file.path(cs_dir_lrg, "dem_lrg.tif"),
  cs_cha = file.path(cs_dir_lrg, "cha_lrg.tif"),
  sp_msk = terra::vect(file.path(cs_dir_lrg, "msk.shp")),
  sp_olp = terra::vect(file.path(cs_dir_lrg, "olp.shp")),
  sp_sds = terra::vect(file.path(cs_dir_lrg, "sds.shp")),
  cs_rds = file.path(cs_dir_lrg, "rds_lrg.tif"),
  ls_tmp = TRUE
)} # }
```
