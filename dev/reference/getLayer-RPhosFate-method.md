# Get layer

Obtains a project raster layer for further analysis.

## Usage

``` r
# S4 method for class 'RPhosFate'
getLayer(x, i, j = NULL)

# S4 method for class 'RPhosFate,ANY,ANY'
x[i, j]
```

## Arguments

- x:

  An S4
  [`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
  river catchment object.

- i:

  A character string specifying a layer name. Substance related layers
  whose names start with *xx* are treated differently. They have to be
  queried by their name (not filename), for example, `"xxc"` in
  combination with `"PP"` in argument `j` queries the particulate
  phosphorus concentrations in top soils. See subdirectory sections for
  further information.

- j:

  A character string specifying a substance if applicable.

## Value

A
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object.

## *Input* subdirectory

This directory holds all possible user input raster data (flow obstacles
like roads must be considered during the generation of the flow
accumulation layer and must also be cut out from it in order to be
properly respected):

- *acc_inf:* D-infinity flow accumulations in number of upslope grid
  cells required for everything.

- *CFa:* (R)USLE C-factors required for
  [`erosion`](https://gisler.github.io/RPhosFate/dev/reference/erosion-RPhosFate-method.md).

- *cha:* Channel cells required for everything (`1`: channel cell, `NA`:
  no channel cell).

- *clc:* Clay contents of top soils in % required for substance
  [`emission`](https://gisler.github.io/RPhosFate/dev/reference/emission-RPhosFate-method.md)s.

- *dem:* Digital elevation model in m a.s.l. (optional).

- *dir_inf:* D-infinity flow directions in azimuth degrees measured from
  north (0 to 360 clockwise) required for
  [`transportPrerequisites`](https://gisler.github.io/RPhosFate/dev/reference/transportPrerequisites-RPhosFate-method.md)
  and substance
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md).

- *fid:* Field IDs (optional).

- *KFa:* (R)USLE K-factors required for
  [`erosion`](https://gisler.github.io/RPhosFate/dev/reference/erosion-RPhosFate-method.md).

- *lue:* Land use classes (optional).

- *man:* Manning's roughness coefficients required for substance
  [`transport`](https://gisler.github.io/RPhosFate/dev/reference/transport-RPhosFate-method.md).

- *xxc:* Substance contents of top soils in mg/kg required for substance
  [`emission`](https://gisler.github.io/RPhosFate/dev/reference/emission-RPhosFate-method.md)s,
  for example, *ppc* for PP top soil contents.

- *rds:* Road cells required for
  [`transportPrerequisites`](https://gisler.github.io/RPhosFate/dev/reference/transportPrerequisites-RPhosFate-method.md)
  (`0`: road cell without subsurface drainage, `1`: road cell with
  subsurface drainage, `NA`: no road cell).

- *RFa:* (R)USLE R-factors required for
  [`erosion`](https://gisler.github.io/RPhosFate/dev/reference/erosion-RPhosFate-method.md).

- *slp_inf:* D-infinity slopes in % required for everything.

- *wsh:* Watershed (optional).

## *Intermediate* subdirectory

This directory holds intermediate calculations:

- *inl:* Cells representing inlets at roads (storm drains).

- *LFa:* L-factors.

- *rip:* Cells representing the riparian zones within channel cells.

- *SFa:* RUSLE S-factors.

- *slp_cap:* Capped slopes in %.

## *Result* subdirectory

This directory holds the model results:

- *ero:* Erosion in t/cell/yr.

- *xxe:* Substance emissions in kg/cell/yr, for example, *ppe* for PP
  emissions.

- *xxr:* Substance retentions in t/cell/yr (SS) or kg/cell/yr, for
  example, *ppr* for PP retentions.

- *xxt:* Substance transports in t/cell/yr (SS) or kg/cell/yr, for
  example, *ppt* for PP transports.

- *xxt_cld:* Substance cell loads in t/cell/yr (SS) or kg/cell/yr, for
  example, *ppt_cld* for PP cell loads.

- *xxt_ctf:* Substance cell transfers in t/cell/yr (SS) or kg/cell/yr,
  for example, *ppt_ctf* for PP transfers.

- *xxt_inp:* Substance inputs into surface waters in t/cell/yr (SS) or
  kg/cell/yr, for example, *ppt_inp* for PP inputs into surface waters.

- *xxt_out:* Substance outlet loads of subsurface drainages in t/cell/yr
  (SS) or kg/cell/yr, for example, *ppt_out* for PP outlet loads.

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

getLayer(x, "dir_inf")
#> class       : SpatRaster 
#> size        : 92, 155, 1  (nrow, ncol, nlyr)
#> resolution  : 10, 10  (x, y)
#> extent      : 4702990, 4704540, 2795190, 2796110  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs 
#> source      : dir_inf.tif 
#> name        :     dir_inf 
#> min value   :   0.2138069 
#> max value   : 360.0000000 
getLayer(x, "xxt", "SS")
#> class       : SpatRaster 
#> size        : 92, 155, 1  (nrow, ncol, nlyr)
#> resolution  : 10, 10  (x, y)
#> extent      : 4702990, 4704540, 2795190, 2796110  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs 
#> source      : sst.tif 
#> name        :          sst 
#> min value   : 3.027626e-06 
#> max value   : 1.754717e+01 
getLayer(x, "xxe", "PP")# }
#> class       : SpatRaster 
#> size        : 92, 155, 1  (nrow, ncol, nlyr)
#> resolution  : 10, 10  (x, y)
#> extent      : 4702990, 4704540, 2795190, 2796110  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs 
#> source      : ppe.tif 
#> name        :          ppe 
#> min value   : 5.648591e-05 
#> max value   : 9.866220e-01 
```
