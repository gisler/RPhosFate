# RPhosFate v1.0.4.9006

This version introduces a new flow method, but also comes with several breaking changes.

## Overview

* Use of the D-infinity (Tarboton, [1997](https://doi.org/10.1029/96WR03137)) instead of the D8 flow method, however, the D8 flow method can still be mimicked by rounding D-infinity flow directions to the nearest multiple of 45 degrees.
* Weighted flow accumulations are no longer supported for the time being.
* The basis for the calculation of the D-infinity slopes is the original digital elevation model and not the breached one with reversed stream burning as it used to be for the calculation of the D8 slopes.
* The default value for the parameter holding the minimum bounding slope _(ns\_slp\_min)_ is now 1 instead of 0.001%.
* The L factor is now calculated using equation\ (9) as published by Desmet and Govers ([1996](https://www.jswconline.org/content/51/5/427)) and not as modified by Kovacs ([2013](http://hdl.handle.net/20.500.12708/9468)).
* The channel retentions of the substance outlet loads of subsurface drainages are now calculated analogous to the overland retentions of the local emissions, i.e. by using half of the respective flow path length.
* `RPhosFate` now utilises the `terra` instead of the `raster` package.
* Ceased the support for _ERDAS IMAGINE_ (\*.img) raster files and added `img2tif()` function can be used to convert all _ERDAS IMAGINE_ raster files in a directory and its subdirectories into _GeoTIFF_ raster files.
* Dropped backward compatibility to major version zero.

## Technical details

* Input data changes:
  * Removed the layers holding the (weighted) D8 flow accumulations (_acc_ and _acc\_wtd_) and added _acc\_inf_ holding the D-infinity flow accumulations.
  * Removed the layer holding the D8 flow directions _(dir)_ and added _dir\_inf_ holding the D-infinity flow directions.
  * Removed the layer holding the D8 slopes _(slp)_ and added _slp\_inf_ holding the D-infinity slopes.
* `DEMrelatedInput()` function changes:
  * Considerable adjustments to reflect the aforementioned input data changes.
  * Addition of `ns_cha` argument: allows for specifying the minimum D8 flow accumulation determining a channel.
  * Addition of `ls_mD8` argument: allows for mimicking D8 flow directions by rounding the D-infinity flow directions to the nearest multiple of 45 degrees. Please note that this treatment is always applied to channel cells independently of this argument.
  * Removal of the `cs_wgs` and `cs_dir` arguments: These input data are no longer supported for the time being.
* Switched to utilising the `SpatRaster` and `SpatVector` classes from the `terra` package instead of the `RasterLayer` class from the `raster` and the `Spatial*DataFrame` classes from the `sp` packages.
* Added `is_ths` argument to the `RPhosFate()` and `catchment()` constructors: allows for specifying the number of threads to use for processing, where applicable.
* Removed the layer holding the hydraulic radii _(rhy),_ as the calculation of the hydraulic radii is now integrated into the `transport()` method. This implies that the `transportPrerequisites()` method does not save it to disk any longer.
* Removed the parameter holding the D8 outflow direction vector _(iv\_fDo)._ Existing parameter files containing it can still be used, but the parameter will definitely be lost upon saving the project's state.
* Removed the `transportCalcOrder()` method, as the determination of the cell transport calculation order is now integrated into the `transport()` method. This implies that the `firstRun()` as well as `subsequentRun()` methods do not call it and the `saveState()` method does not save it to disk any longer.
* The `calibrationQuality()` method now returns its return value invisibly.
* Considerably revised the internal `RPhosFateHelpers` class.
* Removed `spatstat.geom` from the imported packages list (utilised functionality is now also provided by `terra`).
* Bumped the minimum tested R version from 4.2.3 to 4.3.2 using the corresponding _Posit_ public package manager snapshot.
* Slightly improved the documentation.
* Major internal code improvements.

# RPhosFate v1.0.4

* Removed `hydroGOF` from the imported packages list (`maptools` is retiring and `hydroGOF` depends on it via `hydroTSM`). Thanks to Roger Bivand for raising this issue (#17).
* Bumped the minimum tested R version from 4.2.2 to 4.2.3 using the corresponding _Posit_ public package manager snapshot.
* Slightly improved the documentation.

# RPhosFate v1.0.3

* Removed `rgdal` from the suggested packages list and set the minimum required version of the `raster` package to ≥ 3.6.3 (`rgdal` is retiring and `raster` ≥ 3.6.3 does not depend on it any longer). Thanks to Roger Bivand for raising this issue (#17).
* Bumped the minimum tested R version from 4.1.2 to 4.2.2 using the corresponding _MRAN_ repository snapshot.
* Slightly improved the documentation.
* Minor internal code improvements.

# RPhosFate v1.0.2

* Removed _NRMSE_ from calibration quality metrics and added _KGE_ as well as _RCV._

# RPhosFate v0.12.0

* Added `autoCalibrate2()` method: allows for calibrating the overland and channel deposition rate in one go.
* The `DEMrelatedInput()` function can handle multiple catchment outlets now and so does the `calibrationQuality()` method.
* The `DEMrelatedInput()` function now returns the breached DEM with reversed stream burning if applicable instead of the original one.
* The `DEMrelatedInput()` function now calculates “correct” slopes even if the channels used for stream burning contain gaps.
* `RPhosFate` now makes sure that the x- and y-coordinates of gauges used for calibration lie within the extent of the river catchment object.
* Renamed the `"inChannelRetention"` output of the `calibrationQuality()` method to `"inChannelRetentionRatio"`.
* `RPhosFate` requires R ≥ 3.5.0 now.
* Improved the documentation.

# RPhosFate v0.11.0

* Monte Carlo simulation mode can now make use of repeated random samples, i.e. raster data, of distributions of about all input data.
* Added the following arguments to the `subsequentRun()` method, which allow for calling the respective methods:
  * `erosionPrerequisites`
  * `erosion`
  * `emission`
  * `transportPrerequisites`
  * `transportCalcOrder`
* Added `cv_MCl` argument to the `RPhosFate()` and `catchment()` constructors: allows for specifying the names of the layers, which shall be written to disk with the associated Monte Carlo iteration in their filenames upon calling the appropriate methods.
* Initialising a project in Monte Carlo simulation mode now also reads model results produced by a possible earlier run associated with the specified iteration. This implies that Monte Carlo input data additionally can reside in the project root subdirectories and not only in a separate directory.
* The plot produced by the `calibrationQuality()` method is now prettier.
* Added a vignette describing the basic modelling workflow.
* `RPhosFate` now depends on the `spatstat.geom` instead of the `spatstat` package.
* Fixed the minimum required version of the `whitebox` package (≥ 2.0.0).
* Added test for the standard use case of the `DEMrelatedInput()` function.
* Slightly improved the documentation.

# RPhosFate v0.10.0

* The `DEMrelatedInput()` function now calculates the slopes from the breached DEM (stream burning is undone beforehand).
* Switched to utilising _GeoTIFF_ (\*.tif) instead of _ERDAS IMAGINE_ (\*.img) raster files.
* Added `cs_fex` argument to the `DEMrelatedInput()` function: allows for using _ERDAS IMAGINE_ raster files for backward compatibility.
* Added `cs_dir` argument to the `DEMrelatedInput()` function: allows for utilising an existing D8 flow directions raster using _ArcGIS_ codes.
* Added `demoProject()` function providing training data.
* Added examples to the documentation.
* Fixed a backward incompatibility.
* Added tests utilising the unit testing framework of the `tinytest` package.
* Added means to measure code coverage with the help of `covr`.
