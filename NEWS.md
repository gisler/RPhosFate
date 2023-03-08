# RPhosFate v1.0.3.9000

* The current `raster` package does not compress _GeoTIFF_ raster files any longer by default. This is probably due to the switch from `rgdal` to `terra` and has been fixed by generally using the _LZW_ algorithm.
* Fixed warning "GDAL Message 6: driver GTiff does not support creation option COMPRESSED" curiously only occurring in the reference on GitHub Pages.

# RPhosFate v1.0.3

* Removed `rgdal` from suggested packages list and set minimum required version of the `raster` package to ≥ 3.6.3 (`rgdal` is retiring and `raster` ≥ 3.6.3 does not depend on it any longer). Thanks to Roger Bivand for raising this issue (#17).
* Bumped minimum tested R version from 4.1.2 to 4.2.2 using the corresponding MRAN repository snapshot.
* Slightly improved documentation.
* Minor internal code improvements.

# RPhosFate v1.0.2

* Removed _NRMSE_ from calibration quality metrics and added _KGE_ as well as _RCV._

# RPhosFate v0.12.0

* Added `autoCalibrate2()` method: allows for calibrating the overland and channel deposition rate in one go.
* `DEMrelatedInput()` can handle multiple catchment outlets now and so does `calibrationQuality()`.
* `DEMrelatedInput()` now returns the breached DEM with reversed stream burning if applicable instead of the original one.
* `DEMrelatedInput()` now calculates “correct” slopes even if the channels used for stream burning contain gaps.
* `RPhosFate` now makes sure that the x- and y-coordinates of gauges used for calibration lie within the extent of the river catchment object.
* Renamed `"inChannelRetention"` output of `calibrationQuality()` to `"inChannelRetentionRatio"`.
* `RPhosFate` requires R ≥ 3.5.0 now.
* Improved documentation.

# RPhosFate v0.11.0

* Monte Carlo simulation mode can now make use of repeated random samples, i.e. raster data, of distributions of about all input data.
* Added the following arguments to the `subsequentRun()` method, which allow for calling the respective methods:
  * `erosionPrerequisites`
  * `erosion`
  * `emission`
  * `transportPrerequisites`
  * `transportCalcOrder`
* Added `cv_MCl` argument to `RPhosFate()` and `catchment()` methods: allows for specifying the names of the layers, which shall be written to disk with the associated Monte Carlo iteration in their filenames upon calling the appropriate methods.
* Initialising a project in Monte Carlo simulation mode now also reads model results produced by a possible earlier run associated with the specified iteration. This implies that Monte Carlo input data additionally can reside in the project root subdirectories and not only in a separate directory.
* Plot produced by `calibrationQuality()` is now prettier.
* Added a vignette describing the basic modelling workflow.
* `RPhosFate` now depends on the `spatstat.geom` instead of the `spatstat` package.
* Fixed minimum required version of the `whitebox` package (≥ 2.0.0).
* Added test for the standard use case of the `DEMrelatedInput()` function.
* Slightly improved documentation.

# RPhosFate v0.10.0

* `DEMrelatedInput()` function now calculates the slopes from the breached DEM (stream burning is undone beforehand).
* Switched to utilising _GeoTIFF_ (\*.tif) instead of _Erdas Imagine_ (\*.img) raster files.
* Added `cs_fex` argument to `DEMrelatedInput()` function: allows for creating _Erdas Imagine_ raster files for backward compatibility.
* Added `cs_dir` argument to `DEMrelatedInput()` function: allows for utilising an existing D8 flow directions raster using _ArcGIS_ codes.
* Added `demoProject()` function providing training data.
* Added examples to documentation.
* Fixed backward incompatibility.
* Added tests utilising the unit testing framework of the `tinytest` package.
* Added means to measure code coverage with the help of `covr`.
