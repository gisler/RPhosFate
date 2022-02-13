# RPhosFate v0.11.0.9000

* Added `autoCalibrate2()` method: allows for calibrating the overland and channel deposition rate in one go
* `DEMrelatedInput()` can handle multiple catchment outlets now and so does `calibrationQuality()`
* Renamed `"inChannelRetention"` output of `calibrationQuality()` to `"inChannelRetentionRatio"`
* `RPhosFate` requires R ≥ 3.5.0 now
* Improved documentation

# RPhosFate v0.11.0

* Monte Carlo simulation mode can now make use of repeated random samples, i.e. raster data, of distributions of about all input data
* Added the following arguments to the `subsequentRun()` method, which allow for calling the respective methods:
  * `erosionPrerequisites`
  * `erosion`
  * `emission`
  * `transportPrerequisites`
  * `transportCalcOrder`
* Added `cv_MCl` argument to `RPhosFate()` and `catchment()` methods: allows for specifying the names of the layers, which shall be written to disk with the associated Monte Carlo iteration in their filenames upon calling the appropriate methods
* Initialising a project in Monte Carlo simulation mode now also reads model results produced by a possible earlier run associated with the specified iteration. This implies that Monte Carlo input data additionally can reside in the project root subdirectories and not only in a separate directory.
* Plot produced by `calibrationQuality()` is now prettier
* Added a vignette describing the basic modelling workflow
* `RPhosFate` now depends on the `spatstat.geom` instead of the `spatstat` package
* Fixed minimum required version of the `whitebox` package (≥ 2.0.0)
* Added test for the standard use case of the `DEMrelatedInput()` function
* Slightly improved documentation

# RPhosFate v0.10.0

* `DEMrelatedInput()` function now calculates the slopes from the breached DEM (stream burning is undone beforehand)
* Switched to utilising _GeoTIFF_ (\*.tif) instead of _Erdas Imagine_ (\*.img) raster files
* Added `cs_fex` argument to `DEMrelatedInput()` function: allows for creating _Erdas Imagine_ raster files for backward compatibility
* Added `cs_dir` argument to `DEMrelatedInput()` function: allows for utilising an existing D8 flow directions raster using _ArcGIS_ codes
* Added `demoProject()` function providing training data
* Added examples to documentation
* Fixed backward incompatibility
* Added tests utilising the unit testing framework of the `tinytest` package
* Added means to measure code coverage with the help of `covr`
