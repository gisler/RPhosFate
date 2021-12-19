# RPhosFate v0.10.9000

* CI: check main sub-architecture only

# RPhosFate v0.10.0

* `DEMrelatedInput()` function now calculates the slopes from the breached DEM (stream burning is undone beforehand)
* Switched to utilising _GeoTIFF_ (\*.tif) instead of _Erdas Imagine_ (\*.img) raster files
* Added `cs_fex` argument to `DEMrelatedInput()` function: allows for creating _Erdas Imagine_ raster files for backward compatibility
* Added `cs_dir` argument to `DEMrelatedInput()` function: allows for utilising an existing D8 flow directions raster using _ArcGIS_ codes
* Added `demoProject()` function providing training data
* Added examples to documentation
* Added tests utilising the unit testing framework of the `tinytest` package
* Added means to measure code coverage with the help of `covr`
* Fixed backward incompatibility
