# RPhosFate v0.9.9000

* Switched to utilising _GeoTIFF_ (\*.tif) instead of _Erdas Imagine_ (\*.img) raster files
* Added `cs_fex` argument to `DEMrelatedInput`: allows for creating _Erdas Imagine_ raster files for backward compatibility
* Added `cs_dir` argument to `DEMrelatedInput`: allows for utilising an existing D8 flow directions raster using _ArcGIS_ codes
* Added function `demoProject()` providing training data
* Added examples to documentation
* Added tests utilising the unit testing framework of the `tinytest` package
* Added means to measure code coverage with the help of `covr`
