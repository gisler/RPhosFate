# Patch release

* The current `raster` package does not compress _GeoTIFF_ raster files any longer by default. This is probably due to the switch from `rgdal` to `terra` and has been fixed by generally using the _LZW_ algorithm.
* Fixed warning "GDAL Message 6: driver GTiff does not support creation option COMPRESSED" curiously only occurring in the reference on GitHub Pages.
* Removed `hydroGOF` from imported packages list (`maptools` is retiring and `hydroGOF` depends on it via `hydroTSM`). Thanks to Roger Bivand for raising this issue (#17).
* Bumped minimum tested R version from 4.2.2 to 4.2.3 using the corresponding Posit Public Package Manager snapshot.
* Slightly improved documentation.

# Test environments

* Windows on GitHub Actions (4.2.3)
* Local Windows (release)
* Windows on GitHub Actions (devel)
* Linux on GitHub Actions (oldrel)
* Linux on GitHub Actions (release)

# R CMD check results

There were no ERRORs, WARNINGs or NOTEs.

# Downstream dependencies

There are currently no downstream dependencies for this package.
