# Patch release

* Removed `rgdal` from suggested packages list and set minimum required version of the `raster` package to ≥ 3.6.3 (`rgdal` is retiring and `raster` ≥ 3.6.3 does not depend on it any longer). Thanks to Roger Bivand for raising this issue (#17).
* Bumped minimum tested R version from 4.1.2 to 4.2.2 using the corresponding MRAN repository snapshot.
* Slightly improved documentation.
* Minor internal code improvements.

# Test environments

* Windows on GitHub Actions (4.2.2)
* Local Windows (release)
* Windows on GitHub Actions (devel)
* Linux on GitHub Actions (oldrel)
* Linux on GitHub Actions (release)

# R CMD check results

There were no ERRORs, WARNINGs or NOTEs.

# Downstream dependencies

There are currently no downstream dependencies for this package.
