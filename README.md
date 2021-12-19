An enhanced version of the PhosFate model implemented in R currently supporting suspended solids (SS) and particulate phosphorus (PP).

Copyright (C) 2021 RPhosFate authors

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

[![R build status](https://github.com/gisler/RPhosFate/workflows/R-CMD-check/badge.svg)](https://github.com/gisler/RPhosFate/actions?query=workflow%3AR-CMD-check) [![codecov](https://codecov.io/gh/gisler/RPhosFate/branch/main/graph/badge.svg?token=SrmeMV0ohC)](https://app.codecov.io/gh/gisler/RPhosFate) [![GitHub Super-Linter](https://github.com/gisler/RPhosFate/workflows/Lint%20Code%20Base/badge.svg)](https://github.com/gisler/RPhosFate/actions?query=workflow%3A%22Lint+Code+Base%22)

## Installation

Currently, no [CRAN](https://cran.r-project.org/) release is planned. A Windows binary package built for [Microsoft R Open](https://mran.microsoft.com/release-history) 3.5.3 is nonetheless available from a [`drat`](https://github.com/eddelbuettel/drat) repository. It can be installed along with its most important dependencies via the following command:

``` r
install.packages(
  c("RPhosFate", "raster", "whitebox"),
  repos = c("https://gisler.github.io/drat", options("repos")),
  type = "win.binary",
  dependencies = TRUE
)
```

## Semantic versioning

Releases of this project are versioned following the rules of [SemVer](https://semver.org).
