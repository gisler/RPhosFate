An enhanced version of the semi-empirical, spatially distributed emission and transport model PhosFate implemented in R and C++. It is based on the D-infinity, but also supports the D8 flow method. The currently available substances are suspended solids (SS) and particulate phosphorus (PP). A major feature is the allocation of substance loads entering surface waters to their sources of origin, which is a basic requirement for the identification of critical source areas and in consequence a cost-effective implementation of mitigation measures.

Copyright (C) 2021 RPhosFate authors

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

## Installation

Install the latest release from CRAN:

`install.packages("RPhosFate")`

[![CRAN Version](https://www.r-pkg.org/badges/version/RPhosFate)](https://cran.r-project.org/package=RPhosFate) [![CRAN Checks](https://badges.cranchecks.info/worst/RPhosFate.svg)](https://cran.r-project.org/web/checks/check_results_RPhosFate.html)

Install the [development version](https://gisler.github.io/RPhosFate/dev/) from GitHub (requires [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows and the `remotes` package):

`remotes::install_github("gisler/RPhosFate")`

[![R build status](https://github.com/gisler/RPhosFate/workflows/R-CMD-check/badge.svg)](https://github.com/gisler/RPhosFate/actions?query=workflow%3AR-CMD-check)

## Semantic versioning

Releases of this project are versioned following the rules of [SemVer](https://semver.org).
