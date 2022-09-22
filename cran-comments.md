# Major release

* First CRAN release

# Test environments

* Windows on GitHub Actions (4.1.2)
* Local Windows (4.2.0)
* Windows on GitHub Actions (devel)
* macOS on GitHub Actions (oldrel)
* macOS on GitHub Actions (release)
* Linux on GitHub Actions (oldrel)
* Linux on GitHub Actions (release)

# R CMD check results

Thank you for the constructive feedback. I fixed all the mentioned issues.

There were no ERRORs and WARNINGs.

## NOTES

  Found the following (possibly) invalid URLs:
    URL: https://cran.r-project.org/web/checks/check_results_RPhosFate.html
      From: README.md
      Status: 404
      Message: Not Found

> The URL will be available once the package is published.

NOTE only on Linux:

  checking installed package size ... NOTE
    installed size is  5.9Mb
    sub-directories of 1Mb or more:
      libs       3.6Mb
      tinytest   1.7Mb

> I hope the size of libs does not matter, as it merely contains compiled binaries.
>
> The tinytest sub-directory mainly contains many small GeoTIFFs, serving two purposes:
>
> * Demonstration data for the user
> * Test data
>
> Since the purpose of the package is to process many different sources of data (different aspects of elevation, soil, land use etc.) with a geographic reference, all the GeoTIFFs are required for the given purposes. I tried to compress them as much as possible though.

# Downstream dependencies

There are currently no downstream dependencies for this package.
