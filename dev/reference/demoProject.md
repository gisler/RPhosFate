# Demonstration project

Copies a demonstration project to an existing or a temporary directory.

The demonstration project data are a derivative of the

- *[Geoland.at](https://www.data.gv.at/katalog/dataset/d88a1246-9684-480b-a480-ff63286b35b7)*
  (digital elevation model),

- *AMA* (field data; utilised dataset with the ID
  35e36014-ec69-439b-8629-389f52ffaa92 was removed from [Offene Daten
  Österreich](https://www.data.gv.at)),

- *[BMLRT](https://www.data.gv.at/katalog/dataset/c2287ccb-f44c-48cd-bf7c-ac107b771246)*
  (channel data) and

- *[GIP.at](https://www.data.gv.at/katalog/dataset/3fefc838-791d-4dde-975b-a4131a54e7c5)*
  (road data)

data sets, used and licensed under *[(CC BY
4.0)](https://creativecommons.org/licenses/by/4.0/)* by Gerold Hepp.

While the data represent a real catchment
*[(HOAL)](https://hoal.hydrology.at/),* some of them are fictitious, but
plausible. These are, among others, R- and C-factors, soil and related
data, existence of subsurface drainage at road embankments as well as
substance river loads.

## Usage

``` r
demoProject(cs_dir = tempdir(TRUE))
```

## Arguments

- cs_dir:

  An optional character string specifying an existing directory.

## Value

A character string containing the demonstration project root directory.

## See also

[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md),
[`catchment`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md)

## Examples

``` r
demoProject()
#> Warning: A folder called "demoProject" already exists and is left as is.
#> [1] "/tmp/Rtmphy0jp5/demoProject"
```
