# One dimensional automatic model calibration

Automatically calibrates the model with the help of a combination of
golden section search and successive parabolic interpolation.

## Usage

``` r
# S4 method for class 'RPhosFate'
autoCalibrate(
  x,
  substance,
  col,
  interval,
  metric,
  tol = min(interval) * 0.1,
  parameter = NULL
)
```

## Arguments

- x:

  An S4
  [`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
  river catchment object.

- substance:

  A character string specifying the substance to calculate.

- col:

  A character string specifying the calibration data column with the
  respective substance river loads.

- interval:

  A numeric vector specifying the end-points of the interval to be
  searched.

- metric:

  A character string specifying the metric to optimise. See
  [`calibrationQuality`](https://gisler.github.io/RPhosFate/dev/reference/calibrationQuality-RPhosFate-method.md)
  for available metrics.

- tol:

  A numeric scalar specifying the desired accuracy of the parameter used
  for optimisation (not the metric).

- parameter:

  By default, SS are calibrated utilising the overland deposition rate
  and all other substances are calibrated utilising their respective
  enrichment ratio. This argument can be used to specify a dedicated
  parameter utilised for calibration via a character string:
  `"ns_dep_ovl"` for overland or `"ns_dep_cha"` for channel deposition
  rate.

## Value

An S4
[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
river catchment object and side effects in the form of raster files.

## See also

[`snapGauges`](https://gisler.github.io/RPhosFate/dev/reference/snapGauges-RPhosFate-method.md),
[`optimize`](https://rdrr.io/r/stats/optimize.html)

## Examples

``` r
# \donttest{
# temporary demonstration project copy
cv_dir <- demoProject()
# load temporary demonstration project
x <- RPhosFate(
  cv_dir = cv_dir,
  ls_ini = TRUE
)
# presupposed method calls
x <- firstRun(x, "SS")
x <- snapGauges(x)

x <- autoCalibrate(
  x,
  "SS",
  col = "SS_load",
  interval = c(1e-3, 2e-3),
  metric = "KGE"
)# }
#> NSE:   0.9692195
#> mNSE:  0.8404272
#> KGE:   0.8993928
#> RMSE:  1.872744
#> PBIAS: 6.3
#> RSR:   0.1432493
#> RCV:   0.9220409
#> GMRAE: 0.1541012
#> MdRAE: 0.06190472
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.9171912
#> mNSE:  0.7903329
#> KGE:   0.8649407
#> RMSE:  3.071698
#> PBIAS: -11.6
#> RSR:   0.2349593
#> RCV:   0.9316332
#> GMRAE: 0.0867775
#> MdRAE: 0.2116875
#> 
#> In-channel retention ratio: 3.330669e-16
#> 

#> NSE:   0.8659469
#> mNSE:  0.640928
#> KGE:   0.7769389
#> RMSE:  3.908217
#> PBIAS: 20.5
#> RSR:   0.298946
#> RCV:   0.9127025
#> GMRAE: 0.3730323
#> MdRAE: 0.2531637
#> 
#> In-channel retention ratio: -8.881784e-16
#> 

#> NSE:   0.972409
#> mNSE:  0.8372972
#> KGE:   0.9246256
#> RMSE:  1.773064
#> PBIAS: -1.2
#> RSR:   0.1356246
#> RCV:   0.9263376
#> GMRAE: 0.1616979
#> MdRAE: 0.2276317
#> 
#> In-channel retention ratio: -6.661338e-16
#> 

#> NSE:   0.9653981
#> mNSE:  0.8360763
#> KGE:   0.9180164
#> RMSE:  1.985594
#> PBIAS: -3.7
#> RSR:   0.1518814
#> RCV:   0.9277124
#> GMRAE: 0.1410067
#> MdRAE: 0.2846086
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.975304
#> mNSE:  0.8384925
#> KGE:   0.9226028
#> RMSE:  1.677466
#> PBIAS: 1.5
#> RSR:   0.1283122
#> RCV:   0.9248487
#> GMRAE: 0.1723862
#> MdRAE: 0.1682604
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.972409
#> mNSE:  0.8372972
#> KGE:   0.9246256
#> RMSE:  1.773064
#> PBIAS: -1.2
#> RSR:   0.1356246
#> RCV:   0.9263376
#> GMRAE: 0.1616979
#> MdRAE: 0.2276317
#> 
#> In-channel retention ratio: -6.661338e-16
#> 

#> $maximum
#> [1] 0.00147171
#> 
#> $objective
#>       KGE 
#> 0.9246256 
#> 
```
