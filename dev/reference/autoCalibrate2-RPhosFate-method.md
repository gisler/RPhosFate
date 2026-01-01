# Two dimensional automatic model calibration

Automatically calibrates the model with the help of a general-purpose
optimisation function. In contrast to
[`autoCalibrate`](https://gisler.github.io/RPhosFate/dev/reference/autoCalibrate-RPhosFate-method.md),
this method always utilises the overland and channel deposition rate at
the same time and never the respective enrichment ratio for calibration.
Beware of local optima and parameters approximately within the
convergence tolerance of interval end-points.

## Usage

``` r
# S4 method for class 'RPhosFate'
autoCalibrate2(
  x,
  substance,
  col,
  metric,
  method = "Nelder-Mead",
  lower = 0,
  upper = 0.1,
  control = list(fnscale = if (metric %in% c("NSE", "mNSE", "KGE")) -1 else 1)
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

- metric:

  A character string specifying the metric to optimise. See
  [`calibrationQuality`](https://gisler.github.io/RPhosFate/dev/reference/calibrationQuality-RPhosFate-method.md)
  for available metrics.

- method:

  A character string specifying the utilised optimisation method. See
  [`optim`](https://rdrr.io/r/stats/optim.html) for further information
  (use
  [`autoCalibrate`](https://gisler.github.io/RPhosFate/dev/reference/autoCalibrate-RPhosFate-method.md)
  instead of method `"Brent"`).

- lower:

  A numeric scalar or vector specifying the lower end-point(s) of the
  interval(s) to be searched.

- upper:

  A numeric scalar or vector specifying the upper end-point(s) of the
  interval(s) to be searched.

- control:

  A [`list`](https://rdrr.io/r/base/list.html) of control parameters
  passed on to [`optim`](https://rdrr.io/r/stats/optim.html). See
  [`optim`](https://rdrr.io/r/stats/optim.html) for further information.

## Value

An S4
[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/RPhosFate-class.md)
river catchment object and side effects in the form of raster files.

## See also

[`snapGauges`](https://gisler.github.io/RPhosFate/dev/reference/snapGauges-RPhosFate-method.md)

## Examples

``` r
# \donttest{
# temporary demonstration project copy
cv_dir <- demoProject()
#> Warning: A folder called "demoProject" already exists and is left as is.
# load temporary demonstration project
x <- RPhosFate(
  cv_dir = cv_dir,
  ls_ini = TRUE
)
# presupposed method calls
x <- firstRun(x, "SS")
x <- snapGauges(x)

x <- autoCalibrate2(
  x,
  "SS",
  col = "SS_load",
  metric = "KGE",
  method = "L-BFGS-B",
  lower = c(1e-3, 0),
  upper = c(2e-3, 2e-3),
  control = list(fnscale = -1, parscale = c(1e-3, 1e-3), factr = 1e12)
)# }
#> NSE:   0.6071161
#> mNSE:  0.4420663
#> KGE:   0.6753664
#> RMSE:  6.690713
#> PBIAS: -31.9
#> RSR:   0.5117838
#> RCV:   0.9389487
#> GMRAE: 0.3970503
#> MdRAE: 0.9200186
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.6071161
#> mNSE:  0.4420663
#> KGE:   0.6753664
#> RMSE:  6.690713
#> PBIAS: -31.9
#> RSR:   0.5117838
#> RCV:   0.9389487
#> GMRAE: 0.3970503
#> MdRAE: 0.9200186
#> 
#> In-channel retention ratio: -2.220446e-16
#> 
#> NSE:   0.6080637
#> mNSE:  0.4428216
#> KGE:   0.6757878
#> RMSE:  6.682639
#> PBIAS: -31.8
#> RSR:   0.5111662
#> RCV:   0.9389382
#> GMRAE: 0.3963558
#> MdRAE: 0.919039
#> 
#> In-channel retention ratio: 1.221245e-15
#> 

#> NSE:   0.6064526
#> mNSE:  0.4416167
#> KGE:   0.675093
#> RMSE:  6.69636
#> PBIAS: -31.9
#> RSR:   0.5122157
#> RCV:   0.9388427
#> GMRAE: 0.3972695
#> MdRAE: 0.9208093
#> 
#> In-channel retention ratio: 0.0005159496
#> 

#> NSE:   0.6071161
#> mNSE:  0.4420663
#> KGE:   0.6753664
#> RMSE:  6.690713
#> PBIAS: -31.9
#> RSR:   0.5117838
#> RCV:   0.9389487
#> GMRAE: 0.3970503
#> MdRAE: 0.9200186
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.9377421
#> mNSE:  0.8298658
#> KGE:   0.885968
#> RMSE:  2.663406
#> PBIAS: -9
#> RSR:   0.2037284
#> RCV:   0.9303839
#> GMRAE: 0.04548961
#> MdRAE: 0.01767349
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.9372657
#> mNSE:  0.8288419
#> KGE:   0.8854508
#> RMSE:  2.673578
#> PBIAS: -9
#> RSR:   0.2045064
#> RCV:   0.930417
#> GMRAE: 0.04929736
#> MdRAE: 0.02270224
#> 
#> In-channel retention ratio: 8.881784e-16
#> 

#> NSE:   0.9382161
#> mNSE:  0.8308907
#> KGE:   0.886484
#> RMSE:  2.653249
#> PBIAS: -8.9
#> RSR:   0.2029514
#> RCV:   0.9303506
#> GMRAE: 0.04080259
#> MdRAE: 0.01340401
#> 
#> In-channel retention ratio: 1.110223e-16
#> 

#> NSE:   0.9374155
#> mNSE:  0.8292696
#> KGE:   0.8856287
#> RMSE:  2.670384
#> PBIAS: -9
#> RSR:   0.2042621
#> RCV:   0.9302727
#> GMRAE: 0.04700307
#> MdRAE: 0.01944738
#> 
#> In-channel retention ratio: 0.0005129603
#> 

#> NSE:   0.9377421
#> mNSE:  0.8298658
#> KGE:   0.885968
#> RMSE:  2.663406
#> PBIAS: -9
#> RSR:   0.2037284
#> RCV:   0.9303839
#> GMRAE: 0.04548961
#> MdRAE: 0.01767349
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.4664724
#> mNSE:  0.2591888
#> KGE:   0.5642745
#> RMSE:  7.796841
#> PBIAS: 42.3
#> RSR:   0.5963934
#> RCV:   0.8958417
#> GMRAE: 0.7711987
#> MdRAE: 0.7283918
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.4700652
#> mNSE:  0.261723
#> KGE:   0.565708
#> RMSE:  7.770545
#> PBIAS: 42.2
#> RSR:   0.5943819
#> RCV:   0.8959632
#> GMRAE: 0.7685763
#> MdRAE: 0.7252692
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.4628606
#> mNSE:  0.2566501
#> KGE:   0.5628382
#> RMSE:  7.823188
#> PBIAS: 42.4
#> RSR:   0.5984087
#> RCV:   0.89572
#> GMRAE: 0.7738261
#> MdRAE: 0.7315197
#> 
#> In-channel retention ratio: -8.881784e-16
#> 

#> NSE:   0.467853
#> mNSE:  0.2601182
#> KGE:   0.5647587
#> RMSE:  7.786747
#> PBIAS: 42.2
#> RSR:   0.5956212
#> RCV:   0.8957155
#> GMRAE: 0.7704267
#> MdRAE: 0.726758
#> 
#> In-channel retention ratio: 0.0005136488
#> 

#> NSE:   0.4664724
#> mNSE:  0.2591888
#> KGE:   0.5642745
#> RMSE:  7.796841
#> PBIAS: 42.3
#> RSR:   0.5963934
#> RCV:   0.8958417
#> GMRAE: 0.7711987
#> MdRAE: 0.7283918
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.9707589
#> mNSE:  0.836942
#> KGE:   0.9234894
#> RMSE:  1.825311
#> PBIAS: -1.9
#> RSR:   0.1396211
#> RCV:   0.9267516
#> GMRAE: 0.1566302
#> MdRAE: 0.2445573
#> 
#> In-channel retention ratio: -8.881784e-16
#> 

#> NSE:   0.9705705
#> mNSE:  0.8369055
#> KGE:   0.9233325
#> RMSE:  1.831183
#> PBIAS: -2
#> RSR:   0.1400703
#> RCV:   0.9267934
#> GMRAE: 0.1560648
#> MdRAE: 0.2462799
#> 
#> In-channel retention ratio: 3.330669e-16
#> 

#> NSE:   0.9709438
#> mNSE:  0.8369785
#> KGE:   0.923639
#> RMSE:  1.819532
#> PBIAS: -1.8
#> RSR:   0.139179
#> RCV:   0.9267096
#> GMRAE: 0.1571871
#> MdRAE: 0.2428326
#> 
#> In-channel retention ratio: 4.440892e-16
#> 

#> NSE:   0.9705897
#> mNSE:  0.8366748
#> KGE:   0.9232838
#> RMSE:  1.830586
#> PBIAS: -2
#> RSR:   0.1400246
#> RCV:   0.9266386
#> GMRAE: 0.15667
#> MdRAE: 0.2456889
#> 
#> In-channel retention ratio: 0.0005125375
#> 

#> NSE:   0.9707589
#> mNSE:  0.836942
#> KGE:   0.9234894
#> RMSE:  1.825311
#> PBIAS: -1.9
#> RSR:   0.1396211
#> RCV:   0.9267516
#> GMRAE: 0.1566302
#> MdRAE: 0.2445573
#> 
#> In-channel retention ratio: -8.881784e-16
#> 

#> NSE:   0.9752522
#> mNSE:  0.8384166
#> KGE:   0.9230168
#> RMSE:  1.679224
#> PBIAS: 1.3
#> RSR:   0.1284466
#> RCV:   0.924948
#> GMRAE: 0.1720404
#> MdRAE: 0.1721504
#> 
#> In-channel retention ratio: 2.220446e-16
#> 

#> NSE:   0.9752215
#> mNSE:  0.8383812
#> KGE:   0.9231968
#> RMSE:  1.680266
#> PBIAS: 1.2
#> RSR:   0.1285264
#> RCV:   0.9249941
#> GMRAE: 0.1718613
#> MdRAE: 0.1739629
#> 
#> In-channel retention ratio: -4.440892e-16
#> 

#> NSE:   0.9752788
#> mNSE:  0.838452
#> KGE:   0.9228283
#> RMSE:  1.678322
#> PBIAS: 1.4
#> RSR:   0.1283777
#> RCV:   0.9249017
#> GMRAE: 0.1722083
#> MdRAE: 0.1703357
#> 
#> In-channel retention ratio: 1.110223e-16
#> 

#> NSE:   0.9751637
#> mNSE:  0.8381407
#> KGE:   0.9229651
#> RMSE:  1.682224
#> PBIAS: 1.3
#> RSR:   0.1286761
#> RCV:   0.9248341
#> GMRAE: 0.1722736
#> MdRAE: 0.1733189
#> 
#> In-channel retention ratio: 0.0005124075
#> 

#> NSE:   0.9752522
#> mNSE:  0.8384166
#> KGE:   0.9230168
#> RMSE:  1.679224
#> PBIAS: 1.3
#> RSR:   0.1284466
#> RCV:   0.924948
#> GMRAE: 0.1720404
#> MdRAE: 0.1721504
#> 
#> In-channel retention ratio: 2.220446e-16
#> 

#> NSE:   0.9736669
#> mNSE:  0.8376379
#> KGE:   0.925005
#> RMSE:  1.732172
#> PBIAS: -0.4
#> RSR:   0.1324967
#> RCV:   0.9259289
#> GMRAE: 0.1657955
#> MdRAE: 0.211108
#> 
#> In-channel retention ratio: -4.440892e-16
#> 

#> NSE:   0.9735493
#> mNSE:  0.8376018
#> KGE:   0.9249994
#> RMSE:  1.736038
#> PBIAS: -0.5
#> RSR:   0.1327924
#> RCV:   0.9259728
#> GMRAE: 0.1653983
#> MdRAE: 0.2128719
#> 
#> In-channel retention ratio: 0
#> 

#> NSE:   0.9737807
#> mNSE:  0.8376739
#> KGE:   0.9250022
#> RMSE:  1.728424
#> PBIAS: -0.3
#> RSR:   0.1322101
#> RCV:   0.925885
#> GMRAE: 0.1661836
#> MdRAE: 0.209342
#> 
#> In-channel retention ratio: -4.440892e-16
#> 

#> NSE:   0.9735343
#> mNSE:  0.8373667
#> KGE:   0.9248663
#> RMSE:  1.736529
#> PBIAS: -0.5
#> RSR:   0.13283
#> RCV:   0.9258156
#> GMRAE: 0.1659181
#> MdRAE: 0.2122566
#> 
#> In-channel retention ratio: 0.0005124726
#> 

#> NSE:   0.9736669
#> mNSE:  0.8376379
#> KGE:   0.925005
#> RMSE:  1.732172
#> PBIAS: -0.4
#> RSR:   0.1324967
#> RCV:   0.9259289
#> GMRAE: 0.1657955
#> MdRAE: 0.211108
#> 
#> In-channel retention ratio: -4.440892e-16
#> 

#> NSE:   0.9736871
#> mNSE:  0.8376442
#> KGE:   0.9250051
#> RMSE:  1.731508
#> PBIAS: -0.4
#> RSR:   0.1324459
#> RCV:   0.9259213
#> GMRAE: 0.1658641
#> MdRAE: 0.2107992
#> 
#> In-channel retention ratio: 2.220446e-16
#> 

#> NSE:   0.9735701
#> mNSE:  0.8376081
#> KGE:   0.925001
#> RMSE:  1.735353
#> PBIAS: -0.5
#> RSR:   0.13274
#> RCV:   0.9259651
#> GMRAE: 0.1654684
#> MdRAE: 0.2125634
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.9738003
#> mNSE:  0.8376803
#> KGE:   0.9250008
#> RMSE:  1.727781
#> PBIAS: -0.3
#> RSR:   0.1321608
#> RCV:   0.9258773
#> GMRAE: 0.1662506
#> MdRAE: 0.2090327
#> 
#> In-channel retention ratio: -2.220446e-16
#> 

#> NSE:   0.9735548
#> mNSE:  0.8373729
#> KGE:   0.9248671
#> RMSE:  1.735855
#> PBIAS: -0.4
#> RSR:   0.1327784
#> RCV:   0.9258079
#> GMRAE: 0.1659875
#> MdRAE: 0.2119479
#> 
#> In-channel retention ratio: 0.0005124721
#> 

#> NSE:   0.9736871
#> mNSE:  0.8376442
#> KGE:   0.9250051
#> RMSE:  1.731508
#> PBIAS: -0.4
#> RSR:   0.1324459
#> RCV:   0.9259213
#> GMRAE: 0.1658641
#> MdRAE: 0.2107992
#> 
#> In-channel retention ratio: 2.220446e-16
#> 

#> $par
#> [1] 0.001462119 0.000000000
#> 
#> $value
#> [1] 0.9250051
#> 
#> $counts
#> function gradient 
#>        7        7 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
#> 
```
