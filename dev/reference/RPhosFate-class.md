# RPhosFate class

An S4 object representing a river catchment.

## Slots

- `cv_dir`:

  A character vector holding the project root (first position) and
  optionally the Monte Carlo input data directory (second position).

- `ls_ini`:

  A logical scalar specifying if the state of an existing project was
  loaded from disk.

- `is_ths`:

  An integer scalar holding the number of threads to use for processing,
  where applicable.

- `is_MCi`:

  An integer scalar holding the current Monte Carlo iteration if
  applicable.

- `cv_MCl`:

  A character vector holding the names of the layers, which shall be
  written to disk with the associated Monte Carlo iteration in their
  filenames upon calling the appropriate methods.

- `parameters`:

  An S4 object holding the model parameters.

- `topo`:

  An S4 object holding the raster layers related to topography in the
  broader sense.

- `erosion`:

  An S4 object holding the raster layers related to erosion.

- `transport`:

  An S4 object holding raster layers required for modelling transport.

- `substances`:

  An S4 object holding the substance raster layer containers.

- `helpers`:

  An S4 object holding helper data.

## See also

[`RPhosFate`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md),
[`catchment`](https://gisler.github.io/RPhosFate/dev/reference/catchment.md)
