# sfcopula 0.1-0

* First release of sfcopula, the successor of spcopula 0.2-5. sf and stars replace sp and spacetime.
* All functions and arguments use snake_case names, except the copula families and copula helpers, which keep the
  copula package's style on purpose (see README for the mapping). Deprecated: the spcopula function names (aliases)
  and the spcopula camelCase argument names (accepted with a warning).
* Functions and S4 classes use snake_case names, and spatial and spatio-temporal variants are merged. Each merged
  function or class dispatches on its input (see README for the full mapping):
  * `spatial_copula()` covers convex and geometric combinations (`combination` argument).
  * `distance_vine_copula()` covers spatial and spatio-temporal vines, with an optional top copula.
    Spatio-temporal vines can now have several trees.
  * The `neighbourhood` class and `neighbours()` cover `sf` and `stars` data.
  * `drop_tree()`, `tree_dists()`, `cond_vine()` and `loglik_by_lags()` cover both kinds of data.
  * `predict()` methods for distance vines, the new `spatial_gauss_copula` class and covariate vines (new) replace the separate prediction functions and share one implementation.
  The old spcopula function names are deprecated aliases for the merged functions.
* `covariate_vine_copula(coVarCop, tree, topCop)` replaces spcopula's spatio-temporal covariate vine and also accepts
  a spatial tree, so spatial data with a covariate can be modelled as well. New `fitCopula()` method: drops the tree,
  conditions the covariate and fits the top copula in one call. `predict()`, `cond_vine()` and `dCopula()` work for
  both kinds of tree; the index of the central location is passed as `centre`.
* Spatial neighbourhoods support covariates: `neighbours(..., coVar = )` adds the covariate at the central location,
  taken from the data for fitting and from the target for prediction (as for spatio-temporal data).
  `cond_covariate()` works for spatial and spatio-temporal neighbourhoods.
* `cond_vine()` for covariate vines expects the conditioning values in neighbourhood order (neighbours, then
  covariate); the deprecated `condStCoVarVine()` keeps spcopula's order (covariate first).
* New functions `upgrade_spcopula()` and `load_spcopula()` translate objects saved with spcopula into the new
  classes.
* Spatial data are `sf`/`sfc` POINT objects. Spatio-temporal data are `stars` vector data cubes (geometry x time).
* Spatio-temporal prediction targets can be `stars` cubes or `sf` objects with a time column (new argument `timeCol`
  in `neighbours()`). Covariates for prediction are taken from the target.
* New function `as_spacetime_cube()` converts long `sf` tables into vector data cubes.
* Geographic coordinates use great-circle distances in metres (sp: kilometres).
* Nearest neighbour search is implemented in vectorised R, so the package no longer contains compiled code.
* Data sets `EU_RB` and `EU_RB_2005` are `stars` cubes. `meuse` is shipped as an `sf` object.
* Performance: bivariate copulas of the families known to VineCopula (Gaussian, Student t, Clayton, Gumbel, Frank,
  Joe and the VC2copula families) are evaluated vectorised in compiled code, with parameters calibrated once per
  distinct distance; other families use the copula package as before. Inverse partial derivatives use a vectorised
  bisection (accurate to 1e-12 instead of about 1e-4), and vine simulation draws all samples at once. Measured on the
  meuse and EU_RB data: spatial copulas with a dependence function 27x, fitting distance vines 3x, prediction 2.4x,
  `invdduCopula()`/`invddvCopula()` 16x, `rCopula()` for vines 9x, `calc_bins()` on cubes 5x and
  `loglik_by_lags()` on cubes 8x faster; results agree with spcopula to 1e-8.
* Where the dependence function gives independence, all families now have exactly the same log-likelihood (spcopula
  differed by rounding noise of about 1e-11). Select the best family per lag with
  `apply(loglik, 1, which.max)`; the idiom `which(rank(x) == length(families))` fails on such ties.
* Bug fixes over spcopula:
  * The `show()` method of `spatial_copula` was overwritten by the one for geometric copulas, so every spatial
    copula printed "Spatial Copula based on geometric means".
  * `spatial_copula()`/`spGeomCopula()` failed with their default `unit = NULL` (now `"m"`).
  * Dropping a spatio-temporal tree (`dropStTree()`, now `drop_tree()`) returned the distances between the centre
    and each neighbour instead of between the remaining neighbours, so a second spatio-temporal tree would have used
    wrong distances. The conditioned data were correct.
  * Simulating from spatial vines with several trees used the first tree's copula when updating the
    pseudo-observations of later trees.
  * Fitting a spatial vine whose trees left one or two variables at the top failed (broken `if`/`else` chain and
    wrong argument order in the log-likelihood).
  * `dCopula(..., log = TRUE)` for geometric spatial copulas returned the density instead of its logarithm.
  * `dCopula()` for spatio-temporal covariate vines could not be evaluated (it required two columns and passed its
    arguments in the wrong order); it now evaluates the full density with the columns ordered as in a neighbourhood.
  * Spatial neighbourhoods with a covariate failed the validity check, as the covariate values were never added.
  * Further: `stCopPredict(method = "expectation")` indexed the distance array incorrectly;
  `fitCopula()` for `stVineCopula` failed when `method` was a list; `calcBins()` for spatio-temporal data counted pairs
  with missing values for `cor.method` other than `"fasttau"` and failed for vectors of `instances`; `spCopula()` failed for more than one independence
  copula among the components; the `neighbourhood` validity check used an undefined variable; and the `stVineCopFit` demo used outdated arguments.
