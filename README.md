# sfcopula

The sfcopula package models spatial and spatio-temporal random fields with vine copulas. It is the successor of
[spcopula](https://github.com/bengraeler/spcopula) and handles data with [sf](https://r-spatial.github.io/sf/) and
[stars](https://r-spatial.github.io/stars/) instead of the former sp and spacetime packages. The starting point of 
this new package has been developed with the help of AI based on the previous implementation spcopula. This fresh 
start has also been used to streamline several 'historic' decisions looking at the complete package.

Vine copulas capture the dependence structure, and their bivariate building blocks depend on the distances
separating the locations. A correlogram based on Kendall's tau (comparable to the variogram in geostatistics/kriging)
models how the strength of dependence changes with distance. The package can estimate the dependence structure,
and interpolate and simulate the modelled random fields. It also calculates multivariate return periods based on
bivariate copulas or vine copulas.

## Data structures

| purpose | spcopula (sp/spacetime) | sfcopula (sf/stars) |
|---|---|---|
| spatial data | `SpatialPointsDataFrame` | `sf` with `POINT` geometries |
| locations without data | `SpatialPoints` | `sfc` (`POINT`) |
| spatio-temporal data | `STFDF` | `stars` vector data cube with one geometry and one time dimension |
| spatio-temporal targets | `STF`/`ST*` | `stars` cube, or `sf` with a time column (`timeCol`, or the time column of an `sftime` object) |

Use `as_spacetime_cube()` to turn a long `sf` table (one row per location and time) into a vector data cube.

## Example

```r
library(sf)
library(sfcopula)

data("meuse", package = "sfcopula")  # sf, EPSG:28992
meuse$marZinc <- plnorm(meuse$zinc, mean(log(meuse$zinc)), sd(log(meuse$zinc)))

bins <- calc_bins(meuse, var = "marZinc", nbins = 10, cutoff = 800)
calcKTau <- fit_cor_fun(bins, degree = 3)
neigh <- neighbours(meuse, var = "marZinc", size = 5)
```

The demos cover the full workflows, for example `demo("spatial_copula", package = "sfcopula")` or
`demo("spacetime_covar_vine_copula", package = "sfcopula")`.

## Migrating from spcopula

Functions, classes and arguments follow a snake_case naming scheme. The old spcopula names still work, but emit a
deprecation warning (see `?"sfcopula-deprecated"`): old function names point to the new (partly merged) functions,
and the new functions accept the old camelCase argument names (e.g. `spDepFun`, `coVar`, `min.dist`, `cor.method`).
Two groups keep the style of the copula package on purpose: the copula families (`asCopula()`, `cqsCopula()`,
`tawn3pCopula()`, `empiricalCopula()`, ...) and the copula helpers (`invdduCopula()`, `invddvCopula()`,
`qCopula_u()`, `rCopula_y()`, ...), as well as the methods for the copula package's generics (`dCopula()`,
`fitCopula()`, `dduCopula()`, ...). Slot names (e.g. `@topCop`, `@coVar`) are unchanged.

| spcopula | sfcopula |
|---|---|
| `spCopula()`, `spGeomCopula()` | `spatial_copula(..., combination = "convex" / "geometric")` |
| `stCopula()` | `spacetime_copula()` |
| `spVineCopula()`, `stVineCopula()` | `distance_vine_copula(trees, topCop = NULL)` |
| `stCoVarVineCopula()` | `covariate_vine_copula()` (spatial or spatio-temporal tree) |
| `neighbourhood()`, `stNeighbourhood()` | `neighbourhood()` |
| `getNeighbours()`, `getStNeighbours()` | `neighbours(data, target)` (dispatches on `sf` or `stars`) |
| `reduceNeighbours()` | `reduce_neighbours()` |
| `calcBins()`, `fitCorFun()` | `calc_bins()`, `fit_cor_fun()` |
| `loglikByCopulasLags()`, `loglikByCopulasStLags()` | `loglik_by_lags()` |
| `fitSpCopula()`, `composeSpCopula()` | `fit_spatial_copula()`, `compose_spatial_copula()` |
| `dropSpTree()`, `dropStTree()` | `drop_tree()` |
| `calcSpTreeDists()` | `tree_dists()` |
| `condSpVine()`, `condStVine()`, `condStCoVarVine()` | `cond_vine()` |
| `condCovariate()` | `cond_covariate()` |
| `spCopPredict()`, `stCopPredict()`, `spGaussCopPredict()` | `predict(model, neigh, data, target, margin)` |
| `spGaussLogLik()` | `spatial_gauss_loglik()` (models: `spatial_gauss_copula(cor_fun)`) |
| `kendallRP()`, `criticalLevel()`, `criticalPair()`, `criticalTriple()` | `kendall_rp()`, `critical_level()`, `critical_pair()`, `critical_triple()` |
| `kendallDistribution()`, `getKendallDistr()`, `genEmpKenFun()`, `genInvKenFun()` | `kendall_distribution()`, `get_kendall_distr()`, `gen_emp_ken_fun()`, `gen_inv_ken_fun()` |
| `bivJointDepFun()`, `lowerBivJointDepFun()`, `upperBivJointDepFun()` and the `emp...` versions | `biv_joint_dep_fun()`, `lower_biv_joint_dep_fun()`, `upper_biv_joint_dep_fun()`, `emp_biv_joint_dep_fun()`, ... |
| `rankTransform()`, `dependencePlot()`, `unitScatter()`/`univScatter()` | `rank_transform()`, `dependence_plot()`, `unit_scatter()` |

| spcopula argument | sfcopula argument |
|---|---|
| `spDepFun`, `stDepFun` | `dep_fun` |
| `coVar`, `coVarCop`, `topCop` | `covar`, `covar_cop`, `top_cop` |
| `min.dist`, `timeSteps`, `timeCol`, `spSize` | `min_dist`, `time_steps`, `time_col`, `size` |
| `cor.method`, `calcCor`, `lagSub`, `bestFit` | `cor_method`, `calc_cor`, `lag_sub`, `best_fit` |
| `condVar`, `stInd`, `n.trees`, `corFun`, `dataLocs` | `cond_var`, `centre`, `n_trees`, `cor_fun`, `data` |
| `stNeigh`, `dropEmpty`, `kendallFun`, `KRP`, `kenFun` | `neigh`, `drop_empty`, `kendall_fun`, `krp`, `ken_fun` |

The classes follow the same scheme: `spatial_copula`, `spacetime_copula`, `distance_vine_copula`,
`covariate_vine_copula`, `spatial_gauss_copula` and `neighbourhood`. Spatial and spatio-temporal variants share
one class or function, which dispatches on its input. Objects saved with spcopula can be translated with
`upgrade_spcopula()`, or a whole workspace with `load_spcopula("workspace.RData")`.

Beyond the renaming, only the spatial containers differ (see the table above), and some of these changes affect
results:

* Distances are Euclidean in CRS units for projected data. For geographic (lon/lat) coordinates they are great-circle
  distances in **metres** (`sf::st_distance`), whereas sp returned kilometres. Correlograms and spatial copulas fitted
  on lon/lat data therefore need distances in metres.
* Predictions are added as a column to the `sf` target, or as an attribute to the `stars` target.
* With identical random seeds, the results reproduce those of spcopula 0.2-5.

## Background

The package was initially developed in the DFG project "Developing Spatio-Temporal Copulas" at the Institute for
Geoinformatics, University of Münster. It grew along a series of use cases and should be seen as a proof of concept,
but it applies to many other use cases and has been used in a number of scientific publications.
