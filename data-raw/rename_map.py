# old spcopula-style name -> new name (functions and S4 classes)
MAP = {
 "spCopula": "spatial_copula", "stCopula": "spacetime_copula",
 "spVineCopula": "spatial_vine_copula", "mixedSpVineCopula": "mixed_spatial_vine_copula",
 "pureSpVineCopula": "pure_spatial_vine_copula", "stVineCopula": "spacetime_vine_copula",
 "stCoVarVineCopula": "spacetime_covar_vine_copula", "spGeomCopula": "spatial_geom_copula",
 "stNeighbourhood": "spacetime_neighbourhood",
 "getNeighbours": "spatial_neighbours", "getStNeighbours": "spacetime_neighbours",
 "stCube": "as_spacetime_cube", "calcBins": "calc_bins", "fitCorFun": "fit_cor_fun",
 "spCopPredict": "spatial_predict", "stCopPredict": "spacetime_predict",
 "spGaussCopPredict": "spatial_gauss_predict", "spGaussLogLik": "spatial_gauss_loglik",
 "dropSpTree": "drop_spatial_tree", "dropStTree": "drop_spacetime_tree",
 "calcSpTreeDists": "spatial_tree_dists", "condSpVine": "cond_spatial_vine",
 "condStVine": "cond_spacetime_vine", "condStCoVarVine": "cond_spacetime_covar_vine",
 "reduceNeighbours": "reduce_neighbours", "loglikByCopulasLags": "loglik_by_lags",
 "loglikByCopulasStLags": "loglik_by_lags_spacetime", "fitSpCopula": "fit_spatial_copula",
 "composeSpCopula": "compose_spatial_copula", "condCovariate": "cond_covariate",
}
# "neighbourhood" is also an ordinary word and is handled separately
NEIGH_OLD, NEIGH_NEW = "neighbourhood", "spatial_neighbourhood"
