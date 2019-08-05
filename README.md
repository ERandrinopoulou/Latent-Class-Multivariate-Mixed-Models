# LatentClassMultivariateMixedModels

This repository includes functions to fit a multivariate latent class mixed model.

The function is called mv_lclme and it is included in the mv_lclme function.R file. The file build_JAGS_code.R includes a function that mv_lclme function uses in order to create the jags model.

The file other_functions.R includes further functions that are needed.

An example on the Pbc2 data can be found in the Analysis folder.

The arguments of the mv_lvlme function are:

- formulas: the formulas for the multivariate mixed model (as a list e.g.: 
                                   formulas <- list(serBilir ~ sex + year + (year | id),
                                                    spiders ~ sex + year + drug + (1 | id))
- data: data set
- classes: how many classes we want to assume
- families: the families of the outcomes
- hc: hierarchical centering (does not work yet!)
- RM_method: class selection approach (does not work yet!)
- predicted: whether the predicted values should be obtained in order to investigate the fit of the model
