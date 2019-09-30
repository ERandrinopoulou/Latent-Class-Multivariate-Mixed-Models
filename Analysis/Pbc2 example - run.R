library(JMbayes)
library(jagsUI)
library(label.switching)

formulas <- list(serBilir ~ sex + ns(year, 2) + (ns(year, 2) | id),
                 spiders ~ sex + year + drug + (1 | id))
data <- pbc2

families <- c("gaussian", "binomial")

classes <- 1
hc <- FALSE
RM_method <- FALSE
predicted <- TRUE

source("build_JAGS_code.R")
source("other_functions.R")
source("mv_lclme function.R")


fitLCLME <- mv_lclme(formulas, data, classes, families, hc, RM_method, predicted)

summary.mvlclme(fitLCLME, classes)
