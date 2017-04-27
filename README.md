# LEGIT

[![](http://cranlogs.r-pkg.org/badges/LEGIT)](http://cran.rstudio.com/web/packages/LEGIT/index.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/LEGIT)](http://cran.rstudio.com/web/packages/LEGIT/index.html)

This is a R implementation of the Latent Environmental &amp; Genetic InTeraction (LEGIT) model. 

The LEGIT model is an interaction model with two latent variables: a weighted sum of genetic variants (genetic score) and a weighted sum of environmental variables (environmental score). Alternating optimization is used to estimate the model parameters (https://arxiv.org/abs/1703.08111). This approach has greatly enhanced predictive power over traditional  gene x environment (GxE) models which include only a single genetic variant and a single environmental exposure. Gene x gene interactions (GxG) and environment x environment interactions (ExE) can be incorporated into the genetic and environmental scores. Any standard model formula can be used, it doesn't have to be a two-way interaction model. Although this approach was made for GxE modelling, the genetic score doesn't need to contain genetic variables and the environmental score doesn't need to contain environmental variables, therefore this approach is highly flexible. 

The vignette is available here: https://rawgit.com/AlexiaJM/LEGIT/master/inst/doc/LEGIT.html

The paper discussing the alternating optimization approach to estimating the LEGIT parameters is available here: https://arxiv.org/abs/1703.08111

To install the latest GitHub development version, run in R :

install.packages("devtools")

devtools::install_github("AlexiaJM/LEGIT")
