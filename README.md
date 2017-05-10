# LEGIT

[![](http://cranlogs.r-pkg.org/badges/LEGIT)](http://cran.rstudio.com/web/packages/LEGIT/index.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/LEGIT)](http://cran.rstudio.com/web/packages/LEGIT/index.html)

This is a R implementation of the *Latent Environmental &amp; Genetic InTeraction* (LEGIT) model. 

![](https://raw.githubusercontent.com/AlexiaJM/LEGIT/master/images/LEGIT.png)

**Description**

The LEGIT model is an interaction model with two latent variables: a weighted sum of genetic variants (genetic score) and a weighted sum of environmental variables (environmental score). Alternating optimization is used to estimate the model parameters (https://arxiv.org/abs/1703.08111). This approach has greatly enhanced predictive power over traditional GxE models which include only a single genetic variant and a single environmental exposure. Although this approach was originally made for GxE modelling, it is flexible and does not require the use of genetic and environmental variables. It can also handle more than 2 latent variables (rather than just G and E) and 3-way interactions or more. The LEGIT model produces highly interpretable results and is very parameter-efficient thus it can even be used with small sample sizes (n < 250).

**How to use**

A vignette explaining how to use the software is available here : https://rawgit.com/AlexiaJM/LEGIT/master/inst/doc/LEGIT.html

**How to install**

To install the latest stable version, run in R :

* install.packages("LEGIT")

To install the latest GitHub development version, run in R :

* install.packages("devtools")
* devtools::install_github("AlexiaJM/LEGIT")

**Example**

Here is an example with 2 latent variables and a 2-way interaction (see: https://arxiv.org/abs/1703.08111) :
![](https://raw.githubusercontent.com/AlexiaJM/LEGIT/master/images/LEGIT_2way.png)

Here is an example with 3 latent variables and a 3-way interaction (see: https://arxiv.org/abs/1703.08111) :
![](https://raw.githubusercontent.com/AlexiaJM/LEGIT/master/images/LEGIT_3way.png)

**References**

* https://arxiv.org/abs/1703.08111
* https://ajolicoeur.wordpress.com/legit/
