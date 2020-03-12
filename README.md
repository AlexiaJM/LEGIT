# LEGIT

[![](http://cranlogs.r-pkg.org/badges/LEGIT)](http://cran.rstudio.com/web/packages/LEGIT/index.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/LEGIT)](http://cran.rstudio.com/web/packages/LEGIT/index.html)

This is a R implementation of the *Latent Environmental &amp; Genetic InTeraction* (LEGIT) model. 

![](https://raw.githubusercontent.com/AlexiaJM/LEGIT/master/images/LEGIT.png)

## Citation

If you use this software, please cite:

``
Jolicoeur-Martineau, A., Wazana, A., Szekely, E., Steiner, M., Fleming, A. S., Kennedy, J. L., ... & Greenwood, C. M. (2018). Alternating optimization for G× E modelling with weighted genetic and environmental scores: Examples from the MAVAN study. Psychological methods.
``

If you use the function for GxE interaction testing, please also cite:

``
Jolicoeur-Martineau, A., Belsky, J., Szekely, E., Widaman, K., Pluess, M., Greenwood, C., & Wazana, A. (2020). Distinguishing differential susceptibility, diathesis-stress, and vantage sensitivity: Beyond the single gene and environment model. Development and Psychopathology, 32(1), 73-83. doi:10.1017/S0954579418001438
``

**Description**

The LEGIT model is an interaction model with two latent variables: a weighted sum of genetic variants (genetic score) and a weighted sum of environmental variables (environmental score). Alternating optimization is used to estimate the model parameters (https://arxiv.org/abs/1703.08111). This approach has greatly enhanced predictive power over traditional GxE models which include only a single genetic variant and a single environmental exposure. Although this approach was originally made for GxE modelling, it is flexible and does not require the use of genetic and environmental variables. It can also handle more than 2 latent variables (rather than just G and E) and 3-way interactions or more. The LEGIT model produces highly interpretable results and is very parameter-efficient thus it can even be used with small sample sizes (n < 250). Tools to determine the type of interaction (vantage sensitivity, diathesis-stress or differential susceptibility), with any number of genetic variants or environments, are available (https://psyarxiv.com/27uw8).

**How to use**

A vignette explaining how to use the software is available here : https://rawgit.com/AlexiaJM/LEGIT/master/inst/doc/LEGIT.html

An additional vignette explaining how it can be used for GxE testing as per [Belsky et al. (2013)](https://www.researchgate.net/publication/256600905_FormalGXEtestJCPP2013) is available here: https://rawgit.com/AlexiaJM/LEGIT/master/inst/doc/GxE_testing.html

**How to install**

To install the latest stable version, run in R :

* install.packages("LEGIT")

To install the latest GitHub development version which could contain new or experimental features, run in R :

* install.packages("devtools")
* devtools::install_github("AlexiaJM/LEGIT")

Regarding the GitHub installation. When it ask "Enter one or more numbers, or an empty line to skip updates:", just write nothing and press enter. Otherwise, you risk getting a "Cannot remove prior package X" error.

**Examples**

Here is an example with 2 latent variables and a 2-way interaction (see https://arxiv.org/abs/1703.08111) :
![](https://raw.githubusercontent.com/AlexiaJM/LEGIT/master/images/LEGIT_2way.png)

Here is an example with 3 latent variables and a 3-way interaction (see https://arxiv.org/abs/1703.08111) :
![](https://raw.githubusercontent.com/AlexiaJM/LEGIT/master/images/LEGIT_3way.png)

**References**

* https://arxiv.org/abs/1703.08111
* https://psyarxiv.com/27uw8
* https://ajolicoeur.wordpress.com/legit/
