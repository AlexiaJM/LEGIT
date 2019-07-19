#' @title Simulated example of a 2 way interaction GxE model.
#' @description Simulated example of a 2 way interaction GxE model (where G and E are latent variables). 
#' \deqn{g_j \sim Binomial(n=1,p=.30)}
#' \deqn{j = 1, 2, 3, 4}
#' \deqn{e_l \sim Normal(\mu=0,\sigma=1.5)}
#' \deqn{l = 1, 2, 3}
#' \deqn{g = .2g_1 + .15g_2 - .3g_3 + .1g_4 + .05g_1g_3 + .2g_2g_3}
#' \deqn{e = -.45e_1 + .35e_2 + .2e_3}
#' \deqn{\mu = -1 + 2g + 3e + 4ge}
#' \tabular{cc}{
#' \eqn{y \sim Normal(\mu=\mu,\sigma=\code{sigma})} if \code{logit}=FALSE \cr
#' \eqn{y \sim Binomial(n=1,p=logit(\mu))} if \code{logit}=TRUE
#' }
#' @param N Sample size.
#' @param sigma Standard deviation of the gaussian noise (if \code{logit}=FALSE).
#' @param logit If TRUE, the outcome is transformed to binary with a logit link.
#' @param seed RNG seed.
#' @return Returns a list containing, in the following order: data.frame with the observed outcome (with noise) and the true outcome (without noise), data.frame of the genetic variants (G), data.frame of the environments (E), vector of the true genetic coefficients, vector of the true environmental coefficients, vector of the true main model coefficients
#' @examples
#'	example_2way(5,1,logit=FALSE)
#'	example_2way(5,0,logit=TRUE)
#' @export
"example_2way"

#' @title Simulated example of a 2 way interaction GxE model with crossover point.
#' @description Simulated example of a 2 way interaction GxE model with crossover point (where G and E are latent variables). 
#' \deqn{g_j \sim Binomial(n=1,p=.30)}
#' \deqn{j = 1, 2, 3, 4}
#' \deqn{e_l \sim 10 Beta(\alpha,\beta))}
#' \deqn{l = 1, 2, 3}
#' \deqn{g = .30g_1 + .10g_2 + .20g_3 + .40g_4}
#' \deqn{e = .45e_1 + .35e_2 + .2e_3}
#' \deqn{\mu = coef[1] + coef[2]e + coef[3]ge}
#' \tabular{cc}{
#' \eqn{y \sim Normal(\mu=\mu,\sigma=\code{sigma})} if \code{logit}=FALSE \cr
#' \eqn{y \sim Binomial(n=1,p=logit(\mu))} if \code{logit}=TRUE
#' }
#' @param N Sample size.
#' @param sigma Standard deviation of the gaussian noise (if \code{logit}=FALSE).
#' @param c crossover point
#' @param coef_main Coefficients of the main model, must be a vector of size 3 for intercept, E main effect and GxE effect (Default = c(0,1,2)).
#' @param coef_G Coefficients of the 4 genes, must be a vector of size 4 (Default = c(.30, .10, .20, .40)).
#' @param coef_E Coefficients of the 3 environments, must be a vector of size 3 (Default = c(.45, .35, .2)).
#' @param logit If TRUE, the outcome is transformed to binary with a logit link.
#' @param seed RNG seed.
#' @param beta_param Vector of size two for the parameters of the beta distribution of the environmental variables (Default = c(2,2)).
#' @return Returns a list containing, in the following order: data.frame with the observed outcome (with noise) and the true outcome (without noise), data.frame of the genetic variants (G), data.frame of the environments (E), vector of the true genetic coefficients, vector of the true environmental coefficients, vector of the true main model coefficients, the crossover point.
#' @examples
#' ## Examples
#' # Diathesis Stress WEAK
#' ex_dia = example_with_crossover(250, c=10, coef_main = c(3,1,2), sigma=1)
#' # Diathesis Stress STRONG
#' ex_dia_s = example_with_crossover(250, c=10, coef_main = c(3,0,2), sigma=1)
#' # Differential Susceptibility WEAK
#' ex_ds = example_with_crossover(250, c=5, coef_main = c(3+5,1,2), sigma=1)
#' # Differential Susceptibility STRONG
#' ex_ds_s = example_with_crossover(250, c=5, coef_main = c(3+5,0,2), sigma=1)
#' @export
"example_with_crossover"

#' @title Simulated example of a 3 way interaction GxExz model
#' @description Simulated example of a 3 way interaction GxExz model (where G and E are latent variables). 
#' \deqn{g_j \sim Binomial(n=1,p=.30)}
#' \deqn{j = 1, 2, 3, 4}
#' \deqn{e_l \sim Normal(\mu=0,\sigma=1.5)}
#' \deqn{l = 1, 2, 3}
#' \deqn{z \sim Normal(\mu=3,\sigma=1)}
#' \deqn{g = .2g_1 + .15g_2 - .3g_3 + .1g_4 + .05g_1g_3 + .2g_2g_3}
#' \deqn{e = -.45e_1 + .35e_2 + .2e_3}
#' \deqn{\mu = -2 + 2g + 3e + z + 5ge - 1.5ez + 2gz + 2gez}
#' \tabular{cc}{
#' \eqn{y \sim Normal(\mu=\mu,\sigma=\code{sigma})} if \code{logit}=FALSE \cr
#' \eqn{y \sim Binomial(n=1,p=logit(\mu))} if \code{logit}=TRUE
#' }
#' @param N Sample size.
#' @param sigma Standard deviation of the gaussian noise (if \code{logit}=FALSE).
#' @param logit If TRUE, the outcome is transformed to binary with a logit link.
#' @param seed RNG seed.
#' @return Returns a list containing, in the following order: data.frame with the observed outcome (with noise), the true outcome (without noise) and \eqn{z}, data.frame of the genetic variants (G), data.frame of the environments (E), vector of the true genetic coefficients, vector of the true environmental coefficients, vector of the true main model coefficients
#' @examples
#'	example_3way(5,2.5,logit=FALSE)
#'	example_3way(5,0,logit=TRUE)
#' @export
"example_3way"

#' @title Simulated example of a 3 way interaction GxExZ model
#' @description Simulated example of a 3 way interaction GxExZ model (where G, E and Z are latent variables). 
#' \deqn{g_j \sim Binomial(n=1,p=.30)}
#' \deqn{j = 1, 2, 3, 4}
#' \deqn{e_k \sim Normal(\mu=0,\sigma=1.5)}
#' \deqn{k = 1, 2, 3}
#' \deqn{z_l \sim Normal(\mu=3,\sigma=1)}
#' \deqn{l = 1, 2, 3}
#' \deqn{g = .2g_1 + .15g_2 - .3g_3 + .1g_4 + .05g_1g_3 + .2g_2g_3}
#' \deqn{e = -.45e_1 + .35e_2 + .2e_3}
#' \deqn{z = .15z_1 + .60z_2 + .25z_3}
#' \deqn{\mu = -2 + 2g + 3e + z + 5ge - 1.5ez + 2gz + 2gez}
#' \tabular{cc}{
#' \eqn{y \sim Normal(\mu=\mu,\sigma=\code{sigma})} if \code{logit}=FALSE \cr
#' \eqn{y \sim Binomial(n=1,p=logit(\mu))} if \code{logit}=TRUE
#' }
#' @param N Sample size.
#' @param sigma Standard deviation of the gaussian noise (if \code{logit}=FALSE).
#' @param logit If TRUE, the outcome is transformed to binary with a logit link.
#' @param seed RNG seed.
#' @return Returns a list containing, in the following order: data.frame with the observed outcome (with noise) and the true outcome (without noise), list containing the data.frame of the genetic variants (G), the data.frame of the \eqn{e} environments (E) and the data.frame of the \eqn{z} environments (Z), vector of the true genetic coefficients, vector of the true \eqn{e} environmental coefficients, vector of the true \eqn{z} environmental coefficients, vector of the true main model coefficients
#' @examples
#'	example_3way_3latent(5,1,logit=FALSE)
#'	example_3way_3latent(5,0,logit=TRUE)
#' @export
"example_3way_3latent"

#' @title Longitudinal folds
#' @description Function to create folds adequately for longitudinal datasets by forcing every observation with the same id to be in the same fold. Can be used with LEGIT_cv to make sure that the cross-validation folds are appropriate when using longitudinal data.
#' @param cv_iter Number of cross-validation iterations (Default = 1).
#' @param cv_folds Number of cross-validation folds (Default = 10).
#' @param id Factor vector containing the id number of each observation.
#' @param formula Optional Model formula. If data and formula are provided, only the non-missing observations will be used when creating the folds (Put "formula" here if you have missing data).
#' @param data Optional data.frame used for the formula. If data and formula are provided, only the non-missing observations will be used when creating the folds (Put "data" here if you have missing data).
#' @param data_needed Optional data.frame with variables that have to be included (Put "cbind(genes,env)"" or "latent_var" here if you have missing data).
#' @param print If FALSE, nothing except warnings will be printed. (Default = TRUE).
#' @return Returns a list of vectors containing the fold number for each observation
#' @examples
#'	train = example_2way(500, 1, seed=777)
#'	# Assuming it's longitudinal with 4 timepoints, even though it's not
#'	id = factor(rep(1:125,each=4))
#'	fit_cv = LEGIT_cv(train$data, train$G, train$E, y ~ G*E, folds=longitudinal_folds(1,10, id))
#' @export
"longitudinal_folds"

#' @title Latent Environmental & Genetic InTeraction (LEGIT) model
#' @description Constructs a generalized linear model (glm) with a weighted latent environmental score and weighted latent genetic score using alternating optimization.
#' @param data data.frame of the dataset to be used. 
#' @param genes data.frame of the variables inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic).
#' @param env data.frame of the variables inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental).
#' @param formula Model formula. Use \emph{E} for the environmental score and \emph{G} for the genetic score. Do not manually code interactions, write them in the formula instead (ex: G*E*z or G:E:z).
#' @param start_genes Optional starting points for genetic score (must be the same length as the number of columns of \code{genes}).
#' @param start_env Optional starting points for environmental score (must be the same length as the number of columns of \code{env}).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param print If FALSE, nothing except warnings will be printed (Default = TRUE).
#' @param print_steps If TRUE, print the parameters at all iterations, good for debugging (Default = FALSE).
#' @param crossover If not NULL, estimates the crossover point of \emph{E} using the provided value as starting point (To test for diathesis-stress vs differential susceptibility).
#' @param crossover_fixed If TRUE, instead of estimating the crossover point of E, we force/fix it to the value of "crossover". (Used when creating a diathes-stress model) (Default = FALSE).
#' @param reverse_code If TRUE, after fitting the model, the genes with negative weights are reverse coded (ex: \eqn{g_rev} = 1 - \eqn{g}). It assumes that the original coding is in [0,1]. The purpose of this option is to prevent genes with negative weights which cause interpretation problems (ex: depression normally decreases attention but with a negative genetic score, it increases attention). Warning, using this option with GxG interactions could cause nonsensical results since GxG could be inverted. Also note that this may fail with certain models (Default=FALSE).
#' @param rescale If TRUE, the environmental variables are automatically rescaled to the range [-1,1]. This improves interpretability (Default=FALSE).
#' @return Returns an object of the class "LEGIT" which is list containing, in the following order: a glm fit of the main model, a glm fit of the genetic score, a glm fit of the environmental score, a list of the true model parameters (AIC, BIC, rank, df.residual, null.deviance) for which the individual model parts (main, genetic, environmental) don't estimate properly and the formula.
#' @examples
#'	train = example_2way(500, 1, seed=777)
#'	fit_best = LEGIT(train$data, train$G, train$E, y ~ G*E, train$coef_G, train$coef_E)
#'	fit_default = LEGIT(train$data, train$G, train$E, y ~ G*E)
#'	summary(fit_default)
#'	summary(fit_best)
#'	train = example_3way(500, 2.5, seed=777)
#'	fit_best = LEGIT(train$data, train$G, train$E, y ~ G*E*z, train$coef_G, train$coef_E)
#'	fit_default = LEGIT(train$data, train$G, train$E, y ~ G*E*z)
#'	summary(fit_default)
#'	summary(fit_best)
#' @import formula.tools stats
#' @references Alexia Jolicoeur-Martineau, Ashley Wazana, Eszter Szekely, Meir Steiner, Alison S. Fleming, James L. Kennedy, Michael J. Meaney, Celia M.T. Greenwood and the MAVAN team. \emph{Alternating optimization for GxE modelling with weighted genetic and environmental scores: examples from the MAVAN study} (2017). arXiv:1703.08111.
#' @export
"LEGIT"

#' @title Testing of the GxE interaction
#' @description Testing of the GxE interaction using the competitive-confirmatory approach adapted from Belsky, Pluess et Widaman (2013). Reports the different hypotheses (diathesis-stress, vantage-sensitivity, or differential susceptibility), assuming or not assuming a main effect for \emph{E} (WEAK vs STRONG) using the LEGIT model.
#' @param data data.frame of the dataset to be used. 
#' @param genes data.frame of the variables inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic).
#' @param env data.frame of the variables inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental).
#' @param formula_noGxE formula WITHOUT \emph{G} or \emph{E} (y ~ covariates). \emph{G} and \emph{E} will automatically be added properly based on the hypotheses tested.
#' @param crossover A tuple containting the minimum and maximum of the environment used as crossover point of \emph{E} used in the vantage sensitivity and diathesis-stress models. Instead of providing two number, you can also write c("min","max") to automatically choose the expected minimum or maximum of the environmental score which is calculated based on the min/max of the environments and the current weights.
#' @param include_noGxE_models If True, we test for models with only G, only E, both G and E, neither G and E (four models without a GxE). This is to verify for false positives, if one of those models has the best fit, then it is possible that there is no GxE, thus no type of GxE. With a single gene and environment, simply looking at the p-value of the GxE is good enough to get around 5-10 percent false positive rate, but with multiple genes and environments, we need to compare model fits to get a low false positive rate. Use your own judgment when using this because if you have multiple genes and environments and small/moderate N, a model without GxE could have a lower BIC but still not be the actual best model. However, if you see little difference in BIC between all 4 GxE models and the non-GxE models have much lower BIC, than it is likely that there is no GxE. Note that this is only implemented for AIC, AICc and BIC. (Default = True)
#' @param reverse_code If TRUE, after fitting the model, the genes with negative weights are reverse coded (ex: \eqn{g_{rev}} = 1 - \eqn{g}). It assumes that the original coding is in [0,1]. The purpose of this option is to prevent genes with negative weights which cause interpretation problems (ex: depression normally decreases attention but with a negative genetic score, it increases attention). Warning, using this option with GxG interactions could cause nonsensical results since GxG could be inverted. Also note that this may fail with certain models (Default=FALSE).
#' @param rescale If TRUE, the environmental variables are automatically rescaled to the range [-1,1]. This improves interpretability (Default=FALSE).
#' @param boot Optional number of bootstrap samples. If not NULL, we use bootstrap to find the confidence interval of the crossover point. This provides more realistic confidence intervals. Make sure to use a bigger number (>= 1000) to get good precision; also note that a too small number could return an error ("estimated adjustment 'a' is NA").
#' @param criterion Criterion used to assess which model is the best. It can be set to "AIC", "AICc", "BIC", "cv", "cv_AUC", "cv_Huber" (Default="BIC").
#' @param start_genes Optional starting points for genetic score (must be the same length as the number of columns of \code{genes}).
#' @param start_env Optional starting points for environmental score (must be the same length as the number of columns of \code{env}).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param folds Optional list of vectors containing the fold number for each observation. Bypass cv_iter and cv_folds. Setting your own folds could be important for certain data types like time series or longitudinal data.
#' @param id Optional id of observations, can be a vector or data.frame (only used when returning list of possible outliers).
#' @param classification Set to TRUE if you are doing classification (binary outcome).
#' @param seed Seed for cross-validation folds.
#' @param Huber_p Parameter controlling the Huber cross-validation error (Default = 1.345).
#' @return Returns a list containing 1) the six models ordered from best to worse (vantage sensitivity WEAK/STRONG, diathesis-stress WEAK/STRONG, differential susceptibility WEAK/STRONG) and 2) a data frame with the criterion, the crossover, 95\% coverage of the crossover, whether the crossover 95\% interval is within the observable range and the percentage of observations below the crossover point in order from best to worst based on the selected criterion. Models not within the observable range should be rejected even if the criterion is slightly better. An extremely low percentage of observations below the crossover point is also evidence toward diathesis-stress. Note that we assume that the environmental score is from bad to good but if this is not the case, then the models labelled as "diathesis-stress" could actually reflect vantage sensitivity and vice-versa. If outcome is Good-to-Bad: C=min(E) is diathesis-stress, C=max(E) is vantage sensitivity. If outcome is Bad-to-Good: C=max(E) is diathesis-stress, C=min(E) is vantage sensitivity.
#' @examples
#' \dontrun{
#' ## Examples where x is in [0, 10]
#' # Diathesis Stress WEAK
#' ex_dia = example_with_crossover(250, c=10, coef_main = c(3,1,2), sigma=1)
#' # Diathesis Stress STRONG
#' ex_dia_s = example_with_crossover(250, c=10, coef_main = c(3,0,2), sigma=1)
#' ## Assuming there is a crossover point at x=5
#' # Differential Susceptibility WEAK
#' ex_ds = example_with_crossover(250, c=5, coef_main = c(3+5,1,2), sigma=1)
#' # Differential Susceptibility STRONG
#' ex_ds_s = example_with_crossover(250, c=5, coef_main = c(3+5,0,2), sigma=1)
#' 
#' ## If true model is "Diathesis Stress WEAK"
#' GxE_test_BIC = GxE_interaction_test(ex_dia$data, ex_dia$G, ex_dia$E, 
#' formula_noGxE = y ~ 1, start_genes = ex_dia$coef_G, start_env = ex_dia$coef_E, 
#' criterion="BIC")
#' GxE_test_BIC$results
#' 
#' ## If true model is "Diathesis Stress STRONG"
#' GxE_test_BIC = GxE_interaction_test(ex_dia_s$data, ex_dia_s$G, ex_dia_s$E, 
#' formula_noGxE = y ~ 1, start_genes = ex_dia_s$coef_G, start_env = ex_dia_s$coef_E, 
#' criterion="BIC")
#' GxE_test_BIC$results
#' 
#' ## If true model is "Differential susceptibility WEAK"
#' GxE_test_BIC = GxE_interaction_test(ex_ds$data, ex_ds$G, ex_ds$E, 
#' formula_noGxE = y ~ 1, start_genes = ex_ds$coef_G, start_env = ex_ds$coef_E, 
#' criterion="BIC")
#' GxE_test_BIC$results
#' 
#' ## If true model is "Differential susceptibility STRONG"
#' GxE_test_BIC = GxE_interaction_test(ex_ds_s$data, ex_ds_s$G, ex_ds_s$E, 
#' formula_noGxE = y ~ 1, start_genes = ex_ds_s$coef_G, start_env = ex_ds_s$coef_E,
#' criterion="BIC")
#' GxE_test_BIC$results
#' 
#' # Example of plots
#' plot(GxE_test_BIC$fits$diff_suscept_STRONG, xlim=c(0,10), ylim=c(3,13))
#' plot(GxE_test_BIC$fits$diff_suscept_WEAK, xlim=c(0,10), ylim=c(3,13))
#' plot(GxE_test_BIC$fits$diathesis_stress_STRONG, xlim=c(0,10), ylim=c(3,13))
#' plot(GxE_test_BIC$fits$diathesis_stress_WEAK, xlim=c(0,10), ylim=c(3,13))
#' }
#' @import formula.tools stats
#' @references Alexia Jolicoeur-Martineau, Jay Belsky, Eszter Szekely, Keith F. Widaman, Michael Pluess, Celia Greenwood and Ashley Wazana. \emph{Distinguishing differential susceptibility, diathesis-stress and vantage sensitivity: beyond the single gene and environment model} (2017). psyarxiv.com/27uw8. 10.17605/OSF.IO/27UW8.
#' @references Alexia Jolicoeur-Martineau, Ashley Wazana, Eszter Szekely, Meir Steiner, Alison S. Fleming, James L. Kennedy, Michael J. Meaney, Celia M.T. Greenwood and the MAVAN team. \emph{Alternating optimization for GxE modelling with weighted genetic and environmental scores: examples from the MAVAN study} (2017). arXiv:1703.08111.
#' @references Jay Belsky, Michael Pluess and Keith F. Widaman. \emph{Confirmatory and competitive evaluation of alternative gene-environment interaction hypotheses} (2013). Journal of Child Psychology and Psychiatry, 54(10), 1135-1143.
#' @export
"GxE_interaction_test"

#' @title Regions of significance using Johnson-Neyman technique
#' @description Constructs a LEGIT model and returns the regions of significance (RoS) with the predicted type of interaction (diathesis-stress, vantage-sensitivity, or differential susceptibility). RoS is not recommended due to poor accuracy with small samples and small effect sizes, GxE_interaction_test has much better accuracy overall. Only implemented for family=gaussian.
#' @param data data.frame of the dataset to be used. 
#' @param genes data.frame of the variables inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic).
#' @param env data.frame of the variables inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental).
#' @param formula_noGxE formula WITHOUT \emph{G} or \emph{E} (y ~ covariates). \emph{G} and \emph{E} will automatically be added.
#' @param t_alpha Alpha level of the student-t distribution for the regions of significance (Default = .05)
#' @param start_genes Optional starting points for genetic score (must be the same length as the number of columns of \code{genes}).
#' @param start_env Optional starting points for environmental score (must be the same length as the number of columns of \code{env}).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param reverse_code If TRUE, after fitting the model, the genes with negative weights are reverse coded (ex: \eqn{g_rev} = 1 - \eqn{g}). It assumes that the original coding is in [0,1]. The purpose of this option is to prevent genes with negative weights which cause interpretation problems (ex: depression normally decreases attention but with a negative genetic score, it increases attention). Warning, using this option with GxG interactions could cause nonsensical results since GxG could be inverted. Also note that this may fail with certain models (Default=FALSE).
#' @param rescale If TRUE, the environmental variables are automatically rescaled to the range [-1,1]. This improves interpretability (Default=FALSE).
#' @return Returns a list containing the RoS and the predicted type of interaction.
#' @examples
#'	train = example_2way(500, 1, seed=777)
#'	ros = GxE_interaction_RoS(train$data, train$G, train$E, y ~ 1)
#'	ros
#' @import formula.tools stats
#' @references Alexia Jolicoeur-Martineau, Jay Belsky, Eszter Szekely, Keith F. Widaman, Michael Pluess, Celia Greenwood and Ashley Wazana. \emph{Distinguishing differential susceptibility, diathesis-stress and vantage sensitivity: beyond the single gene and environment model} (2017). psyarxiv.com/27uw8. 10.17605/OSF.IO/27UW8.
#' @references  Daniel J. Bauer & Patrick J. Curran. \emph{Probing Interactions in Fixed and Multilevel Regression: Inferential and Graphical Techniques} (2005). Multivariate Behavioral Research, 40:3, 373-400, DOI: 10.1207/s15327906mbr4003_5.
#' @export
"GxE_interaction_RoS"

#' @title Plot
#' @description Plot of LEGIT models. By default, variables that are not in \emph{G} or \emph{E} are fixed to the mean.
#' @param x An object of class "LEGIT", usually, a result of a call to LEGIT.
#' @param cov_values Vector of the values, for each covariate, that will be used in the plotting, if there are any covariates. It must contain the names of the variables. Covariates are the variables that are not \emph{G} nor \emph{E} but still are adjusted for in the model. By default, covariates are fixed to the mean.
#' @param gene_quant Vector of the genes quantiles used to make the plot. We use quantiles instead of fixed values because genetic scores can vary widely depending on the weights, thus looking at quantiles make this simpler. (Default = c(.025,.50,.975))
#' @param env_quant Vector of the environments quantiles used to make the plot. We use quantiles instead of fixed values because environmental scores can vary widely depending on the weights, thus looking at quantiles make this simpler. (Default = c(.025,.50,.975))
#' @param outcome_quant Vector of the outcome quantiles used to make the plot. We use quantiles instead of fixed values because environmental scores can vary widely depending on the weights, thus looking at quantiles make this simpler. (Default = c(.025,.50,.975))
#' @param cols Colors for the slopes with different genetic score. Must be a vector same length as "gene_range". (Default = c("#3288BD", "#CAB176", #D53E4F"))
#' @param ylab Y-axis label (Default = "Outcome")
#' @param xlab X-axis label (Default = "Environment")
#' @param legtitle Title of the Legend for the genes slopes label (Default = "Genetic score")
#' @param leglab Optional vector of labels of the Legend for the genes slopes label
#' @param cex.axis relative scale of axis (Default = 1.9)
#' @param cex.lab relative scale of labels (Default = 2)
#' @param cex.main relative scale overall (Default = 2.2)
#' @param cex.leg relative scale of legend (Default = 2.2)
#' @param xlim X-axis vector of size two with min and max (Default = NULL which leads to min="2.5 percentile" and max="97.5 percentile").
#' @param ylim Y-axis vector of size two with min and max (Default = NULL which leads to min="2.5 percentile" and max="97.5 percentile").
#' @param x_at specific ticks for the X-axis, first and last will be min and max respectively (Default = NULL which leads to 2.5, 50 and 97.5 percentiles).
#' @param y_at specific ticks for the Y-axis, first and last will be min and max respectively (Default = NULL which leads to 2.5, 50 and 97.5 percentiles).
#' @param legend The location may of the legend be specified by setting legend to a single keyword from the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center" (Default = "topleft").
#' @return Returns a list containing the different models (diathesis-stress, differential susceptibility and vantage sensitivity WEAK or STRONG) in order from best to worst for each selected criterion.
#' @param ... Further arguments passed to or from other methods.
#' @examples
#'	train = example_2way(500, 1, seed=777)
#'	fit = LEGIT(train$data, train$G, train$E, y ~ G*E, train$coef_G, train$coef_E)
#'	plot(fit)
#' @import formula.tools graphics
#' @references Alexia Jolicoeur-Martineau, Ashley Wazana, Eszter Szekely, Meir Steiner, Alison S. Fleming, James L. Kennedy, Michael J. Meaney, Celia M.T. Greenwood and the MAVAN team. \emph{Alternating optimization for GxE modelling with weighted genetic and environmental scores: examples from the MAVAN study} (2017). arXiv:1703.08111.
#' @export
"plot.LEGIT"

#' @title Independent Multiple Latent Environmental & Genetic InTeraction (IMLEGIT) model
#' @description Constructs a generalized linear model (glm) with latent variables using alternating optimization. This is an extension of the LEGIT model to accommodate more than 2 latent variables.
#' @param data data.frame of the dataset to be used. 
#' @param latent_var list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise, the latent variables will be named L1, L2, ... (See examples below for more details)
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in the formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param start_latent_var Optional list of starting points for each latent variable (The list must have the same length as the number of latent variables and each element of the list must have the same length as the number of variables of the corresponding latent variable).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param print If FALSE, nothing except warnings will be printed. (Default = TRUE).
#' @return Returns an object of the class "IMLEGIT" which is list containing, in the following order: a glm fit of the main model, a list of the glm fits of the latent variables and a list of the true model parameters (AIC, BIC, rank, df.residual, null.deviance) for which the individual model parts (main, genetic, environmental) don't estimate properly.
#' @examples
#'	train = example_2way(500, 1, seed=777)
#'	fit_best = IMLEGIT(train$data, list(G=train$G, E=train$E), y ~ G*E, 
#'	list(train$coef_G, train$coef_E))
#'	fit_default = IMLEGIT(train$data, list(G=train$G, E=train$E), y ~ G*E)
#'	summary(fit_default)
#'	summary(fit_best)
#'	train = example_3way_3latent(500, 1, seed=777)
#'	fit_best = IMLEGIT(train$data, train$latent_var, y ~ G*E*Z, 
#'	list(train$coef_G, train$coef_E, train$coef_Z))
#'	fit_default = IMLEGIT(train$data, train$latent_var, y ~ G*E*Z)
#'	summary(fit_default)
#'	summary(fit_best)
#' @import formula.tools stats
#' @references Alexia Jolicoeur-Martineau, Ashley Wazana, Eszter Szekely, Meir Steiner, Alison S. Fleming, James L. Kennedy, Michael J. Meaney, Celia M.T. Greenwood and the MAVAN team. \emph{Alternating optimization for GxE modelling with weighted genetic and environmental scores: examples from the MAVAN study} (2017). arXiv:1703.08111.
#' @export
"IMLEGIT"

#' @title Elastic net for variable selection in IMLEGIT model
#' @description [Fast and accurate, highly recommended] Apply Elastic Net (from the glmnet package) with IMLEGIT to obtain the order of variable removal that makes the most sense. The output shows the information criterion at every step, so you can decide which variable to retain. It is significantly faster (seconds/minutes instead of hours) than all other variable selection approaches (except for stepwise) and it is very accurate. Note that, as opposed to LEGIT/IMLEGIT, the parameters of variables inside the latent variables are not L1-normalized; instead, its the main model parameters which are L1-normalized. This is needed to make elastic net works. It doesn't matter in the end, because we only care about which variables were removed and we only output the IMLEGIT models without elastic net penalization.
#' @param data data.frame of the dataset to be used. 
#' @param latent_var list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise, the latent variables will be named L1, L2, ... (See examples below for more details)
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in the formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param latent_var_searched Optional If not null, you must specify a vector containing all indexes of the latent variables you want to use elastic net on. Ex: If latent_var=list(G=genes, E=env), specifying latent_var_search=c(1,2) will use both, latent_var_search=1 will only do it for G, and latent_var_search=2 will only do it for E.
#' @param cross_validation (Optional) If TRUE, will return cross-validation criterion (slower, but very good criterion).
#' @param alpha The elasticnet mixing parameter (between 0 and 1). 1 leads to lasso, 0 leads to ridge. See glmnet package manual for more information. We recommend somewhere betwen .50 and 1.
#' @param standardize If TRUE, standardize all variables inside every latent_var component. Note that if FALSE, glmnet will still standardize and unstandardize, but it will do so for each model (i.e., when at the step of estimating the parameters of latent variable G it standardize them, apply glmnet, then unstandarize them). This means that fixed parameters in the alternating steps are not standardized when standardize=FALSE. In practice, we found that standardize=FALSE leads to weird paths that do not always make sense. In the end, we only care about the order of the variable removal from the glmnet. We highly recommend standardize=TRUE for best results.
#' @param lambda_path Optional vector of all lambda (penalty term for elastic net, see glmnet package manual). By default, we automatically determine it.
#' @param lambda_mult scalar which multiplies the maximum lambda (penalty term for elastic net, see glmnet package manual) from the lambda path determined automatically. Sometimes, the maximum lambda found automatically is too big or too small and you may not want to spend the effort to manually set your own lambda path. This is where this comes in, you can simply scale lambda max up or down. (Default = 1)
#' @param lambda_min minimum lambda (penalty term for elastic net, see glmnet package manual) from the lambda path. (Default = .0001)
#' @param n_lambda Number of lambda (penalty term for elastic net, see glmnet package manual) in lambda path. Make lower for faster training, or higher for more precision. If you have many variables, make it bigger than 100 (Default = 100).
#' @param start_latent_var Optional list of starting points for each latent variable (The list must have the same length as the number of latent variables and each element of the list must have the same length as the number of variables of the corresponding latent variable).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param folds Optional list of vectors containing the fold number for each observation. Bypass cv_iter and cv_folds. Setting your own folds could be important for certain data types like time series or longitudinal data.
#' @param classification Set to TRUE if you are doing classification (binary outcome).
#' @param Huber_p Parameter controlling the Huber cross-validation error (Default = 1.345).
#' @param print If FALSE, nothing except warnings will be printed. (Default = TRUE).
#' @return Returns an object of the class "elastic_net_var_select" which is list containing, in the following order: the criterion at each lambda, the coefficients of the latent variables at each lambda, the fits of each IMLEGIT models for each variable retained at each lambda, and the vector of lambda used.
#' @examples
#'	\dontrun{
#'	N = 1000
#'	train = example_3way(N, sigma=1, logit=FALSE, seed=7)
#'	g1_bad = rbinom(N,1,.30)
#'	g2_bad = rbinom(N,1,.30)
#'	g3_bad = rbinom(N,1,.30)
#'	g4_bad = rbinom(N,1,.30)
#'	g5_bad = rbinom(N,1,.30)
#'	train$G = cbind(train$G, g1_bad, g2_bad, g3_bad, g4_bad, g5_bad)
#'	lv = list(G=train$G, E=train$E)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E)
#'	summary(fit)
#'	best_model(fit, criterion="BIC")
#'  # Instead of taking the best, if you want the model with "Model index"=17 from summary, do
#   fit_mychoice = fit$fit[[17]]
#'	plot(fit)
#'	# With Cross-validation
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, cross_validation=TRUE, cv_iter=1, cv_folds=5)
#'	best_model(fit, criterion="cv_R2")
#'	# Elastic net only applied on G
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(1))
#'	# Elastic net only applied on E
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(2))
#'	# Most E variables not removed, use lambda_mult > 1 to remove more
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(2), lambda_mult=5)
#'	# Lasso (only L1 regularization)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, alpha=1)
#'	# Want more lambdas (useful if # of variables is large)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, n_lambda = 200)
#'	}
#' @import formula.tools stats
#' @references Alexia Jolicoeur-Martineau, Ashley Wazana, Eszter Szekely, Meir Steiner, Alison S. Fleming, James L. Kennedy, Michael J. Meaney, Celia M.T. Greenwood and the MAVAN team. \emph{Alternating optimization for GxE modelling with weighted genetic and environmental scores: examples from the MAVAN study} (2017). arXiv:1703.08111.
#' @export
"elastic_net_var_select"


#' @title Independent Multiple Latent Environmental & Genetic InTeraction (IMLEGIT) model with Elastic Net on the latent variables. Do not use on it's own, use elastic_net_var_select instead.
#' @description Constructs a generalized linear model (glm) with latent variables using alternating optimization. This is an extension of the LEGIT model to accommodate more than 2 latent variables. Note that, as opposed to LEGIT/IMLEGIT, the parameters of variables inside the latent variables are not L1-normalized; instead, its the main model parameters which are L1-normalized. This is needed to make elastic net works. It doesn't matter in the end, because we only care about which variables were removed and we only give the IMLEGIT models without elastic net penalization.
#' @param data data.frame of the dataset to be used. 
#' @param latent_var list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise, the latent variables will be named L1, L2, ... (See examples below for more details)
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in the formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param latent_var_searched Optional If not null, you must specify a vector containing all indexes of the latent variables you want to use elastic net on. Ex: If latent_var=list(G=genes, E=env), specifying latent_var_search=c(1,2) will use both, latent_var_search=1 will only do it for G, and latent_var_search=2 will only do it for E.
#' @param cross_validation If TRUE, will return cross-validation criterion (slower)
#' @param alpha The elasticnet mixing parameter (between 0 and 1). 1 leads to lasso, 0 leads to ridge. See glmnet package manual for more information. We recommend somewhere betwen .50 and 1.
#' @param lambda Lambda (penalty term for elastic net, see glmnet package manual) (Default = .0001)
#' @param start_latent_var Optional list of starting points for each latent variable (The list must have the same length as the number of latent variables and each element of the list must have the same length as the number of variables of the corresponding latent variable).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param folds Optional list of vectors containing the fold number for each observation. Bypass cv_iter and cv_folds. Setting your own folds could be important for certain data types like time series or longitudinal data.
#' @param classification Set to TRUE if you are doing classification (binary outcome).
#' @param Huber_p Parameter controlling the Huber cross-validation error (Default = 1.345).
#' @param print If FALSE, nothing except warnings will be printed. (Default = TRUE).
#' @param warn If FALSE, it will not show warnings when all variables inside a latent variable are removed. This serves to prevent lots of warning when running elastic_net_var_select (Default = TRUE).
#' @param family_string Optional String version of the family (gaussian leads to "gaussian"). This is only needed when using elastic_net_var_select. Please ignore this.
#' @return Returns a list containing, in the following order: a IMLEGIT model, the coefficients of the variables in the latent variables from glmnet models, and the cross-validation results (if asked).
#' @import formula.tools stats
#' @references Alexia Jolicoeur-Martineau, Ashley Wazana, Eszter Szekely, Meir Steiner, Alison S. Fleming, James L. Kennedy, Michael J. Meaney, Celia M.T. Greenwood and the MAVAN team. \emph{Alternating optimization for GxE modelling with weighted genetic and environmental scores: examples from the MAVAN study} (2017). arXiv:1703.08111.
#' @export
"IMLEGIT_net"

#' @title Plot function for the output of elastic_net_var_select
#' @description Plot of the coefficients of variables inside the latent variables with respect to the log(lambda). This is your typical elastic-net plot.
#' @param x An object of class "elastic_net_var_select", usually, a result of a call to elastic_net_var_select.
#' @param lwd Thickness of the lines (Default = 2)
#' @param start At which lambda to start (from large lambda to small lambda). If start is not 1, we remove some of the large lambda, this can make plot easier to visualize (Default = 1).
#' @param ... Further arguments passed to or from other methods.
#' @return Returns the plot of the coefficients of variables inside the latent variables with respect to the log(lambda).
#' @examples
#'	\dontrun{
#'	N = 1000
#'	train = example_3way(N, sigma=1, logit=FALSE, seed=7)
#'	g1_bad = rbinom(N,1,.30)
#'	g2_bad = rbinom(N,1,.30)
#'	g3_bad = rbinom(N,1,.30)
#'	g4_bad = rbinom(N,1,.30)
#'	g5_bad = rbinom(N,1,.30)
#'	train$G = cbind(train$G, g1_bad, g2_bad, g3_bad, g4_bad, g5_bad)
#'	lv = list(G=train$G, E=train$E)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E)
#'	summary(fit)
#'	best_model(fit, criterion="BIC")
#'  # Instead of taking the best, if you want the model with "Model index"=17 from summary, do
#   fit_mychoice = fit$fit[[17]]
#'	plot(fit)
#'	# With Cross-validation
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, cross_validation=TRUE, cv_iter=1, cv_folds=5)
#'	best_model(fit, criterion="cv_R2")
#'	# Elastic net only applied on G
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(1))
#'	# Elastic net only applied on E
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(2))
#'	# Most E variables not removed, use lambda_mult > 1 to remove more
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(2), lambda_mult=5)
#'	# Lasso (only L1 regularization)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, alpha=1)
#'	# Want more lambdas (useful if # of variables is large)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, n_lambda = 200)
#'	}
#' @import formula.tools stats
#' @references Alexia Jolicoeur-Martineau, Ashley Wazana, Eszter Szekely, Meir Steiner, Alison S. Fleming, James L. Kennedy, Michael J. Meaney, Celia M.T. Greenwood and the MAVAN team. \emph{Alternating optimization for GxE modelling with weighted genetic and environmental scores: examples from the MAVAN study} (2017). arXiv:1703.08111.
#' @export
"plot.elastic_net_var_select"

#' @title Summary function for the output of elastic_net_var_select
#' @description Summary function for the output of elastic_net_var_select
#' @param object An object of class "elastic_net_var_select", usually, a result of a call to elastic_net_var_select.
#' @param ... Further arguments passed to or from other methods.
#' @return Returns the unique IMLEGIT models resulting from the glmnet path with associated information. Also gives the cross-validation information if asked.
#' @examples
#'	\dontrun{
#'	N = 1000
#'	train = example_3way(N, sigma=1, logit=FALSE, seed=7)
#'	g1_bad = rbinom(N,1,.30)
#'	g2_bad = rbinom(N,1,.30)
#'	g3_bad = rbinom(N,1,.30)
#'	g4_bad = rbinom(N,1,.30)
#'	g5_bad = rbinom(N,1,.30)
#'	train$G = cbind(train$G, g1_bad, g2_bad, g3_bad, g4_bad, g5_bad)
#'	lv = list(G=train$G, E=train$E)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E)
#'	summary(fit)
#'	best_model(fit, criterion="BIC")
#'  # Instead of taking the best, if you want the model with "Model index"=17 from summary, do
#   fit_mychoice = fit$fit[[17]]
#'	plot(fit)
#'	# With Cross-validation
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, cross_validation=TRUE, cv_iter=1, cv_folds=5)
#'	best_model(fit, criterion="cv_R2")
#'	# Elastic net only applied on G
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(1))
#'	# Elastic net only applied on E
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(2))
#'	# Most E variables not removed, use lambda_mult > 1 to remove more
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(2), lambda_mult=5)
#'	# Lasso (only L1 regularization)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, alpha=1)
#'	# Want more lambdas (useful if # of variables is large)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, n_lambda = 200)
#'	}
#' @import formula.tools stats
#' @references Alexia Jolicoeur-Martineau, Ashley Wazana, Eszter Szekely, Meir Steiner, Alison S. Fleming, James L. Kennedy, Michael J. Meaney, Celia M.T. Greenwood and the MAVAN team. \emph{Alternating optimization for GxE modelling with weighted genetic and environmental scores: examples from the MAVAN study} (2017). arXiv:1703.08111.
#' @export
"summary.elastic_net_var_select"

#' @title Best model
#' @description Best model
#' @param object An object
#' @param ... Further arguments passed to or from other methods.
#' @return Best model
#' @export
"best_model"

#' @title Best model from elastic net variable selection
#' @description Best model from elastic net variable selection (based on selected criteria)
#' @param object An object of class "elastic_net_var_select", usually, a result of a call to elastic_net_var_select.
#' @param criterion Criteria used to determine which model is the best. If \code{search_criterion="AIC"}, uses the AIC, if \code{search_criterion="AICc"}, uses the AICc, if \code{search_criterion="BIC"}, uses the BIC, if \code{search_criterion="cv_R2"}, uses the cross-validation R-squared, if \cr \code{search_criterion="cv_AUC"}, uses the cross-validated AUC, if \code{search_criterion="cv_Huber"}, uses the Huber cross-validation error, if \code{search_criterion="cv_L1"}, uses the L1-norm cross-validation error (Default = "AIC"). The Huber and L1-norm cross-validation errors are alternatives to the usual cross-validation L2-norm error (which the \eqn{R^2} is based on) that are more resistant to outliers. For all criterion, lower is better, with the exception of \code{search_criterion="cv_R2"} and \code{search_criterion="cv_AUC"}.
#' @param ... Further arguments passed to or from other methods.
#' @return Returns the best IMLEGIT model resulting from the glmnet path with associated information.
#' @examples
#'	\dontrun{
#'	N = 1000
#'	train = example_3way(N, sigma=1, logit=FALSE, seed=7)
#'	g1_bad = rbinom(N,1,.30)
#'	g2_bad = rbinom(N,1,.30)
#'	g3_bad = rbinom(N,1,.30)
#'	g4_bad = rbinom(N,1,.30)
#'	g5_bad = rbinom(N,1,.30)
#'	train$G = cbind(train$G, g1_bad, g2_bad, g3_bad, g4_bad, g5_bad)
#'	lv = list(G=train$G, E=train$E)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E)
#'	summary(fit)
#'	best_model(fit, criterion="BIC")
#'  # Instead of taking the best, if you want the model with "Model index"=17 from summary, do
#   fit_mychoice = fit$fit[[17]]
#'	plot(fit)
#'	# With Cross-validation
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, cross_validation=TRUE, cv_iter=1, cv_folds=5)
#'	best_model(fit, criterion="cv_R2")
#'	# Elastic net only applied on G
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(1))
#'	# Elastic net only applied on E
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(2))
#'	# Most E variables not removed, use lambda_mult > 1 to remove more
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, c(2), lambda_mult=5)
#'	# Lasso (only L1 regularization)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, alpha=1)
#'	# Want more lambdas (useful if # of variables is large)
#'	fit = elastic_net_var_select(train$data, lv, y ~ G*E, n_lambda = 200)
#'	}
#' @import formula.tools stats
#' @references Alexia Jolicoeur-Martineau, Ashley Wazana, Eszter Szekely, Meir Steiner, Alison S. Fleming, James L. Kennedy, Michael J. Meaney, Celia M.T. Greenwood and the MAVAN team. \emph{Alternating optimization for GxE modelling with weighted genetic and environmental scores: examples from the MAVAN study} (2017). arXiv:1703.08111.
#' @export
"best_model.elastic_net_var_select"

#' @title Predictions of LEGIT fits
#' @description Predictions of LEGIT fits.
#' @param object An object of class "LEGIT", usually, a result of a call to LEGIT.
#' @param data data.frame of the dataset to be used.
#' @param genes data.frame of the variables inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic).
#' @param env data.frame of the variables inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental).
#' @param ... Further arguments passed to or from other methods.
#' @return Returns a vector with the predicted values.
#' @examples
#'	train = example_2way(250, 1, seed=777)
#'	test = example_2way(100, 1, seed=666)
#'	fit = LEGIT(train$data, train$G, train$E, y ~ G*E)
#'	ssres = sum((test$data$y - predict(fit, test$data, test$G, test$E))^2)
#'	sstotal = sum((test$data$y - mean(test$data$y))^2)
#'	R2 = 1 - ssres/sstotal
#' @export
"predict.LEGIT"

#' @title Predictions of IMLEGIT fits
#' @description Predictions of IMLEGIT fits.
#' @param object An object of class "IMLEGIT", usually, a result of a call to IMLEGIT.
#' @param data data.frame of the dataset to be used.
#' @param latent_var list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise, the latent variables will be named L1, L2, ...
#' @param ... Further arguments passed to or from other methods.
#' @return Returns a vector with the predicted values.
#' @examples
#'	train = example_2way(250, 1, seed=777)
#'	test = example_2way(100, 1, seed=666)
#'	fit = IMLEGIT(train$data, list(G=train$G, E=train$E), y ~ G*E)
#'	ssres = sum((test$data$y - predict(fit, test$data, list(G=test$G, E=test$E)))^2)
#'	sstotal = sum((test$data$y - mean(test$data$y))^2)
#'	R2 = 1 - ssres/sstotal
#'	R2
#' @export
"predict.IMLEGIT"

#' @title Summarizing LEGIT fits
#' @description Shows the summary for all parts (main, genetic, environmental) of the LEGIT model.
#' @param object An object of class "LEGIT", usually, a result of a call to LEGIT.
#' @param ... Further arguments passed to or from other methods.
#' @return Returns a list of objects of class "summary.glm" containing the summary of each parts (main, genetic, environmental) of the model.
#' @examples
#' 	train = example_2way(250, 1, seed=777)
#'	fit_default = LEGIT(train$data, train$G, train$E, y ~ G*E)
#'	summary(fit_default)
#' @export
"summary.LEGIT"

#' @title Summarizing IMLEGIT fits
#' @description Shows the summary for all parts (main and latent variables) of the LEGIT model.
#' @param object An object of class "IMLEGIT", usually, a result of a call to IMLEGIT.
#' @param ... Further arguments passed to or from other methods.
#' @return Returns a list of objects of class "summary.glm" containing the summary of each parts (main and latent variables) of the model.
#' @examples
#' 	train = example_2way(250, 1, seed=777)
#'	fit_default = IMLEGIT(train$data, list(G=train$G, E=train$E), y ~ G*E)
#'	summary(fit_default)
#' @export
"summary.IMLEGIT"

#' @title Cross-validation for the LEGIT model
#' @description Uses cross-validation on the LEGIT model. Note that this is not a very fast implementation since it was written in R.
#' @param data data.frame of the dataset to be used.
#' @param genes data.frame of the variables inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic).
#' @param env data.frame of the variables inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental).
#' @param formula Model formula. Use \emph{E} for the environmental score and \emph{G} for the genetic score. Do not manually code interactions, write them in the formula instead (ex: G*E*z or G:E:z).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param folds Optional list of vectors containing the fold number for each observation. Bypass cv_iter and cv_folds. Setting your own folds could be important for certain data types like time series or longitudinal data.
#' @param classification Set to TRUE if you are doing classification (binary outcome).
#' @param start_genes Optional starting points for genetic score (must be the same length as the number of columns of \code{genes}).
#' @param start_env Optional starting points for environmental score (must be the same length as the number of columns of \code{env}).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param seed Seed for cross-validation folds.
#' @param Huber_p Parameter controlling the Huber cross-validation error (Default = 1.345).
#' @param id Optional id of observations, can be a vector or data.frame (only used when returning list of possible outliers).
#' @param crossover If not NULL, estimates the crossover point of \emph{E} using the provided value as starting point (To test for diathesis-stress vs differential susceptibility).
#' @param crossover_fixed If TRUE, instead of estimating the crossover point of E, we force/fix it to the value of "crossover". (Used when creating a diathes-stress model) (Default = FALSE).
#' @return If \code{classification} = FALSE, returns a list containing, in the following order: a vector of the cross-validated \eqn{R^2} at each iteration, a vector of the Huber cross-validation error at each iteration, a vector of the L1-norm cross-validation error at each iteration, a matrix of the possible outliers (standardized residuals > 2.5 or < -2.5) and their corresponding standardized residuals and standardized pearson residuals. If \code{classification} = TRUE, returns a list containing, in the following order: a vector of the cross-validated \eqn{R^2} at each iteration, a vector of the Huber cross-validation error at each iteration, a vector of the L1-norm cross-validation error at each iteration, a vector of the AUC at each iteration, a matrix of the best choice of threshold (based on Youden index) and the corresponding specificity and sensitivity at each iteration, and a list of objects of class "roc" (to be able to make roc curve plots) at each iteration. The Huber and L1-norm cross-validation errors are alternatives to the usual cross-validation L2-norm error (which the \eqn{R^2} is based on) that are more resistant to outliers, the lower the values the better.
#' @examples
#'	\dontrun{
#'	train = example_3way(250, 2.5, seed=777)
#'	# Cross-validation 4 times with 5 Folds
#'	cv_5folds = LEGIT_cv(train$data, train$G, train$E, y ~ G*E*z, cv_iter=4, cv_folds=5)
#'	cv_5folds
#'	# Leave-one-out cross-validation (Note: very slow)
#'	cv_loo = LEGIT_cv(train$data, train$G, train$E, y ~ G*E*z, cv_iter=1, cv_folds=250)
#'	cv_loo
#'	# Cross-validation 4 times with 5 Folds (binary outcome)
#'	train_bin = example_2way(500, 2.5, logit=TRUE, seed=777)
#'	cv_5folds_bin = LEGIT_cv(train_bin$data, train_bin$G, train_bin$E, y ~ G*E, 
#'	cv_iter=4, cv_folds=5, classification=TRUE, family=binomial)
#'	cv_5folds_bin
#'	par(mfrow=c(2,2))
#'	pROC::plot.roc(cv_5folds_bin$roc_curve[[1]])
#'	pROC::plot.roc(cv_5folds_bin$roc_curve[[2]])
#'	pROC::plot.roc(cv_5folds_bin$roc_curve[[3]])
#'	pROC::plot.roc(cv_5folds_bin$roc_curve[[4]])
#'	}
#' @references Denis Heng-Yan Leung. \emph{Cross-validation in nonparametric regression with outliers.} Annals of Statistics (2005): 2291-2310.
#' @export
"LEGIT_cv"

#' @title Cross-validation for the IMLEGIT model
#' @description Uses cross-validation on the IMLEGIT model. Note that this is not a very fast implementation since it was written in R.
#' @param data data.frame of the dataset to be used.
#' @param latent_var list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise, the latent variables will be named L1, L2, ...
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in the formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param folds Optional list of vectors containing the fold number for each observation. Bypass cv_iter and cv_folds. Setting your own folds could be important for certain data types like time series or longitudinal data.
#' @param classification Set to TRUE if you are doing classification (binary outcome).
#' @param start_latent_var Optional list of starting points for each latent variable (The list must have the same length as the number of latent variables and each element of the list must have the same length as the number of variables of the corresponding latent variable).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param seed Seed for cross-validation folds.
#' @param Huber_p Parameter controlling the Huber cross-validation error (Default = 1.345).
#' @param id Optional id of observations, can be a vector or data.frame (only used when returning list of possible outliers).
#' @return If \code{classification} = FALSE, returns a list containing, in the following order: a vector of the cross-validated \eqn{R^2} at each iteration, a vector of the Huber cross-validation error at each iteration, a vector of the L1-norm cross-validation error at each iteration, a matrix of the possible outliers (standardized residuals > 2.5 or < -2.5) and their corresponding standardized residuals and standardized pearson residuals. If \code{classification} = TRUE, returns a list containing, in the following order: a vector of the cross-validated \eqn{R^2} at each iteration, a vector of the Huber cross-validation error at each iteration, a vector of the L1-norm cross-validation error at each iteration, a vector of the AUC at each iteration, a matrix of the best choice of threshold (based on Youden index) and the corresponding specificity and sensitivity at each iteration, and a list of objects of class "roc" (to be able to make roc curve plots) at each iteration. The Huber and L1-norm cross-validation errors are alternatives to the usual cross-validation L2-norm error (which the \eqn{R^2} is based on) that are more resistant to outliers, the lower the values the better.
#' @examples
#'	\dontrun{
#'	train = example_3way_3latent(250, 1, seed=777)
#'	# Cross-validation 4 times with 5 Folds
#'	cv_5folds = IMLEGIT_cv(train$data, train$latent_var, y ~ G*E*Z, cv_iter=4, cv_folds=5)
#'	cv_5folds
#'	# Leave-one-out cross-validation (Note: very slow)
#'	cv_loo = IMLEGIT_cv(train$data, train$latent_var, y ~ G*E*Z, cv_iter=1, cv_folds=250)
#'	cv_loo
#'	# Cross-validation 4 times with 5 Folds (binary outcome)
#'	train_bin = example_2way(500, 2.5, logit=TRUE, seed=777)
#'	cv_5folds_bin = IMLEGIT_cv(train_bin$data, list(G=train_bin$G, E=train_bin$E), y ~ G*E, 
#'	cv_iter=4, cv_folds=5, classification=TRUE, family=binomial)
#'	cv_5folds_bin
#'	par(mfrow=c(2,2))
#'	pROC::plot.roc(cv_5folds_bin$roc_curve[[1]])
#'	pROC::plot.roc(cv_5folds_bin$roc_curve[[2]])
#'	pROC::plot.roc(cv_5folds_bin$roc_curve[[3]])
#'	pROC::plot.roc(cv_5folds_bin$roc_curve[[4]])
#'	}
#' @references Denis Heng-Yan Leung. \emph{Cross-validation in nonparametric regression with outliers.} Annals of Statistics (2005): 2291-2310.
#' @export
"IMLEGIT_cv"

#' Internal function that does the forward step for the stepwise function.
#' @param empty_start_dataset If TRUE, the initial dataset is empty.
#' @param fit Current best fit.
#' @param ... Same parameters as in the stepwise function.
#' @return Returns fit, start_genes, start_env and genes_current, genes_toadd if search="genes" or env_current and env_toadd if search="env".
#' @keywords internal
"forward_step"

#' Internal function that does the forward step for the stepwise function.
#' @param empty_start_dataset If TRUE, the initial dataset is empty.
#' @param fit Current best fit.
#' @param ... Same parameters as in the stepwise function.
#' @return Returns fit, start_latent_var, latent_var_current and latent_var_toadd.
#' @keywords internal
"forward_step_IM"

#' Internal function that does the backward step for the stepwise IM function.
#' @param empty_start_dataset If TRUE, the initial dataset is empty.
#' @param fit Current best fit.
#' @param ... Same parameters as in the stepwise function.
#' @return Returns fit, start_genes, start_env and genes_current, genes_dropped if search="genes" or env_current and env_dropped if search="env".
#' @keywords internal
"backward_step"

#' Internal function that does the backward step for the stepwise IM function.
#' @param empty_start_dataset If TRUE, the initial dataset is empty.
#' @param fit Current best fit.
#' @param ... Same parameters as in the stepwise function.
#' @return Returns fit, start_latent_var, latent_var_current and latent_var_dropped.
#' @keywords internal
"backward_step_IM"

#' @title Stepwise search for the best subset of genetic variants or environments with the LEGIT model
#' @description [Fast, recommended for small number of variables] Adds the best variable or drops the worst variable one at a time in the genetic (if \code{search="genes"}) or environmental score (if \code{search="env"}). You can select the desired search criterion (AIC, BIC, cross-validation error, cross-validation AUC) to determine which variable is the best/worst and should be added/dropped. Note that when the number of variables in \emph{G} and \emph{E} is large, this does not generally converge to the optimal subset, this function is only recommended when you have a small number of variables (e.g. 2 environments, 6 genetic variants). If using cross-validation (\code{search_criterion="cv"} or \code{search_criterion="cv_AUC"}), to prevent cross-validating with each variable (extremely slow), we recommend setting a p-value threshold (\code{p_threshold}) and forcing the algorithm not to look at models with bigger AIC (\code{exclude_worse_AIC=TRUE}).
#' @param data data.frame of the dataset to be used.
#' @param formula Model formula. Use \emph{E} for the environmental score and \emph{G} for the genetic score. Do not manually code interactions, write them in the formula instead (ex: G*E*z or G:E:z).
#' @param interactive_mode If TRUE, uses interactive mode. In interactive mode, at each iteration, the user is shown the AIC, BIC, p-value and also the cross-validation \eqn{R^2} if \code{search_criterion="cv"} and the cross-validation AUC if \code{search_criterion="cv_AUC"} for the best 5 variables. The user must then enter a number between 1 and 5 to select the variable to be added, entering anything else will stop the search.
#' @param genes_original data.frame of the variables inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic).
#' @param env_original data.frame of the variables inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental).
#' @param genes_extra data.frame of the additionnal variables to try including inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic). Set to NULL if using a backward search.
#' @param env_extra data.frame of the variables to try including inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental). Set to NULL if using a backward search.
#' @param search_type If \code{search_type="forward"}, uses a forward search. If \code{search_type="backward"}, uses backward search. If \code{search_type="bidirectional-forward"}, uses bidirectional search (that starts as a forward search). If \code{search_type="bidirectional-backward"}, uses bidirectional search (that starts as a backward search).
#' @param search If \code{search="genes"}, uses a stepwise search for the genetic score variables. If \code{search="env"}, uses a stepwise search for the environmental score variables. If \code{search="both"}, uses a stepwise search for both the gene and environmental score variables (Default = "both").
#' @param search_criterion Criterion used to determine which variable is the best to add or worst to drop. If \code{search_criterion="AIC"}, uses the AIC, if \code{search_criterion="AICc"}, uses the AICc, if \code{search_criterion="BIC"}, uses the BIC, if \code{search_criterion="cv"}, uses the cross-validation error, if \cr \code{search_criterion="cv_AUC"}, uses the cross-validated AUC, if \code{search_criterion="cv_Huber"}, uses the Huber cross-validation error, if \code{search_criterion="cv_L1"}, uses the L1-norm cross-validation error (Default = "AIC"). The Huber and L1-norm cross-validation errors are alternatives to the usual cross-validation L2-norm error (which the \eqn{R^2} is based on) that are more resistant to outliers, the lower the values the better.
#' @param forward_exclude_p_bigger If p-value > \code{forward_exclude_p_bigger}, we do not consider the variable for inclusion in the forward steps (Default = .20). This is an exclusion option which purpose is skipping variables that are likely not worth looking to make the algorithm faster, especially with cross-validation. Set to 1 to prevent any exclusion here.
#' @param backward_exclude_p_smaller If p-value < \code{backward_exclude_p_smaller}, we do not consider the variable for removal in the backward steps (Default = .01). This is an exclusion option which purpose is skipping variables that are likely not worth looking to make the algorithm faster, especially with cross-validation. Set to 0 to prevent any exclusion here.
#' @param exclude_worse_AIC If AIC with variable > AIC without variable, we ignore the variable (Default = TRUE). This is an exclusion option which purpose is skipping variables that are likely not worth looking to make the algorithm faster, especially with cross-validation. Set to FALSE to prevent any exclusion here.
#' @param max_steps Maximum number of steps taken (Default = 50).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param folds Optional list of vectors containing the fold number for each observation. Bypass cv_iter and cv_folds. Setting your own folds could be important for certain data types like time series or longitudinal data.
#' @param classification Set to TRUE if you are doing classification (binary outcome).
#' @param start_genes Optional starting points for genetic score (must be the same length as the number of columns of \code{genes}).
#' @param start_env Optional starting points for environmental score (must be the same length as the number of columns of \code{env}).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param seed Seed for cross-validation folds.
#' @param print If TRUE, print all the steps and notes/warnings. Highly recommended unless you are batch running multiple stepwise searchs. (Default=TRUE).
#' @param Huber_p Parameter controlling the Huber cross-validation error (Default = 1.345).
#' @param remove_miss If TRUE, remove missing data completely, otherwise missing data is only removed when adding or dropping a variable (Default = FALSE).
#' @return Returns an object of the class "LEGIT" which is list containing, in the following order: a glm fit of the main model, a glm fit of the genetic score, a glm fit of the environmental score, a list of the true model parameters (AIC, BIC, rank, df.residual, null.deviance) for which the individual model parts (main, genetic, environmental) don't estimate properly.
#' @examples
#'	\dontrun{
#'	## Continuous example
#'	train = example_3way(250, 2.5, seed=777)
#'	# Forward search for genes based on BIC (in interactive mode)
#'	forward_genes_BIC = stepwise_search(train$data, genes_extra=train$G, env_original=train$E,
#'	formula=y ~ E*G*z,search_type="forward", search="genes", search_criterion="BIC",
#'	interactive_mode=TRUE)
#'	# Bidirectional-backward search for environments based on cross-validation error
#'	bidir_backward_env_cv = stepwise_search(train$data, genes_original=train$G, env_original=train$E,
#'	formula=y ~ E*G*z,search_type="bidirectional-backward", search="env", search_criterion="cv")
#'	## Binary example
#'	train_bin = example_2way(500, 2.5, logit=TRUE, seed=777)
#'	# Forward search for genes based on cross-validated AUC (in interactive mode)
#'	forward_genes_AUC = stepwise_search(train_bin$data, genes_extra=train_bin$G, 
#'	env_original=train_bin$E, formula=y ~ E*G,search_type="forward", search="genes", 
#'	search_criterion="cv_AUC", classification=TRUE, family=binomial, interactive_mode=TRUE)
#'	# Forward search for genes based on AIC
#'	bidir_forward_genes_AIC = stepwise_search(train_bin$data, genes_extra=train_bin$G, 
#'	env_original=train_bin$E, formula=y ~ E*G,search_type="bidirectional-forward", search="genes", 
#'	search_criterion="AIC", classification=TRUE, family=binomial)
#'	}
#' @export
"stepwise_search"

#' @title Stepwise search for the best subset of elements in the latent variables with the IMLEGIT model
#' @description [Fast, recommended  when the number of variables is small] Adds the best variable or drops the worst variable one at a time in the latent variables. You can select the desired search criterion (AIC, BIC, cross-validation error, cross-validation AUC) to determine which variable is the best/worst and should be added/dropped. Note that when the number of variables in \emph{G} and \emph{E} is large, this does not generally converge to the optimal subset, this function is only recommended when you have a small number of variables (e.g. 2 environments, 6 genetic variants). If using cross-validation (\code{search_criterion="cv"} or \code{search_criterion="cv_AUC"}), to prevent cross-validating with each variable (extremely slow), we recommend setting a p-value threshold (\code{p_threshold}) and forcing the algorithm not to look at models with bigger AIC (\code{exclude_worse_AIC=TRUE}).
#' @param data data.frame of the dataset to be used.
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in the formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param interactive_mode If TRUE, uses interactive mode. In interactive mode, at each iteration, the user is shown the AIC, BIC, p-value and also the cross-validation \eqn{R^2} if \code{search_criterion="cv"} and the cross-validation AUC if \code{search_criterion="cv_AUC"} for the best 5 variables. The user must then enter a number between 1 and 5 to select the variable to be added, entering anything else will stop the search.
#' @param latent_var_original list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise, the latent variables will be named L1, L2, ...
#' @param latent_var_extra list of data.frame (with the same structure as latent_var_original) containing the additionnal elements to try including inside the latent variables. Set to NULL if using a backward search.
#' @param search_type If \code{search_type="forward"}, uses a forward search. If \code{search_type="backward"}, uses backward search. If \code{search_type="bidirectional-forward"}, uses bidirectional search (that starts as a forward search). If \code{search_type="bidirectional-backward"}, uses bidirectional search (that starts as a backward search).
#' @param search If \code{search=0}, uses a stepwise search for all latent variables. Otherwise, if search = i, uses a stepwise search on the i-th latent variable (Default = 0).
#' @param search_criterion Criterion used to determine which variable is the best to add or worst to drop. If \code{search_criterion="AIC"}, uses the AIC, if \code{search_criterion="AICc"}, uses the AICc, if \code{search_criterion="BIC"}, uses the BIC, if \code{search_criterion="cv"}, uses the cross-validation error, if \cr \code{search_criterion="cv_AUC"}, uses the cross-validated AUC, if \code{search_criterion="cv_Huber"}, uses the Huber cross-validation error, if \code{search_criterion="cv_L1"}, uses the L1-norm cross-validation error (Default = "AIC"). The Huber and L1-norm cross-validation errors are alternatives to the usual cross-validation L2-norm error (which the \eqn{R^2} is based on) that are more resistant to outliers, the lower the values the better.
#' @param forward_exclude_p_bigger If p-value > \code{forward_exclude_p_bigger}, we do not consider the variable for inclusion in the forward steps (Default = .20). This is an exclusion option which purpose is skipping variables that are likely not worth looking to make the algorithm faster, especially with cross-validation. Set to 1 to prevent any exclusion here.
#' @param backward_exclude_p_smaller If p-value < \code{backward_exclude_p_smaller}, we do not consider the variable for removal in the backward steps (Default = .01). This is an exclusion option which purpose is skipping variables that are likely not worth looking to make the algorithm faster, especially with cross-validation. Set to 0 to prevent any exclusion here.
#' @param exclude_worse_AIC If AIC with variable > AIC without variable, we ignore the variable (Default = TRUE). This is an exclusion option which purpose is skipping variables that are likely not worth looking to make the algorithm faster, especially with cross-validation. Set to FALSE to prevent any exclusion here.
#' @param max_steps Maximum number of steps taken (Default = 50).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param folds Optional list of vectors containing the fold number for each observation. Bypass cv_iter and cv_folds. Setting your own folds could be important for certain data types like time series or longitudinal data.
#' @param classification Set to TRUE if you are doing classification (binary outcome).
#' @param start_latent_var Optional list of starting points for each latent variable (The list must have the same length as the number of latent variables and each element of the list must have the same length as the number of variables of the corresponding latent variable).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param seed Seed for cross-validation folds.
#' @param print If TRUE, print all the steps and notes/warnings. Highly recommended unless you are batch running multiple stepwise searchs. (Default=TRUE).
#' @param Huber_p Parameter controlling the Huber cross-validation error (Default = 1.345).
#' @param remove_miss If TRUE, remove missing data completely, otherwise missing data is only removed when adding or dropping a variable (Default = FALSE).
#' @return Returns an object of the class "IMLEGIT" which is list containing, in the following order: a glm fit of the main model, a list of the glm fits of the latent variables and a list of the true model parameters (AIC, BIC, rank, df.residual, null.deviance) for which the individual model parts (main, genetic, environmental) don't estimate properly.
#' @examples
#'	\dontrun{
#'	## Example
#'	train = example_3way_3latent(250, 1, seed=777)
#'	# Forward search for genes based on BIC (in interactive mode)
#'	forward_genes_BIC = stepwise_search_IM(train$data, 
#'	latent_var_original=list(G=NULL, E=train$latent_var$E, Z=train$latent_var$Z),
#'	latent_var_extra=list(G=train$latent_var$G,E=NULL,Z=NULL), 
#'	formula=y ~ E*G*Z,search_type="forward", search=1, search_criterion="BIC",
#'	interactive_mode=TRUE)
#'	# Bidirectional-backward search for everything based on AIC
#'	bidir_backward_AIC = stepwise_search_IM(train$data, latent_var_extra=NULL, 
#'	latent_var_original=train$latent_var,
#'	formula=y ~ E*G*Z,search_type="bidirectional-backward", search=0, search_criterion="AIC")
#'	}
#' @export
"stepwise_search_IM"

#' @title Bootstrap variable selection (for IMLEGIT)
#' @description [Very slow, not recommended] Creates bootstrap samples, runs a stepwise search on all of them and then reports the percentage of times that each variable was selected. This is very computationally demanding. With small sample sizes, variable selection can be unstable and bootstrap can be used to give us an idea of the degree of certitude that a variable should be included or not.
#' @param data data.frame of the dataset to be used.
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in the formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param boot_iter number of bootstrap samples (Default = 1000).
#' @param boot_size Optional size of the bootstrapped samples (Default = number of observations).
#' @param boot_group Optional vector which represents the group associated with each observation. Sampling will be done by group instead of by observations (very important if you have longitudinal data). The sample sizes of the bootstrap samples might differ by up to "\code{boot_size} - maximum group size" observations.
#' @param latent_var_original list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise, the latent variables will be named L1, L2, ...
#' @param latent_var_extra list of data.frame (with the same structure as latent_var_original) containing the additional elements to try including inside the latent variables. Set to NULL if using a backward search.
#' @param search_type If \code{search_type="forward"}, uses a forward search. If \code{search_type="backward"}, uses backward search. If \code{search_type="bidirectional-forward"}, uses bidirectional search (that starts as a forward search). If \code{search_type="bidirectional-backward"}, uses bidirectional search (that starts as a backward search).
#' @param search If \code{search=0}, uses a stepwise search for all latent variables. Otherwise, if search = i, uses a stepwise search on the i-th latent variable (Default = 0).
#' @param search_criterion Criterion used to determine which variable is the best to add or worst to drop. If \code{search_criterion="AIC"}, uses the AIC, if \code{search_criterion="AICc"}, uses the AICc, if \code{search_criterion="BIC"}, uses the BIC (Default = "AIC").
#' @param forward_exclude_p_bigger If p-value > \code{forward_exclude_p_bigger}, we do not consider the variable for inclusion in the forward steps (Default = .20). This is an exclusion option which purpose is skipping variables that are likely not worth looking to make the algorithm faster, especially with cross-validation. Set to 1 to prevent any exclusion here.
#' @param backward_exclude_p_smaller If p-value < \code{backward_exclude_p_smaller}, we do not consider the variable for removal in the backward steps (Default = .01). This is an exclusion option which purpose is skipping variables that are likely not worth looking to make the algorithm faster, especially with cross-validation. Set to 0 to prevent any exclusion here.
#' @param exclude_worse_AIC If AIC with variable > AIC without variable, we ignore the variable (Default = TRUE). This is an exclusion option which purpose is skipping variables that are likely not worth looking to make the algorithm faster, especially with cross-validation. Set to FALSE to prevent any exclusion here.
#' @param max_steps Maximum number of steps taken (Default = 50).
#' @param start_latent_var Optional list of starting points for each latent variable (The list must have the same length as the number of latent variables and each element of the list must have the same length as the number of variables of the corresponding latent variable).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param seed Optional seed for bootstrap.
#' @param progress If TRUE, shows the progress done (Default=TRUE).
#' @param n_cluster Number of parallel clusters, I recommend using the number of CPU cores - 1 (Default = 1).
#' @param best_subsets If \code{best_subsets = k}, the output will show the k most frequently chosen subsets of variables (Default = 5)
#' @return Returns a list of vectors containing the percentage of times that each variable was selected within each latent variable.
#' @examples
#'	\dontrun{
#'	## Example
#'	train = example_3way_3latent(250, 2, seed=777)
#'	# Bootstrap with Bidirectional-backward search for everything based on AIC
#'	# Normally you should use a lot more than 10 iterations and extra CPUs (n_cluster)
#'	boot = bootstrap_var_select(train$data, latent_var_extra=NULL, 
#'	latent_var_original=train$latent_var,
#'	formula=y ~ E*G*Z,search_type="bidirectional-backward", search=0, 
#'	search_criterion="AIC", boot_iter=10, n_cluster=1)
#'	# Assuming it's longitudinal with 5 timepoints, even though it's not
#'	id = factor(rep(1:50,each=5))
#'	boot_longitudinal = bootstrap_var_select(train$data, latent_var_extra=NULL, 
#'	latent_var_original=train$latent_var,
#'	formula=y ~ E*G*Z,search_type="bidirectional-backward", search=0, 
#'	search_criterion="AIC", boot_iter=10, n_cluster=1, boot_group=id)
#'	}
#' @import foreach snow doSNOW utils iterators
#' @references Peter C Austin and Jack V Tu. \emph{Bootstrap Methods for Developing Predictive Models} (2012). dx.doi.org/10.1198/0003130043277.
#' @references Mark Reiser, Lanlan Yao, Xiao Wang, Jeanne Wilcox and Shelley Gray. \emph{A Comparison of Bootstrap Confidence Intervals for Multi-level Longitudinal Data Using Monte-Carlo Simulation} (2017). 10.1007/978-981-10-3307-0_17.
#' @export
"bootstrap_var_select"

#' @title Parallel genetic algorithm variable selection (for IMLEGIT)
#' @description [Very slow, recommended when the number of variables is large] Use a standard genetic algorithm with single-point crossover and a single mutation ran in parallel to find the best subset of variables. The percentage of times that each variable is included the final populations is also given. This is very computationally demanding but this finds much better solutions than either stepwise search or bootstrap variable selection.
#' @param data data.frame of the dataset to be used.
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in the formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param parallel_iter number of parallel genetic algorithms (Default = 10). I recommend using 2-4 times the number of CPU cores used.
#' @param entropy_threshold Entropy threshold for convergence of the population (Default = .10). Note that not reaching the entropy threshold just means that the population has some diversity, this is not necessarily a bad thing. Reaching the threshold is not necessary but if a population reach the threshold, we want it to stop reproducing (rather than continuing until \code{maxgen}) since the future generations won't change much.
#' @param popsize Size of the population (Default = 25). Between 25 and 100 is generally adequate.
#' @param mutation_prob Probability of mutation (Default = .50). A single variable is selected for mutation and it is mutated with probability \code{mutation_prob}. If the mutation causes a latent variable to become empty, no mutation is done. Using a small value (close to .05) will lead to getting more stuck in suboptimal solutions but using a large value (close to 1) will greatly increase the computing time because it will have a hard time reaching the entropy threshold.
#' @param first_pop optional Starting initial population which is used instead of a fully random one. Mutation is also done on the initial population to increase variability.
#' @param latent_var list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise, the latent variables will be named L1, L2, ...
#' @param search_criterion Criterion used to determine which variable is the best to add or worst to drop. If \code{search_criterion="AIC"}, uses the AIC, if \code{search_criterion="AICc"}, uses the AICc, if \code{search_criterion="BIC"}, uses the BIC, if \code{search_criterion="cv"}, uses the cross-validation error, if \cr \code{search_criterion="cv_AUC"}, uses the cross-validated AUC, if \code{search_criterion="cv_Huber"}, uses the Huber cross-validation error, if \code{search_criterion="cv_L1"}, uses the L1-norm cross-validation error (Default = "AIC"). The Huber and L1-norm cross-validation errors are alternatives to the usual cross-validation L2-norm error (which the \eqn{R^2} is based on) that are more resistant to outliers, the lower the values the better.
#' @param maxgen Maximum number of generations (iterations) of the genetic algorithm (Default = 100). Between 50 and 200 generations is generally adequate.
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results). Note that using .001 rather than .01 (default) can more than double or triple the computing time of genetic_var_select.
#' @param maxiter Maximum number of iterations.
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param seed Optional seed.
#' @param progress If TRUE, shows the progress done (Default=TRUE).
#' @param n_cluster Number of parallel clusters, I recommend using the number of CPU cores - 1 (Default = 1).
#' @param best_subsets If \code{best_subsets = k}, the output will show the k best subsets of variables (Default = 5)
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param folds Optional list of vectors containing the fold number for each observation. Bypass cv_iter and cv_folds. Setting your own folds could be important for certain data types like time series or longitudinal data.
#' @param Huber_p Parameter controlling the Huber cross-validation error (Default = 1.345).
#' @param classification Set to TRUE if you are doing classification and cross-validation (binary outcome).
#' @return Returns a list of vectors containing the percentage of times that each variable was included in the final populations, the criterion of the best k models, the starting points of the best k models (with the names of the best variables) and the entropy of the populations.
#' @examples
#'	\dontrun{
#'	## Example
#'	train = example_3way_3latent(250, 2, seed=777)
#'	# Genetic algorithm based on BIC
#'	# Normally you should use a lot more than 2 populations with 10 generations
#'	ga = genetic_var_select(train$data, latent_var=train$latent_var,
#'	formula=y ~ E*G*Z, search_criterion="AIC", parallel_iter=2, maxgen = 10)
#'	}
#' @import foreach snow doSNOW utils iterators
#' @references Mu Zhu, & Hugh Chipman. \emph{Darwinian evolution in parallel universes: A parallel genetic algorithm for variable selection} (2006). Technometrics, 48(4), 491-502.
#' @export
"genetic_var_select"

#' @title Parallel natural evolutionary variable selection assuming bernouilli distribution (for IMLEGIT)
#' @description [Slow, highly recommended when the number of variables is large] Use natural evolution strategy (nes) gradient descent ran in parallel to find the best subset of variables. It is often as good as genetic algorithms but much faster so it is the recommended variable selection function to use as default. Note that this approach assumes that the inclusion of a variable does not depends on whether other variables are included (i.e. it assumes independent bernouilli distributions); this is generally not true but this approach still converge well and running it in parallel increases the probability of reaching the global optimum.
#' @param data data.frame of the dataset to be used.
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in the formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param parallel_iter number of parallel tries (Default = 3). For speed, I recommend using the number of CPU cores.
#' @param alpha vector of the parameter for the Dirichlet distribution of the starting points (Assuming a symmetric Dirichlet distribution with only one parameter). If the vector has size N and parralel_iter=K, we use alpha[1], ..., alpha[N], alpha[1], ... , alpha[N], ... for parallel_iter 1 to K respectively. We assume a dirichlet distribution for the starting points to get a bit more variability and make sure we are not missing on a great subset of variable that doesn't converge to the global optimum with the default starting points. Use bigger values for less variability and lower values for more variability (Default = c(1,5,10)).
#' @param entropy_threshold Entropy threshold for convergence of the population (Default = .10). The smaller the entropy is, the less diversity there is in the population, which means convergence.
#' @param popsize Size of the population, the number of subsets of variables sampled at each iteration (Default = 25). Between 25 and 100 is generally adequate.
#' @param lr learning rate of the gradient descent, higher will converge faster but more likely to get stuck in local optium (Default = .2).
#' @param prop_ignored The proportion of the population that are given a fixed fitness value, thus their importance is greatly reduce. The higher it is, the longer it takes to converge. Highers values makes the algorithm focus more on favorizing the good subsets of variables than penalizing the bad subsets (Default = .50).
#' @param latent_var list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise, the latent variables will be named L1, L2, ...
#' @param search_criterion Criterion used to determine which variable subset is the best. If \code{search_criterion="AIC"}, uses the AIC, if \code{search_criterion="AICc"}, uses the AICc, if \code{search_criterion="BIC"}, uses the BIC, if \code{search_criterion="cv"}, uses the cross-validation error, if \cr \code{search_criterion="cv_AUC"}, uses the cross-validated AUC, if \code{search_criterion="cv_Huber"}, uses the Huber cross-validation error, if \code{search_criterion="cv_L1"}, uses the L1-norm cross-validation error (Default = "AIC"). The Huber and L1-norm cross-validation errors are alternatives to the usual cross-validation L2-norm error (which the \eqn{R^2} is based on) that are more resistant to outliers, the lower the values the better.
#' @param n_cluster Number of parallel clusters, I recommend using the number of CPU cores (Default = 1).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results). Note that using .001 rather than .01 (default) can more than double or triple the computing time of genetic_var_select.
#' @param maxiter Maximum number of iterations.
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param seed Optional seed.
#' @param progress If TRUE, shows the progress done (Default=TRUE).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param folds Optional list of vectors containing the fold number for each observation. Bypass cv_iter and cv_folds. Setting your own folds could be important for certain data types like time series or longitudinal data.
#' @param Huber_p Parameter controlling the Huber cross-validation error (Default = 1.345).
#' @param classification Set to TRUE if you are doing classification and cross-validation (binary outcome).
#' @param print If TRUE, print the parameters of the search distribution and the entropy at each iteration. Note: Only works using Rterm.exe in Windows due to parallel clusters. (Default = FALSE).
#' @return Returns a list containing the best subset's fit, cross-validation output, latent variables and starting points.
#' @examples
#'	\dontrun{
#'	## Example
#'	train = example_3way_3latent(250, 2, seed=777)
#'	nes = nes_var_select(train$data, latent_var=train$latent_var,
#'	formula=y ~ E*G*Z)
#'	}
#' @import foreach snow doSNOW utils iterators
#' @export
"nes_var_select"

#' @title Parallel natural evolutionary variable selection assuming multivariate normal search distribution with a simple covariance matrix parametrization (for IMLEGIT)
#' @description [Slow, highly recommended when the number of variables is large] Use natural evolution strategy (nes) gradient descent ran in parallel to find the best subset of variables. It is often as good as genetic algorithms but much faster so it is the recommended variable selection function to use as default. This is slower than nes_var_select but much less likely to get stuck into local optimum so the parallelization is not really needed.
#' @param data data.frame of the dataset to be used.
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in the formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param parallel_iter number of parallel tries (Default = 3). For speed, I recommend using the number of CPU cores.
#' @param alpha vector of the parameter for the Dirichlet distribution of the starting points (Assuming a symmetric Dirichlet distribution with only one parameter). If the vector has size N and parralel_iter=K, we use alpha[1], ..., alpha[N], alpha[1], ... , alpha[N], ... for parallel_iter 1 to K respectively. We assume a dirichlet distribution for the starting points to get a bit more variability and make sure we are not missing on a great subset of variable that doesn't converge to the global optimum with the default starting points. Use bigger values for less variability and lower values for more variability (Default = c(1,5,10)).
#' @param entropy_threshold Entropy threshold for convergence of the population (Default = .10). The smaller the entropy is, the less diversity there is in the population, which means convergence.
#' @param popsize Size of the population, the number of subsets of variables sampled at each iteration (Default = 25). Between 25 and 100 is generally adequate.
#' @param lr learning rate of the gradient descent, higher will converge faster but more likely to get stuck in local optium (Default = .2).
#' @param prop_ignored The proportion of the population that are given a fixed fitness value, thus their importance is greatly reduce. The higher it is, the longer it takes to converge. Highers values makes the algorithm focus more on favorizing the good subsets of variables than penalizing the bad subsets (Default = .50).
#' @param latent_var list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise, the latent variables will be named L1, L2, ...
#' @param search_criterion Criterion used to determine which variable subset is the best. If \code{search_criterion="AIC"}, uses the AIC, if \code{search_criterion="AICc"}, uses the AICc, if \code{search_criterion="BIC"}, uses the BIC, if \code{search_criterion="cv"}, uses the cross-validation error, if \cr \code{search_criterion="cv_AUC"}, uses the cross-validated AUC, if \code{search_criterion="cv_Huber"}, uses the Huber cross-validation error, if \code{search_criterion="cv_L1"}, uses the L1-norm cross-validation error (Default = "AIC"). The Huber and L1-norm cross-validation errors are alternatives to the usual cross-validation L2-norm error (which the \eqn{R^2} is based on) that are more resistant to outliers, the lower the values the better.
#' @param n_cluster Number of parallel clusters, I recommend using the number of CPU cores (Default = 1).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results). Note that using .001 rather than .01 (default) can more than double or triple the computing time of genetic_var_select.
#' @param maxiter Maximum number of iterations.
#' @param ylim Optional vector containing the known min and max of the outcome variable. Even if your outcome is known to be in [a,b], if you assume a Gaussian distribution, predict() could return values outside this range. This parameter ensures that this never happens. This is not necessary with a distribution that already assumes the proper range (ex: [0,1] with binomial distribution).
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param seed Optional seed.
#' @param progress If TRUE, shows the progress done (Default=TRUE).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param folds Optional list of vectors containing the fold number for each observation. Bypass cv_iter and cv_folds. Setting your own folds could be important for certain data types like time series or longitudinal data.
#' @param Huber_p Parameter controlling the Huber cross-validation error (Default = 1.345).
#' @param classification Set to TRUE if you are doing classification and cross-validation (binary outcome).
#' @param print If TRUE, print the parameters of the search distribution and the entropy at each iteration. Note: Only works using Rterm.exe in Windows due to parallel clusters. (Default = FALSE).
#' @return Returns a list containing the best subset's fit, cross-validation output, latent variables and starting points.
#' @examples
#'	\dontrun{
#'	## Example
#'	train = example_3way_3latent(250, 2, seed=777)
#'	nes = r1nes_var_select(train$data, latent_var=train$latent_var,
#'	formula=y ~ E*G*Z)
#'	}
#' @import foreach snow doSNOW utils iterators
#' @export
"r1nes_var_select"

example_2way = function(N, sigma=1, logit=FALSE, seed=NULL){
	set.seed(seed)
	g1 = rbinom(N,1,.30)
	g2 = rbinom(N,1,.30)
	g3 = rbinom(N,1,.30)
	g4 = rbinom(N,1,.30)
	g1_g3 = g1*g3
	g2_g3 = g2*g3
	e1 = rnorm(N,0,1.5)
	e2 = rnorm(N,0,1.5)
	e3 = rnorm(N,0,1.5)
	g = (.2*g1 + .15*g2 - .3*g3 + .1*g4 + .05*g1_g3 + .2*g2_g3)
	e = -.45*e1 + .35*e2 + .2*e3
	y_true = -1 + 2*g + 3*e + 4*g*e
	if (logit){
		y_true = 1/(1+exp(-(y_true)))
		y = rbinom(N,1,y_true)
	}
	else{
		eps = rnorm(N,0,sigma)
		y = y_true + eps
	}
	return(list(data=data.frame(y,y_true),G=data.frame(g1,g2,g3,g4,g1_g3,g2_g3),E=data.frame(e1,e2,e3),coef_G=c(.2,.15,-.3,.1,.05,.2),coef_E=c(-.45,.35,.2), coef_main=c(-1,2,3,4)))
}

example_with_crossover = function(N, sigma=1, c = 0, coef_main=c(0,1,2), coef_G=c(.30, .10, .20, .40), coef_E=c(.45,.35,.2), logit=FALSE, seed=NULL, beta_param=c(2,2)){
	set.seed(seed)
	g1 = rbinom(N,1,.30)
	g2 = rbinom(N,1,.30)
	g3 = rbinom(N,1,.30)
	g4 = rbinom(N,1,.30)
	e1 = rbeta(N,beta_param[1],beta_param[2])*10
	e2 = rbeta(N,beta_param[1],beta_param[2])*10
	e3 = rbeta(N,beta_param[1],beta_param[2])*10
	g = coef_G[1]*g1 + coef_G[2]*g2 + coef_G[3]*g3 + coef_G[4]*g4
	e = coef_E[1]*e1 + coef_E[2]*e2 + coef_E[3]*e3
	y_true = coef_main[1] + coef_main[2]*(e-c) + coef_main[3]*g*(e-c)
	if (logit){
		y_true = 1/(1+exp(-(y_true)))
		y = rbinom(N,1,y_true)
	}
	else{
		eps = rnorm(N,0,sigma)
		y = y_true + eps
	}
	return(list(data=data.frame(y,y_true),G=data.frame(g1,g2,g3,g4),E=data.frame(e1,e2,e3),coef_G=coef_G,coef_E=coef_E, coef_main=coef_main, c=c))
}

example_3way = function(N, sigma=2.5, logit=FALSE, seed=NULL){
	set.seed(seed)
	g1 = rbinom(N,1,.30)
	g2 = rbinom(N,1,.30)
	g3 = rbinom(N,1,.30)
	g4 = rbinom(N,1,.30)
	g1_g3 = g1*g3
	g2_g3 = g2*g3
	e1 = rnorm(N,0,1.5)
	e2 = rnorm(N,0,1.5)
	e3 = rnorm(N,0,1.5)
	g = (.2*g1 + .15*g2 - .3*g3 + .1*g4 + .05*g1_g3 + .2*g2_g3)
	e = -.45*e1 + .35*e2 + .2*e3
	z = rnorm(N,3,1)
	y_true = -2 + 2*g + 3*e + z + 5*g*e - 1.5*z*e + 2*z*g + 2*z*g*e
	if (logit){
		y_true = 1/(1+exp(-(y_true)))
		y = rbinom(N,1,y_true)
	}
	else{
		eps = rnorm(N,0,sigma)
		y = y_true + eps
	}
	return(list(data=data.frame(y,y_true,z),G=data.frame(g1,g2,g3,g4,g1_g3,g2_g3),E=data.frame(e1,e2,e3),coef_G=c(.2,.15,-.3,.1,.05,.2),coef_E=c(-.45,.35,.2), coef_main=c(-2,2,3,1,5,-1.5,2,2)))
}

example_3way_3latent = function(N, sigma=1, logit=FALSE, seed=NULL){
	set.seed(seed)
	g1 = rbinom(N,1,.30)
	g2 = rbinom(N,1,.30)
	g3 = rbinom(N,1,.30)
	g4 = rbinom(N,1,.30)
	g1_g3 = g1*g3
	g2_g3 = g2*g3
	e1 = rnorm(N,0,1.5)
	e2 = rnorm(N,0,1.5)
	e3 = rnorm(N,0,1.5)
	z1 = rnorm(N,3,1)
	z2 = rnorm(N,3,1)
	z3 = rnorm(N,3,1)
	g = (.2*g1 + .15*g2 - .3*g3 + .1*g4 + .05*g1_g3 + .2*g2_g3)
	e = -.45*e1 + .35*e2 + .2*e3
	z = .15*z1 + .75*z2 + .10*z3
	y_true = -2 + 2*g + 3*e + z + 5*g*e - 1.5*z*e + 2*z*g + 2*z*g*e
	if (logit){
		y_true = 1/(1+exp(-(y_true)))
		y = rbinom(N,1,y_true)
	}
	else{
		eps = rnorm(N,0,sigma)
		y = y_true + eps
	}
	return(list(data=data.frame(y,y_true),latent_var=list(G=data.frame(g1,g2,g3,g4,g1_g3,g2_g3),E=data.frame(e1,e2,e3),Z=data.frame(z1,z2,z3)),coef_G=c(.2,.15,-.3,.1,.05,.2),coef_E=c(-.45,.35,.2),coef_Z=c(.15,.75,.10), coef_main=c(-2,2,3,1,5,-1.5,2,2)))
}

longitudinal_folds = function(cv_iter=1, cv_folds=10, id, formula=NULL, data=NULL, data_needed=NULL, print=TRUE){
	if (cv_folds > length(unique(id))) stop("cv_folds must be smaller than the number of unique id")
	# in IMLEGIT, data_needed would be latent_var which is a list and we need to unlist it if that's the case
	if (!is.null(data_needed)) if(class(data_needed)=="list") data_needed = do.call(cbind.data.frame, data_needed)
	if (!is.null(data) && !is.null(formula)){
		# Extracting only the variables available from the formula
		formula = as.formula(formula)
		formula_full = stats::terms(formula,simplify=TRUE)
		formula_outcome = get.vars(formula)[1]
		formula_elem_ = attributes(formula_full)$term.labels
		vars_names = get.vars(formula)[get.vars(formula) %in% names(data)]
		if (!is.null(data_needed)){
			vars_names = c(vars_names,names(data_needed))
			data = data.frame(data,data_needed)
		}
		vars_names[-length(vars_names)] = paste0(vars_names[-length(vars_names)], " + ")
		formula_n = paste0(formula_outcome, " ~ ", paste0(vars_names,collapse=""))

		data = stats::model.frame(formula_n, data, na.action=na.pass)
		id = id[stats::complete.cases(data)]
	}
	else{
		if (!is.null(data_needed)) id = id[stats::complete.cases(data_needed)]
	}
	if(print) print(paste0("You have ",length(unique(id))," unique observations with ", max(table(factor(id)))," categories/time-points for a total of ", length(id)," observations"))
	folds = vector("list", cv_iter)
	for (i in 1:cv_iter){
		s = sample(sort(unique(id)))
	 	id_new = cut(1:length(s),breaks=cv_folds,labels=FALSE)
	 	folds[[i]] = rep(NA, length(id))
	 	for (j in 1:cv_folds){
	 		folds[[i]][id %in% s[id_new==j]] = j
	 	}
	}
 	return(folds)
}

LEGIT = function(data, genes, env, formula, start_genes=NULL, start_env=NULL, eps=.001, maxiter=100, family=gaussian, ylim=NULL, print=TRUE, print_steps=FALSE, crossover = NULL, crossover_fixed = FALSE, reverse_code=FALSE, rescale=FALSE)
{
	if (!is.null(ylim)){
		if (!is.numeric(ylim) || length(ylim) !=2) stop("ylim must either be NULL or a numeric vector of size two")
	}
	if (maxiter <= 0) warning("maxiter must be > 0")
	if(!is.null(start_genes)){
		if (NCOL(genes)!=length(start_genes)) stop("start_genes must either be NULL or have the same length as the number of genes")
		}
	if(!is.null(start_env)){
		if (NCOL(env)!=length(start_env)) stop("start_env must either be NULL or have the same length as the number of environments")
	}
	if (class(data) != "data.frame" && class(data) != "matrix") stop("data must be a data.frame")

	# getting right formats
	# Retaining only the needed variables from the dataset (need to set G and E variables for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	data$G=0
	data$E=0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	genes = as.matrix(genes, drop=FALSE)
	if (is.null(colnames(genes))){
		if (print) cat("You have not specified column names for genes, they will be named gene1, gene2, ...\n")
		colnames(genes) = paste0("gene",1:NCOL(genes))
	}
	env_ = env
	if (!is.null(crossover) && !crossover_fixed) env = cbind(-1,env)
	env_ = as.matrix(env_, drop=FALSE)
	env = as.matrix(env, drop=FALSE)
	if (is.null(colnames(env))){
		if (print) cat("You have not specified column names for env, they will be named env1, env2, ...\n")
		colnames(env) = paste0("env",1:NCOL(env))
		colnames(env_) = paste0("env_",1:NCOL(env_))
	}
	else if (!is.null(crossover) && !crossover_fixed) colnames(env)[1]="crossover"
	formula = stats::as.formula(formula)

	# Can only reverse genes in [0,1]
	if (reverse_code && (min(genes) !=0 || max(genes) !=1)) stop("Please make sure that genes are in range [0,1].")
	else{
		genes_names_orig = colnames(genes)
		genes_inverted = rep(FALSE, NCOL(genes))
	}

	# Error message about factors
	if (sum(apply(data,2,is.numeric)) != NCOL(data) || sum(apply(genes,2,is.numeric)) != NCOL(genes) || sum(apply(env,2,is.numeric)) != NCOL(env)) stop("All variables used must be numeric, factors are not allowed. Please dummy code all categorical variables inside your datasets (data, gene, env)")

	# remove missing data
	comp = stats::complete.cases(data,genes,env)
	data = data[comp,, drop=FALSE]
	genes = genes[comp,, drop=FALSE]
	env = env[comp,, drop=FALSE]
	env_ = env_[comp,, drop=FALSE]
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	# Rescale environments
	if (rescale){
		if (!is.null(crossover) && !crossover_fixed) env[,-1] = apply(env[,-1, drop=FALSE], 2, function(x) (x - min(x))*2/(max(x)-min(x))-1)
		else env = apply(env, 2, function(x) (x - min(x))*2/(max(x)-min(x))-1)
		env_ = apply(env_, 2, function(x) (x - min(x))*2/(max(x)-min(x))-1)
	}

	#Adding empty variables in main dataset for genes and env
	data[,colnames(genes)]=0
	data[,colnames(env)]=0
	data$R0_b=0
	data$R0_c=0

	# Setting up initial weighted scores
	if (is.null(start_genes)) weights_genes = rep(1/dim(genes)[2],dim(genes)[2])
	else if (sum(abs(start_genes))==0) weights_genes = rep(1/dim(genes)[2],dim(genes)[2])
	else weights_genes = start_genes/sum(abs(start_genes))

	if (is.null(start_env)) weights_env = rep(1/dim(env_)[2],dim(env_)[2])
	else if (sum(abs(start_env))==0) weights_env = rep(1/dim(env_)[2],dim(env_)[2])
	else weights_env = start_env/sum(abs(start_env))
	if (!is.null(crossover) && !crossover_fixed) weights_env = c(crossover, weights_env)

	data$G = genes%*%weights_genes
	if (!is.null(crossover) && crossover_fixed) data$E = env%*%weights_env - crossover
	else data$E = env%*%weights_env

	if (print_steps) cat(paste0("G: ", paste(round(weights_genes,2), collapse = " "), " // E: ", paste(round(weights_env,2), collapse = " "),"\n"))

	# Deconstructing formula into parts (No E or G / only E / only G / both G and E)
	formula_full = stats::terms(formula,simplify=TRUE)
	formula_outcome = get.vars(formula)[1]
	formula_elem_ = attributes(formula_full)$term.labels
	# Adding white spaces before and after to recognize a "E" as opposed to another string like "Elephant"
	formula_elem = paste("", formula_elem_,"")
	index_with_G = grepl(" G ",formula_elem, fixed=TRUE) | grepl(" G:",formula_elem, fixed=TRUE) | grepl(":G:",formula_elem, fixed=TRUE) | grepl(":G ",formula_elem, fixed=TRUE)
	index_with_E = grepl(" E ",formula_elem, fixed=TRUE) | grepl(" E:",formula_elem, fixed=TRUE) | grepl(":E:",formula_elem, fixed=TRUE) | grepl(":E ",formula_elem, fixed=TRUE)
	index_with_GE = index_with_G & index_with_E
	index_with_G = index_with_G & !index_with_GE
	index_with_E = index_with_E & !index_with_GE
	data_expanded = stats::model.matrix(formula, data=data)
	if (colnames(data_expanded)[1] == "(Intercept)"){
		formula_elem = c("1",formula_elem)
		index_with_G = c(FALSE, index_with_G)
		index_with_E = c(FALSE, index_with_E)
		index_with_GE = c(FALSE, index_with_GE)
	}
	index_without_GE = !(index_with_G | index_with_E | index_with_GE) 

	## Formulas for reparametrization in step b (estimating G)
	formula_elem_withoutG = formula_elem[index_without_GE | index_with_E]
	formula_elem_withoutG[-length(formula_elem_withoutG)] = paste0(formula_elem_withoutG[-length(formula_elem_withoutG)], " + ")
	formula_withoutG = paste0(formula_outcome, " ~ ", paste0(formula_elem_withoutG,collapse=""))
	if (formula_elem[1] != "1") formula_withoutG = paste0(formula_withoutG, " - 1")
	formula_withoutG = stats::as.formula(formula_withoutG)

	formula_elem_withG = formula_elem[index_with_G | index_with_GE]
	# Remove G elements from formula because we want (b1 + b2*E + ...)*G rather than b1*G + b2*E*G + ...
	formula_elem_withG = gsub(" G ","1",formula_elem_withG, fixed=TRUE)
	formula_elem_withG = gsub(" G:","",formula_elem_withG, fixed=TRUE)
	formula_elem_withG = gsub(":G:",":",formula_elem_withG, fixed=TRUE)
	formula_elem_withG = gsub(":G ","",formula_elem_withG, fixed=TRUE)
	formula_elem_withG[-length(formula_elem_withG)] = paste0(formula_elem_withG[-length(formula_elem_withG)], " + ")
	formula_withG = paste0(formula_outcome, " ~ ", paste0(formula_elem_withG,collapse=""))
	if (sum(grepl("1",formula_elem_withG, fixed=TRUE))==0) formula_withG = paste0(formula_withG, " - 1")
	formula_withG = stats::as.formula(formula_withG)

	## Formulas for reparametrization in step c (estimating E)
	formula_elem_withoutE = formula_elem[index_without_GE | index_with_G]
	formula_elem_withoutE[-length(formula_elem_withoutE)] = paste0(formula_elem_withoutE[-length(formula_elem_withoutE)], " + ")
	formula_withoutE = paste0(formula_outcome, " ~ ", paste0(formula_elem_withoutE,collapse=""))
	if (formula_elem[1] != "1") formula_withoutE = paste0(formula_withoutE, " - 1")
	formula_withoutE = stats::as.formula(formula_withoutE)

	formula_elem_withE = formula_elem[index_with_E | index_with_GE]
	# Remove E elements from formula because we want (b1 + b2*G + ...)*E rather than b1*E + b2*G*E + ...
	formula_elem_withE = gsub(" E ","1",formula_elem_withE, fixed=TRUE)
	formula_elem_withE = gsub(" E:","",formula_elem_withE, fixed=TRUE)
	formula_elem_withE = gsub(":E:",":",formula_elem_withE, fixed=TRUE)
	formula_elem_withE = gsub(":E ","",formula_elem_withE, fixed=TRUE)
	formula_elem_withE[-length(formula_elem_withE)] = paste0(formula_elem_withE[-length(formula_elem_withE)], " + ")
	formula_withE = paste0(formula_outcome, " ~ ", paste0(formula_elem_withE,collapse=""))
	if (sum(grepl("1",formula_elem_withE, fixed=TRUE))==0) formula_withE = paste0(formula_withE, " - 1")
	formula_withE = stats::as.formula(formula_withE)

	# Making formula for step b (estimating G)
	genes_names = colnames(genes)
	genes_names[-length(genes)] = paste0(colnames(genes)[-length(genes)], " + ")
	formula_b = paste0(formula_outcome, " ~ ", paste0(genes_names,collapse=""))
	formula_b = paste0(formula_b, " offset(R0_b) - 1")
	formula_b = stats::as.formula(formula_b)

	# Making formula for step c (estimating E)
	env_names = colnames(env)
	env_names[-length(env)] = paste0(colnames(env)[-length(env)], " + ")
	formula_c = paste0(formula_outcome, " ~ ", paste0(env_names,collapse=""))
	formula_c = paste0(formula_c, " offset(R0_c) - 1")
	formula_c = stats::as.formula(formula_c)

	weights_genes_old_old = NULL
	weights_genes_old = NULL

	for (i in 1:maxiter){

		## Step a : fit main model
		fit_a = stats::glm(formula, data=data, family=family, y=FALSE, model=FALSE)

		if (NCOL(genes)>1){
			# Reparametrizing variables for step b (estimating G)
			data_expanded_withoutG = stats::model.matrix(formula_withoutG, data=data)
			data$R0_b = data_expanded_withoutG%*%stats::coef(fit_a)[(index_without_GE | index_with_E)]
			data_expanded_withG = stats::model.matrix(formula_withG, data=data)
			R1_b = data_expanded_withG%*%stats::coef(fit_a)[(index_with_G | index_with_GE)]
			R1_b_genes = genes*as.vector(R1_b)
			data[,colnames(genes)]=R1_b_genes

			## Step b : fit model for G
			fit_b = stats::glm(formula_b, data=data, family=family, y=FALSE, model=FALSE)
			weights_genes_ = stats::coef(fit_b)

			# Reverse coding
			if (reverse_code && sum(weights_genes_ < 0) > 0){
				weights_genes_old_old = weights_genes_old
				# Invert gene coding
				genes[,weights_genes_ < 0] = 1 - genes[,weights_genes_ < 0]
				genes_inverted[weights_genes_ < 0] = !genes_inverted[weights_genes_ < 0]
				# changing names for user to know they were inverted
				colnames(genes)[(weights_genes_ < 0) & genes_inverted] = paste0(genes_names_orig[weights_genes_ < 0  & genes_inverted],"_inverted")
				colnames(genes)[(weights_genes_ < 0) & !genes_inverted] = genes_names_orig[weights_genes_ < 0 & !genes_inverted]
				# Remaking formula for step b (estimating G)
				genes_names = colnames(genes)
				genes_names[-length(genes)] = paste0(colnames(genes)[-length(genes)], " + ")
				formula_b = paste0(formula_outcome, " ~ ", paste0(genes_names,collapse=""))
				formula_b = paste0(formula_b, " offset(R0_b) - 1")
				formula_b = stats::as.formula(formula_b)
				# Refit
				R1_b_genes = genes*as.vector(R1_b)
				data[,colnames(genes)]=R1_b_genes
				fit_b = stats::glm(formula_b, data=data, family=family, y=FALSE, model=FALSE)
				weights_genes_ = stats::coef(fit_b)
			}

			# Updating G estimates and checking convergence
			weights_genes_old = weights_genes
			weights_genes = weights_genes_/sum(abs(weights_genes_))
			if (!is.null(weights_genes_old_old)){
				if(sqrt(sum((weights_genes_old_old-weights_genes)^2)) < eps){
					warning("Can't reverse code properly, returning the model with reverse_code=FALSE.")
					return(LEGIT(data=data, genes=genes, env=env, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=print, print_steps=print_steps, crossover = crossover, crossover_fixed = crossover_fixed, reverse_code=FALSE))
				}
			}

			data$G = genes%*%weights_genes
			if(sqrt(sum((weights_genes_old-weights_genes)^2)) < eps) conv_G = TRUE
			else conv_G = FALSE
		}
		else conv_G = TRUE

		if (NCOL(env)>1){
			# Reparametrizing variables for step c (estimating E)
			data_expanded_withoutE = stats::model.matrix(formula_withoutE, data=data)
			data$R0_c = data_expanded_withoutE%*%stats::coef(fit_a)[(index_without_GE | index_with_G)]
			data_expanded_withE = stats::model.matrix(formula_withE, data=data)
			R1_c = data_expanded_withE%*%stats::coef(fit_a)[(index_with_E | index_with_GE)]
			R1_c_env = env*as.vector(R1_c)
			data[,colnames(env)]=R1_c_env
			if (!is.null(crossover) && crossover_fixed) data$R0_c = data$R0_c - crossover*R1_c

			## Step c : fit model for E
			fit_c = stats::glm(formula_c, data=data, family=family, y=FALSE, model=FALSE)
			weights_env_ = stats::coef(fit_c)

			# Updating E estimates and checking convergence
			weights_env_old = weights_env
			if (!is.null(crossover) && !crossover_fixed) weights_env[-1] = weights_env_[-1]/sum(abs(weights_env_[-1]))
			else weights_env = weights_env_/sum(abs(weights_env_))
			if (!is.null(crossover) && crossover_fixed) data$E = env%*%weights_env - crossover
			else data$E = env%*%weights_env
			if(sqrt(sum((weights_env_old-weights_env)^2)) < eps) conv_E = TRUE
			else conv_E = FALSE
		}
		else conv_E = TRUE

		if (print_steps && !is.null(crossover) && crossover_fixed) cat(paste0("Main: ", paste(round(coef(fit_a),2), collapse = " "), " // G: ", paste(round(weights_genes,2), collapse = " "), " // E: ", paste(round(weights_env,2), collapse = " "), " // C: ", round(crossover,2),"\n"))
		else if (print_steps && !is.null(crossover) && !crossover_fixed) cat(paste0("Main: ", paste(round(coef(fit_a),2), collapse = " "), " // G: ", paste(round(weights_genes,2), collapse = " "), " // E: ", paste(round(weights_env[-1],2), collapse = " "), " // C: ", round(weights_env[1],2),"\n"))
		else if (print_steps) cat(paste0("Main: ", paste(round(coef(fit_a),2), collapse = " "), " // G: ", paste(round(weights_genes,2), collapse = " "), " // E: ", paste(round(weights_env,2), collapse = " "),"\n"))


		if (conv_G & conv_E) break
	}
	if (print_steps) cat("\n")

	# Rerunning last time and scaling to return as results
	fit_a = stats::glm(formula, data=data, family=family, y=FALSE, model=FALSE)

	# Reparametrizing variables for step b (estimating G)
	data_expanded_withoutG = stats::model.matrix(formula_withoutG, data=data)
	data$R0_b = data_expanded_withoutG%*%stats::coef(fit_a)[(index_without_GE | index_with_E)]
	data_expanded_withG = stats::model.matrix(formula_withG, data=data)
	R1_b = data_expanded_withG%*%stats::coef(fit_a)[(index_with_G | index_with_GE)]
	R1_b_genes = genes*as.vector(R1_b)
	data[,colnames(genes)]=R1_b_genes

	fit_b = stats::glm(formula_b, data=data, family=family, y=FALSE, model=FALSE)
	data[,colnames(genes)] = data[,colnames(genes)]*sum(abs(stats::coef(fit_b)))
	fit_b = stats::glm(formula_b, data=data, family=family, y=FALSE, model=FALSE)

	# Reparametrizing variables for step c (estimating E)
	data_expanded_withoutE = stats::model.matrix(formula_withoutE, data=data)
	data$R0_c = data_expanded_withoutE%*%stats::coef(fit_a)[(index_without_GE | index_with_G)]
	data_expanded_withE = stats::model.matrix(formula_withE, data=data)
	R1_c = data_expanded_withE%*%stats::coef(fit_a)[(index_with_E | index_with_GE)]
	R1_c_env = env*as.vector(R1_c)
	data[,colnames(env)]=R1_c_env
	if (!is.null(crossover) && crossover_fixed) data$R0_c = data$R0_c - crossover*R1_c

	fit_c = stats::glm(formula_c, data=data, family=family, y=FALSE, model=FALSE)
	if (!is.null(crossover) && !crossover_fixed) data[,colnames(env)] = data[,colnames(env)]*sum(abs(stats::coef(fit_c)[-1]))
	else data[,colnames(env)] = data[,colnames(env)]*sum(abs(stats::coef(fit_c)))
	fit_c = stats::glm(formula_c, data=data, family=family, y=FALSE, model=FALSE)

	if (!is.null(crossover) && !crossover_fixed) crossover = as.numeric(coef(fit_c)[1])

	# Make sure data make sense for user (and for plot function)
	fit_a$data[,colnames(env)] = env
	fit_a$data[,colnames(genes)] = genes
	fit_b$data[,colnames(genes)] = genes
	fit_c$data[,colnames(env)] = env
	weights_genes = stats::coef(fit_b)
	fit_a$data$G = genes%*%weights_genes
	fit_b$data$G = genes%*%weights_genes
	fit_c$data$G = genes%*%weights_genes
	if (!is.null(crossover) && !crossover_fixed) weights_env = stats::coef(fit_c)[-1]
	else weights_env = stats::coef(fit_c)
	fit_a$data$E = env_%*%weights_env
	fit_b$data$E = env_%*%weights_env
	fit_c$data$E = env_%*%weights_env

	if (!(abs((fit_a$deviance-fit_b$deviance)/fit_a$deviance))<.01 && abs(((fit_a$deviance-fit_c$deviance)/fit_c$deviance))<.01) warning("Deviance differs by more than 1% between model parts. Make sure that everything was set up properly and try increasing the number of iterations (maxiter).")

	#Change some arguments so that we get the right AIC, BIC and dispersion for the model
	true_rank = fit_a$rank + (fit_b$rank - 1) + (fit_c$rank - 1)
	# true_rank might be different from df of logLik, this is because glm() assume that variance components should be counted
	# I don't agree but we need to follow the same guideline to make sure it matches with glm()
	logLik_df =  attr(logLik(fit_a), "df") + (fit_b$rank - 1) + (fit_c$rank - 1)
	true_aic = 2*logLik_df - 2*logLik(fit_a)[1]
	true_aicc = true_aic + ((2*logLik_df*(logLik_df+1))/(stats::nobs(fit_a)-logLik_df-1))
	true_bic = log(stats::nobs(fit_a))*logLik_df - 2*logLik(fit_a)[1]
	true_df.residual = stats::nobs(fit_a) - true_rank
	true_null.deviance = fit_a$null.deviance

	# print convergences stuff;
	conv = (conv_G & conv_E)
	if (conv_G & conv_E){
		if (print) cat(paste0("Converged in ",i, " iterations\n"))
	} 
	else{
		warning(paste0("Did not reach convergence in maxiter iterations. Try increasing maxiter or make eps smaller."))
		if (!is.null(crossover) && !crossover_fixed) warning(paste0("The model is not garanteed to converge when a crossover point is estimated. Try different starting points for the crossover point."))
	}
	if (is.null(crossover)) crossover_fixed=NULL
	result = list(fit_main = fit_a, fit_genes = fit_b, fit_env = fit_c, true_model_parameters=list(AIC = true_aic, AICc = true_aicc, BIC = true_bic, rank = true_rank, df.residual = true_df.residual, null.deviance=true_null.deviance), ylim=ylim, formula=formula, crossover=crossover, crossover_fixed=crossover_fixed, conv=conv)
	class(result) <- "LEGIT"
	return(result)
}

GxE_interaction_test = function(data, genes, env, formula_noGxE, crossover=c("min","max"), include_noGxE_models = TRUE, reverse_code=FALSE, rescale = FALSE, boot = NULL, criterion="BIC", start_genes=NULL, start_env=NULL, eps=.001, maxiter=100, family=gaussian, ylim=NULL, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, id=NULL, classification=FALSE, seed=NULL)
{
	formula = stats::as.formula(formula_noGxE)
	formula_WEAK = paste0(formula, " + G*E - G")
	formula_STRONG = paste0(formula, " + G*E - G - E")
	formula_WEAK_nocrossover = paste0(formula, " + G*E")
	formula_STRONG_nocrossover = paste0(formula, " + G*E - E")

	if (include_noGxE_models){
		if (criterion == "cv" || criterion == "cv_AUC" || criterion=="cv_Huber" || criterion=="cv_L1") warning("When using include_noGxE_models=TRUE, the criterion must be one of the following: AIC, AICc, or BIC.")

		genes_formula = paste0(" + ", colnames(genes))
		env_formula = paste0(" + ", colnames(env))
		formula_model_E_only = paste0(formula, env_formula)
		formula_model_intercept_only = formula
		formula_model_GandE_only = paste0(formula, genes_formula, env_formula)
		formula_model_G_only = paste0(formula, genes_formula)

		model_E_only = glm(data=cbind(data, genes, env), formula=formula_model_E_only, family=family)
		model_intercept_only = glm(data=cbind(data, genes, env), formula=formula_model_intercept_only, family=family)
		model_GandE_only = glm(data=cbind(data, genes, env), formula=formula_model_GandE_only, family=family)
		model_G_only = glm(data=cbind(data, genes, env), formula=formula_model_G_only, family=family)
	}

	# 4 Models
	# If min or max, then we must find the min or max E and make it the crossover after
	# We also use the G and E weights too as starting point for even more stability.
	crossover_vantage_sensitivity_WEAK = crossover[1]
	crossover_vantage_sensitivity_STRONG = crossover[1]
	crossover_diathesis_stress_WEAK = crossover[2]
	crossover_diathesis_stress_STRONG = crossover[2]

	# Need to take min if coef > 0 and max if coef < 0 to get lower bound
	# Need to take max if coef > 0 and min if coef < 0 to get upper bound
	get_LU = function(fit){
		coef_env = coef(fit$fit_env)
		env_lower_bound = rep(0, NCOL(env))
		env_upper_bound = rep(0, NCOL(env))
		for (i in 1:NCOL(env)){
			if (coef_env[i] >= 0){
				env_lower_bound[i] = min(fit$fit_main$data[,names(env)[i]])
				env_upper_bound[i] = max(fit$fit_main$data[,names(env)[i]])
			}
			else{
				env_lower_bound[i] = max(fit$fit_main$data[,names(env)[i]])
				env_upper_bound[i] = min(fit$fit_main$data[,names(env)[i]])

			}
		}
		return(list(L=env_lower_bound, U=env_upper_bound))
	}

	# Vantage sensitivity c = min
	if (crossover[1] == "min"){
		vantage_sensitivity_WEAK = LEGIT(data=data, genes=genes, env=env, formula=formula_WEAK, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, reverse_code=reverse_code, rescale=rescale)
		if (crossover[1] == "min") crossover_vantage_sensitivity_WEAK = as.numeric(get_LU(vantage_sensitivity_WEAK)$L %*% coef(vantage_sensitivity_WEAK$fit_env))
		vantage_sensitivity_WEAK = LEGIT(data=data, genes=genes, env=env, formula=formula_WEAK, start_genes=coef(vantage_sensitivity_WEAK$fit_genes), start_env=coef(vantage_sensitivity_WEAK$fit_env), eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover_vantage_sensitivity_WEAK, crossover_fixed = TRUE, reverse_code=reverse_code, rescale=rescale)

		vantage_sensitivity_STRONG = LEGIT(data=data, genes=genes, env=env, formula=formula_STRONG, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, reverse_code=reverse_code, rescale=rescale)
		if (crossover[1] == "min") crossover_vantage_sensitivity_STRONG = as.numeric(get_LU(vantage_sensitivity_STRONG)$L %*% coef(vantage_sensitivity_STRONG$fit_env))
		vantage_sensitivity_STRONG = LEGIT(data=data, genes=genes, env=env, formula=formula_STRONG, start_genes=coef(vantage_sensitivity_STRONG$fit_genes), start_env=coef(vantage_sensitivity_STRONG$fit_env), eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover_vantage_sensitivity_STRONG, crossover_fixed = TRUE, reverse_code=reverse_code, rescale=rescale)
	}
	else{
		vantage_sensitivity_WEAK = LEGIT(data=data, genes=genes, env=env, formula=formula_WEAK, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover[1], crossover_fixed = TRUE, reverse_code=reverse_code, rescale=rescale)
		vantage_sensitivity_STRONG = LEGIT(data=data, genes=genes, env=env, formula=formula_STRONG, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover[1], crossover_fixed = TRUE, reverse_code=reverse_code, rescale=rescale)
	}
	# Diathesis-Stress c = max
	if (crossover[2] == "max"){
		diathesis_stress_WEAK = LEGIT(data=data, genes=genes, env=env, formula=formula_WEAK, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, reverse_code=reverse_code, rescale=rescale)
		if (crossover[2] == "max") crossover_diathesis_stress_WEAK = as.numeric(get_LU(diathesis_stress_WEAK)$U %*% coef(diathesis_stress_WEAK$fit_env))
		diathesis_stress_WEAK = LEGIT(data=data, genes=genes, env=env, formula=formula_WEAK, start_genes=coef(diathesis_stress_WEAK$fit_genes), start_env=coef(diathesis_stress_WEAK$fit_env), eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover_diathesis_stress_WEAK, crossover_fixed = TRUE, reverse_code=reverse_code, rescale=rescale)

		diathesis_stress_STRONG = LEGIT(data=data, genes=genes, env=env, formula=formula_STRONG, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, reverse_code=reverse_code, rescale=rescale)
		if (crossover[2] == "max") crossover_diathesis_stress_STRONG = as.numeric(get_LU(diathesis_stress_STRONG)$U %*% coef(diathesis_stress_STRONG$fit_env))
		diathesis_stress_STRONG = LEGIT(data=data, genes=genes, env=env, formula=formula_STRONG, start_genes=coef(diathesis_stress_STRONG$fit_genes), start_env=coef(diathesis_stress_STRONG$fit_env), eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover_diathesis_stress_STRONG, crossover_fixed = TRUE, reverse_code=reverse_code, rescale=rescale)
	}
	else{
		diathesis_stress_WEAK = LEGIT(data=data, genes=genes, env=env, formula=formula_WEAK, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover[2], crossover_fixed = TRUE, reverse_code=reverse_code, rescale=rescale)
		diathesis_stress_STRONG = LEGIT(data=data, genes=genes, env=env, formula=formula_STRONG, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover[2], crossover_fixed = TRUE, reverse_code=reverse_code, rescale=rescale)
	}

	# We first run the models as normal LEGIT models without crossover, then we use C=-beta_G/beta_GxE as the crossover starting point. 
	# We also use the G and E weights too as starting point for even more stability.
	diff_suscept_WEAK = LEGIT(data=data, genes=genes, env=env, formula=formula_WEAK_nocrossover, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, reverse_code=reverse_code, rescale=rescale)
	crossover_diff_suscept_WEAK = -coef(diff_suscept_WEAK$fit_main)[2]/coef(diff_suscept_WEAK$fit_main)[4]
	diff_suscept_WEAK = LEGIT(data=data, genes=genes, env=env, formula=formula_WEAK, start_genes=coef(diff_suscept_WEAK$fit_genes), start_env=coef(diff_suscept_WEAK$fit_env), eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover_diff_suscept_WEAK, crossover_fixed = FALSE, reverse_code=reverse_code, rescale=rescale)
	diff_suscept_STRONG = LEGIT(data=data, genes=genes, env=env, formula=formula_STRONG_nocrossover, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, reverse_code=reverse_code, rescale=rescale)
	crossover_diff_suscept_STRONG = -coef(diff_suscept_STRONG$fit_main)[2]/coef(diff_suscept_STRONG$fit_main)[3]
	diff_suscept_STRONG = LEGIT(data=data, genes=genes, env=env, formula=formula_STRONG, start_genes=coef(diff_suscept_STRONG$fit_genes), start_env=coef(diff_suscept_STRONG$fit_env), eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover_diff_suscept_STRONG, crossover_fixed = FALSE, reverse_code=reverse_code, rescale=rescale)

	# Criterions
	if (criterion == "cv" || criterion == "cv_AUC" || criterion=="cv_Huber" || criterion=="cv_L1"){

		if (crossover[1] == "min"){
			vantage_sensitivity_WEAK_cv = LEGIT_cv(data=data, genes=genes, env=env, formula=formula_WEAK, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed, id=id, crossover = crossover_vantage_sensitivity_WEAK, crossover_fixed = TRUE)
			vantage_sensitivity_STRONG_cv = LEGIT_cv(data=data, genes=genes, env=env, formula=formula_STRONG, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed, id=id, crossover = crossover_vantage_sensitivity_STRONG, crossover_fixed = TRUE)
		}
		else{
			vantage_sensitivity_WEAK_cv = LEGIT_cv(data=data, genes=genes, env=env, formula=formula_WEAK, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed, id=id, crossover = crossover[1], crossover_fixed = TRUE)
			vantage_sensitivity_STRONG_cv = LEGIT_cv(data=data, genes=genes, env=env, formula=formula_STRONG, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed, id=id, crossover = crossover[1], crossover_fixed = TRUE)
		}
		if (crossover[2] == "max"){
			diathesis_stress_WEAK_cv = LEGIT_cv(data=data, genes=genes, env=env, formula=formula_WEAK, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed, id=id, crossover = crossover_diathesis_stress_WEAK, crossover_fixed = TRUE)
			diathesis_stress_STRONG_cv = LEGIT_cv(data=data, genes=genes, env=env, formula=formula_STRONG, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed, id=id, crossover = crossover_diathesis_stress_STRONG, crossover_fixed = TRUE)
		}
		else{
			diathesis_stress_WEAK_cv = LEGIT_cv(data=data, genes=genes, env=env, formula=formula_WEAK, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed, id=id, crossover = crossover[2], crossover_fixed = TRUE)
			diathesis_stress_STRONG_cv = LEGIT_cv(data=data, genes=genes, env=env, formula=formula_STRONG, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed, id=id, crossover = crossover[2], crossover_fixed = TRUE)
		}
		diff_suscept_WEAK_cv = LEGIT_cv(data=data, genes=genes, env=env, formula=formula_WEAK, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed, id=id, crossover = crossover_diff_suscept_WEAK, crossover_fixed = FALSE)
		diff_suscept_STRONG_cv = LEGIT_cv(data=data, genes=genes, env=env, formula=formula_STRONG, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed, id=id, crossover = crossover_diff_suscept_STRONG, crossover_fixed = FALSE)

		if (criterion == "cv") model_criterion = c(mean(vantage_sensitivity_WEAK_cv$R2_cv), mean(vantage_sensitivity_STRONG_cv$R2_cv),mean(diathesis_stress_WEAK_cv$R2_cv), mean(diathesis_stress_STRONG_cv$R2_cv), mean(diff_suscept_WEAK_cv$R2_cv), mean(diff_suscept_STRONG_cv$R2_cv))
		else if (criterion == "cv_Huber") model_criterion = c(mean(vantage_sensitivity_WEAK_cv$Huber_cv), mean(vantage_sensitivity_STRONG_cv$Huber_cv),mean(diathesis_stress_WEAK_cv$Huber_cv), mean(diathesis_stress_STRONG_cv$Huber_cv), mean(diff_suscept_WEAK_cv$Huber_cv), mean(diff_suscept_STRONG_cv$Huber_cv))
		else if (criterion == "cv_L1") model_criterion = c(mean(vantage_sensitivity_WEAK_cv$L1_cv), mean(vantage_sensitivity_STRONG_cv$L1_cv),mean(diathesis_stress_WEAK_cv$L1_cv), mean(diathesis_stress_STRONG_cv$L1_cv), mean(diff_suscept_WEAK_cv$L1_cv), mean(diff_suscept_STRONG_cv$L1_cv))
		else if (criterion == "cv_AUC") model_criterion = c(mean(vantage_sensitivity_WEAK_cv$AUC), mean(vantage_sensitivity_STRONG_cv$AUC),mean(diathesis_stress_WEAK_cv$AUC), mean(diathesis_stress_STRONG_cv$AUC), mean(diff_suscept_WEAK_cv$AUC), mean(diff_suscept_STRONG_cv$AUC))
		# List in order from best to worse
		ordering = order(model_criterion, decreasing = TRUE)
	}
	else{
		if (criterion == "AIC") model_criterion = c(vantage_sensitivity_WEAK$true_model_parameters$AIC, vantage_sensitivity_STRONG$true_model_parameters$AIC,diathesis_stress_WEAK$true_model_parameters$AIC, diathesis_stress_STRONG$true_model_parameters$AIC, diff_suscept_WEAK$true_model_parameters$AIC, diff_suscept_STRONG$true_model_parameters$AIC)
		else if (criterion == "AICc") model_criterion = c(vantage_sensitivity_WEAK$true_model_parameters$AICc, vantage_sensitivity_STRONG$true_model_parameters$AICc,diathesis_stress_WEAK$true_model_parameters$AICc, diathesis_stress_STRONG$true_model_parameters$AICc, diff_suscept_WEAK$true_model_parameters$AICc, diff_suscept_STRONG$true_model_parameters$AICc)
		else if (criterion == "BIC") model_criterion = c(vantage_sensitivity_WEAK$true_model_parameters$BIC, vantage_sensitivity_STRONG$true_model_parameters$BIC,diathesis_stress_WEAK$true_model_parameters$BIC, diathesis_stress_STRONG$true_model_parameters$BIC, diff_suscept_WEAK$true_model_parameters$BIC, diff_suscept_STRONG$true_model_parameters$BIC)
		if (include_noGxE_models){
			logLik_df1 =  attr(logLik(model_E_only), "df")
			logLik_df2 =  attr(logLik(model_intercept_only), "df")
			logLik_df3 =  attr(logLik(model_GandE_only), "df")
			logLik_df4 =  attr(logLik(model_G_only), "df")
			if (criterion == "AIC") model_criterion = c(model_criterion, 2*logLik_df1 - 2*logLik(model_E_only)[1], 2*logLik_df2 - 2*logLik(model_intercept_only)[1], 2*logLik_df3 - 2*logLik(model_GandE_only)[1], 2*logLik_df4 - 2*logLik(model_G_only)[1])
			else if (criterion == "AICc") model_criterion = c(model_criterion, (2*logLik_df1 - 2*logLik(model_E_only)[1]) + ((2*logLik_df1*(logLik_df1+1))/(stats::nobs(model_E_only)-logLik_df1-1)), (2*logLik_df2 - 2*logLik(model_intercept_only)[1]) + ((2*logLik_df2*(logLik_df2+1))/(stats::nobs(model_intercept_only)-logLik_df2-1)), (2*logLik_df3 - 2*logLik(model_GandE_only)[1]) + ((2*logLik_df3*(logLik_df3+1))/(stats::nobs(model_GandE_only)-logLik_df3-1)),(2*logLik_df4 - 2*logLik(model_G_only)[1]) + ((2*logLik_df4*(logLik_df4+1))/(stats::nobs(model_G_only)-logLik_df4-1)))
			else if (criterion == "BIC") model_criterion = c(model_criterion, log(stats::nobs(model_E_only))*logLik_df1 - 2*logLik(model_E_only)[1], log(stats::nobs(model_intercept_only))*logLik_df2 - 2*logLik(model_intercept_only)[1], log(stats::nobs(model_GandE_only))*logLik_df3 - 2*logLik(model_GandE_only)[1], log(stats::nobs(model_G_only))*logLik_df4 - 2*logLik(model_G_only)[1])
		}
		# List in order from best to worse
		ordering = order(model_criterion)
	}
	model_criterion = round(model_criterion[ordering],2)
	if (include_noGxE_models) crossover_models = round(c(vantage_sensitivity_WEAK$crossover, vantage_sensitivity_STRONG$crossover, diathesis_stress_WEAK$crossover, diathesis_stress_STRONG$crossover, diff_suscept_WEAK$crossover, diff_suscept_STRONG$crossover, NA, NA, NA, NA)[ordering],2)
	else crossover_models = round(c(vantage_sensitivity_WEAK$crossover, vantage_sensitivity_STRONG$crossover, diathesis_stress_WEAK$crossover, diathesis_stress_STRONG$crossover, diff_suscept_WEAK$crossover, diff_suscept_STRONG$crossover)[ordering],2)
	# 95% interval of crossover for differential susceptibility models
	if (!is.null(boot)){
		# We must do bootstrap
		LEGIT_boot = function(d, i){
			data_ = d[i,,drop=FALSE]
			genes_ = genes[i,,drop=FALSE]
			env_ = env[i,,drop=FALSE]
			diff_suscept_WEAK = LEGIT(data=data_, genes=genes_, env=env_, formula=formula_WEAK_nocrossover, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, reverse_code=reverse_code, rescale=rescale)
			diff_suscept_STRONG = LEGIT(data=data_, genes=genes_, env=env_, formula=formula_STRONG_nocrossover, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, reverse_code=reverse_code, rescale=rescale)
			return (c(-coef(diff_suscept_WEAK$fit_main)[2]/coef(diff_suscept_WEAK$fit_main)[4],-coef(diff_suscept_STRONG$fit_main)[2]/coef(diff_suscept_STRONG$fit_main)[3]))
		}
		boot_results = boot::boot(data, LEGIT_boot, boot)
		crossover_interval_WEAK = round(boot::boot.ci(boot_results, index=1, type="bca")$bca[4:5],2)
		crossover_interval_STRONG = round(boot::boot.ci(boot_results, index=2, type="bca")$bca[4:5],2)
	}
	else{
		# Non-bootstrapped
		crossover_interval_WEAK = round(suppressMessages(stats::confint(diff_suscept_WEAK$fit_env,"crossover")),2)
		crossover_interval_STRONG = round(suppressMessages(stats::confint(diff_suscept_STRONG$fit_env,"crossover")),2)
	}
	data = diff_suscept_WEAK$fit_main$data
	inside_WEAK = (crossover_interval_WEAK[1] > min(data$E) && crossover_interval_WEAK[2] < max(data$E))
	inside_STRONG = (crossover_interval_STRONG[1] > min(data$E) && crossover_interval_STRONG[2] < max(data$E))
	crossover_interval_WEAK = paste0("( ",crossover_interval_WEAK[1]," / ",crossover_interval_WEAK[2]," )")
	crossover_interval_STRONG = paste0("( ",crossover_interval_STRONG[1]," / ",crossover_interval_STRONG[2]," )")

	# Proportion of observations below crossover point
	if (include_noGxE_models) prop_below = c(sum(vantage_sensitivity_WEAK$fit_main$data$E <= crossover_vantage_sensitivity_WEAK)/NROW(vantage_sensitivity_WEAK$fit_main$data), sum(vantage_sensitivity_STRONG$fit_main$data$E <= crossover_vantage_sensitivity_STRONG)/NROW(vantage_sensitivity_STRONG$fit_main$data), sum(diathesis_stress_WEAK$fit_main$data$E <= crossover_diathesis_stress_WEAK)/NROW(diathesis_stress_WEAK$fit_main$data), sum(diathesis_stress_STRONG$fit_main$data$E <= crossover_diathesis_stress_STRONG)/NROW(diathesis_stress_STRONG$fit_main$data), sum(diff_suscept_WEAK$fit_main$data$E <= crossover_diff_suscept_WEAK)/NROW(diff_suscept_WEAK$fit_main$data), sum(diff_suscept_STRONG$fit_main$data$E <= crossover_diff_suscept_STRONG)/NROW(diff_suscept_STRONG$fit_main$data), NA, NA, NA, NA)
	else prop_below = c(sum(vantage_sensitivity_WEAK$fit_main$data$E <= crossover_vantage_sensitivity_WEAK)/NROW(vantage_sensitivity_WEAK$fit_main$data), sum(vantage_sensitivity_STRONG$fit_main$data$E <= crossover_vantage_sensitivity_STRONG)/NROW(vantage_sensitivity_STRONG$fit_main$data), sum(diathesis_stress_WEAK$fit_main$data$E <= crossover_diathesis_stress_WEAK)/NROW(diathesis_stress_WEAK$fit_main$data), sum(diathesis_stress_STRONG$fit_main$data$E <= crossover_diathesis_stress_STRONG)/NROW(diathesis_stress_STRONG$fit_main$data), sum(diff_suscept_WEAK$fit_main$data$E <= crossover_diff_suscept_WEAK)/NROW(diff_suscept_WEAK$fit_main$data), sum(diff_suscept_STRONG$fit_main$data$E <= crossover_diff_suscept_STRONG)/NROW(diff_suscept_STRONG$fit_main$data))

	# Aggregating results
	if (include_noGxE_models){
		fits = list(vantage_sensitivity_WEAK = vantage_sensitivity_WEAK, vantage_sensitivity_STRONG = vantage_sensitivity_STRONG, diathesis_stress_WEAK = diathesis_stress_WEAK, diathesis_stress_STRONG = diathesis_stress_STRONG, diff_suscept_WEAK = diff_suscept_WEAK, diff_suscept_STRONG = diff_suscept_STRONG, model_E_only = model_E_only , model_intercept_only = model_intercept_only, model_GandE_only = model_GandE_only, model_G_only = model_G_only)[ordering]
		results = cbind(model_criterion, crossover_models,cbind("","","","",crossover_interval_WEAK,crossover_interval_STRONG,"","","","")[ordering],cbind("","","","",c("No","Yes")[inside_WEAK+1],c("No","Yes")[inside_STRONG+1],"","","","")[ordering], prop_below[ordering])
		rownames(results) = c("Vantage sensitivity WEAK","Vantage sensitivity STRONG","Diathesis-stress WEAK","Diathesis-stress STRONG","Differential susceptibility WEAK","Differential susceptibility STRONG", "E only", "Intercept only", "G + E only", "G only")[ordering]
	} 
	else{
		fits = list(vantage_sensitivity_WEAK = vantage_sensitivity_WEAK, vantage_sensitivity_STRONG = vantage_sensitivity_STRONG, diathesis_stress_WEAK = diathesis_stress_WEAK, diathesis_stress_STRONG = diathesis_stress_STRONG, diff_suscept_WEAK = diff_suscept_WEAK, diff_suscept_STRONG = diff_suscept_STRONG)[ordering]
		results = cbind(model_criterion, crossover_models,cbind("","","","",crossover_interval_WEAK,crossover_interval_STRONG)[ordering],cbind("","","","",c("No","Yes")[inside_WEAK+1],c("No","Yes")[inside_STRONG+1])[ordering], prop_below[ordering])
		rownames(results) = c("Vantage sensitivity WEAK","Vantage sensitivity STRONG","Diathesis-stress WEAK","Diathesis-stress STRONG","Differential susceptibility WEAK","Differential susceptibility STRONG")[ordering]
	}
	colnames(results)[1] = criterion
	colnames(results)[2] = "crossover"
	colnames(results)[3] = "crossover 95%"
	colnames(results)[4] = "Within observable range?"
	colnames(results)[5] = "% of observations below crossover"
	return(list(fits = fits, results = results, E_range=c(min(data$E),max(data$E))))
}

GxE_interaction_RoS = function(data, genes, env, formula_noGxE, t_alpha = .05, start_genes=NULL, start_env=NULL, eps=.001, maxiter=100, ylim=NULL, reverse_code=FALSE, rescale=FALSE){
	formula = stats::as.formula(formula_noGxE)
	formula = paste0(formula, " + G*E")
	fit_LEGIT = LEGIT(data=data, genes=genes, env=env, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=gaussian, ylim=ylim, print=FALSE, print_steps=FALSE, crossover = NULL, crossover_fixed = FALSE, reverse_code=reverse_code, rescale=rescale)
	fit = lm(formula, fit_LEGIT$fit_main$data)

	b1 = coef(fit)["G"]
	b3 = coef(fit)["G:E"]

	b1_var = vcov(fit)["G","G"]
	b3_var = vcov(fit)["G:E","G:E"]
	b1_b3_cov = vcov(fit)["G","G:E"]
  	t_crit = qt(t_alpha/2, fit$df.residual, lower.tail = FALSE)

  	a = (t_crit^2)*b3_var - (b3^2)
	b = 2*((t_crit^2)*b1_b3_cov - b1*b3)
	c = (t_crit^2)*b1_var - (b1^2)	
	if ((b^2) - 4*a*c < 0) return(list(RoS=c(NA,NA), int_type="Undetermined"))
	RoS_results = c((-b + sqrt((b^2) - 4*a*c))/(2*a), (-b - sqrt((b^2) - 4*a*c))/(2*a))
	RoS_results = sort(RoS_results)
	names(RoS_results)=NULL
	Lower_bounded = TRUE
	Upper_bounded = TRUE
	if (RoS_results[1] < min(fit$model$E)) Lower_bounded = FALSE
	if (RoS_results[2] > max(fit$model$E)) Upper_bounded = FALSE
	if (Lower_bounded && Upper_bounded) int_type = "Differential susceptibility"
	else if (!Lower_bounded && Upper_bounded) int_type = "Vantage sensitivity"
	else if (Lower_bounded && !Upper_bounded) int_type = "Diathesis-stress"
	else int_type = "Undetermined"
	return(list(RoS=RoS_results, int_type=int_type))
}

plot.LEGIT = function(x, cov_values = NULL, gene_quant = c(.025,.50,.975), env_quant = c(.025,.50,.975), outcome_quant = c(.025,.50,.975), cols = c("#3288BD", "#CAB176", "#D53E4F"), ylab="Outcome", xlab="Environment", legtitle="Genetic score", leglab=NULL, xlim= NULL, ylim= NULL, x_at = NULL, y_at = NULL, cex.axis = 1.9, cex.lab=2, cex.main=2.2, cex.leg=2.2, legend="topleft", ...){
	
	# Better names (need to use x for S3 class consistency)
	object = x
	formula = object$formula

	# Checks

	# Quantile range
	gene_range = as.numeric(quantile(object$fit_main$data$G, gene_quant))
	env_range = as.numeric(quantile(object$fit_main$data$E, env_quant))
	formula_outcome = get.vars(formula)[1]
	outcome_range = as.numeric(quantile(object$fit_main$data[formula_outcome][,], outcome_quant))
	# Fix this attribute to prevent problems if trying to predict fit_main directly later on
	attr(object$fit_main$terms,"dataClasses")[attr(object$fit_main$terms,"dataClasses")=="nmatrix.1"] = "numeric"

	# Defaults
	if (is.null(ylim)){
		if (is.null(y_at)) ylim = c(outcome_range[1], outcome_range[length(outcome_range)])
		else ylim = c(y_at[1], y_at[length(y_at)])
	}
	else if (is.null(y_at)) y_at = seq(ylim[1], ylim[2], length.out=3) 
	if (is.null(xlim)){
		if (is.null(x_at)) xlim = c(env_range[1], env_range[length(env_range)])
		else xlim = c(x_at[1], x_at[length(x_at)])
	}
	else if (is.null(x_at)) x_at = seq(xlim[1], xlim[2], length.out=3) 

	# Plot
	op <- par(mfrow=c(1,1), mar=c(5.1, 5.1, 5.1, 4.1))
	graphics::plot(x=c(), y=c() ,ylab = ylab, xlab = xlab, ylim = ylim, xlim=xlim, cex.lab=cex.lab, cex.main=cex.main, pch=20, bty="n", xaxt="n", yaxt="n")
	if (is.null(y_at)) graphics::axis(2, at=outcome_range, labels=paste0(outcome_quant*100,"%"), cex.axis = cex.axis)
	else graphics::axis(2, at=y_at, cex.axis = cex.axis)
	if (is.null(x_at)) graphics::axis(1, at=env_range, labels=paste0(env_quant*100,"%"), cex.axis = cex.axis)
	else graphics::axis(1, at=x_at, cex.axis = cex.axis)
	E = seq(xlim[1],xlim[2], length.out=101)
	if (!is.null(object$crossover)) E_ = E - object$crossover
	else E_ = E

	# covariates
	covariates = get.vars(formula)[-1]
	covariates = covariates[covariates != "G" & covariates != "E"]
	if (!is.null(cov_values) && length(cov_values)!=length(covariates)) stop("cov_values doesn't have the correct number of covariates")

	# Predictions and plot
	for (j in 1:length(gene_range)){
		# making data for predictions
		G = gene_range[j]
		newdata_base = data.frame(E=E_, G=G)
		newdata = newdata_base
		for (i in 1:length(covariates)){
			if (identical(covariates, character(0))) newdata = newdata
			else if (is.null(cov_values)){
				newdata = cbind(newdata, rep(apply(object$fit_main$data[covariates[i]], 2, mean),101))
				colnames(newdata)[NCOL(newdata)] = covariates[i]
			}
			else{
				newdata = cbind(newdata, rep(cov_values[i], 101))
				colnames(newdata)[NCOL(newdata)] = names(cov_values)[i]
			}
		}
		# Prediction lines
		preds <- predict(object$fit_main, newdata = newdata, se.fit = TRUE, type = "link")
		ilink <- family(object$fit_main)$linkinv
		graphics::polygon(c(E,rev(E)),c(ilink(preds$fit-2*preds$se.fit),rev(ilink(preds$fit+2*preds$se.fit))),col=grDevices::adjustcolor(cols[j], alpha.f = 0.2),border = NA)
		graphics::lines(E,ilink(preds$fit), col=cols[j], lwd=2)
	}
	if (is.null(leglab)) leglab = paste0(gene_quant*100,"%")
	legend(legend, legend=leglab,col = cols, lty=1, lwd=3, xpd = TRUE, cex = cex.leg, title=legtitle)
}

IMLEGIT = function(data, latent_var, formula, start_latent_var=NULL, eps=.001, maxiter=100, family=gaussian, ylim=NULL, print=TRUE)
{
	if (!is.null(ylim)){
		if (!is.numeric(ylim) || length(ylim) !=2) stop("ylim must either be NULL or a numeric vector of size two")
	}
	# Setting up latent_var and checks
	if (class(latent_var)!="list") stop("latent_var must be a list of datasets")
	k = length(latent_var)
	if (k==0) stop("latent_var cannot be an empty list")
	if (is.null(names(latent_var))){
		if (print) cat("You have not specified names for the latent variables, assigning names to latent_var is highly recommended to prevent confusion. For now, they will be named L1, L2, ...\n")
		names(latent_var) = paste0("L",1:k)
	}
	for (i in 1:k){
		latent_var[[i]] = as.matrix(data.frame(latent_var[[i]],fix.empty.names=FALSE))
		if (sum(colnames(latent_var[[i]])=="") > 0){
			if (print) cat(paste0("You have not specified column names for certain elements in ",names(latent_var)[i], ", elements of this latent variable will be named ",names(latent_var)[i],1,", ",names(latent_var)[i],2," ...\n"))
			colnames(latent_var[[i]]) = paste0(names(latent_var)[i],1:NCOL(latent_var[[i]]))
		}
	}

	# More checks
	if (maxiter <= 0) warning("maxiter must be > 0")
	if (k > 1) for (i in 1:(k-1)) if (NROW(latent_var[[i]]) != NROW(latent_var[[i+1]])) stop("Some datasets in latent_var don't have the same number of observations")
	if(!is.null(start_latent_var)){
		if (class(start_latent_var)!="list") stop("start_latent_var must be a lit of vectors (or NULL)")
		if (k!=length(start_latent_var)) stop("start_latent_var must have the same size as latent_var")
		for (i in 1:k){
			if (!is.null(latent_var[[i]])){
				if (NCOL(latent_var[[i]])!=length(start_latent_var[[i]])) stop("All elements of start_latent_var must either be NULL or have the same length as the number of the elements in its associated latent variable")
			}
		}
	}
	if (class(data) != "data.frame" && class(data) != "matrix") stop("data must be a data.frame")

	# getting right formats
	# Retaining only the needed variables from the dataset (need to set elements in latent_var for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	for (i in 1:k) data[,names(latent_var)[i]] = 0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	formula = stats::as.formula(formula)

	# Error message about factors
	if (sum(apply(data,2,is.numeric)) != NCOL(data)) stop("All variables used must be numeric, factors are not allowed. Please dummy code all categorical variables inside your datasets (data, latent_var[[1]], latent_var[[2]], ...)")
	for (i in 1:k) if (sum(apply(latent_var[[i]],2,is.numeric)) != NCOL(latent_var[[i]])) stop("All variables used must be numeric, factors are not allowed. Please dummy code all categorical variables inside your datasets (data, latent_var[[1]], latent_var[[2]], ...)")

	# remove missing data
	comp = stats::complete.cases(data,latent_var[[1]])
	if (k > 1) for (i in 2:k) comp = comp & stats::complete.cases(latent_var[[i]])
	data = data[comp,, drop=FALSE]
	for (i in 1:k) latent_var[[i]] = latent_var[[i]][comp,, drop=FALSE]
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	#Adding empty variables in main dataset for latent_var
	for (i in 1:k){
		data[,colnames(latent_var[[i]])]=0
		data[,paste0("R0_",i)]=0
	}

	# Setting up initial weighted latent_var
	weights_latent_var_old = vector("list", k)
	weights_latent_var = vector("list", k)
	if (is.null(start_latent_var)){
		for (i in 1:k) weights_latent_var[[i]] = rep(1/dim(latent_var[[i]])[2],dim(latent_var[[i]])[2])
	}
	else{
		for (i in 1:k){
			if (sum(abs(start_latent_var[[i]]))==0) weights_latent_var[[i]] = rep(1/dim(latent_var[[i]])[2],dim(latent_var[[i]])[2])
			else weights_latent_var[[i]] = start_latent_var[[i]]/sum(abs(start_latent_var[[i]]))
		}		
	}
	for (i in 1:k) data[,names(latent_var)[i]] = latent_var[[i]]%*%weights_latent_var[[i]]

	# Lists needed for later
	index_with_latent_var = vector("list", k)
	formula_withoutlatent_var  = vector("list", k)
	formula_withlatent_var  = vector("list", k)
	formula_step  = vector("list", k)
	fit_ = vector("list", k)
	names(fit_) = names(latent_var)

	# Deconstructing formula into parts (With latent_var and without latent_var)
	formula_full = stats::terms(formula,simplify=TRUE)
	formula_outcome = get.vars(formula)[1]
	formula_elem_ = attributes(formula_full)$term.labels
	# Adding white spaces before and after to recognize a "E" as opposed to another string like "Elephant"
	formula_elem = paste("", formula_elem_,"")
	for (i in 1:k) index_with_latent_var[[i]] = grepl(paste0(" ",names(latent_var)[i]," "),formula_elem, fixed=TRUE) | grepl(paste0(" ",names(latent_var)[i],":"),formula_elem, fixed=TRUE) | grepl(paste0(":",names(latent_var)[i],":"),formula_elem, fixed=TRUE) | grepl(paste0(":",names(latent_var)[i]," "),formula_elem, fixed=TRUE)
	data_expanded = stats::model.matrix(formula, data=data)
	if (colnames(data_expanded)[1] == "(Intercept)"){
		formula_elem = c("1",formula_elem)
		for (i in 1:k) index_with_latent_var[[i]] = c(FALSE,index_with_latent_var[[i]])
	}

	for (i in 1:k){
		## Formulas for reparametrizations in each steps
		formula_elem_withoutlatent_var = formula_elem[!index_with_latent_var[[i]]]
		formula_elem_withoutlatent_var[-length(formula_elem_withoutlatent_var)] = paste0(formula_elem_withoutlatent_var[-length(formula_elem_withoutlatent_var)], " + ")
		formula_withoutlatent_var[[i]] = paste0(formula_outcome, " ~ ", paste0(formula_elem_withoutlatent_var,collapse=""))
		if (formula_elem[1] != "1") formula_withoutlatent_var[[i]] = paste0(formula_withoutlatent_var[[i]], " - 1")
		formula_withoutlatent_var[[i]] = stats::as.formula(formula_withoutlatent_var[[i]])

		formula_elem_withlatent_var = formula_elem[index_with_latent_var[[i]]]
		# Remove G elements from formula because we want (b1 + b2*E + ...)*G rather than b1*G + b2*E*G + ...
		formula_elem_withlatent_var = gsub(paste0(" ",names(latent_var)[i]," "),"1",formula_elem_withlatent_var, fixed=TRUE)
		formula_elem_withlatent_var = gsub(paste0(" ",names(latent_var)[i],":"),"",formula_elem_withlatent_var, fixed=TRUE)
		formula_elem_withlatent_var = gsub(paste0(":",names(latent_var)[i],":"),":",formula_elem_withlatent_var, fixed=TRUE)
		formula_elem_withlatent_var = gsub(paste0(":",names(latent_var)[i]," "),"",formula_elem_withlatent_var, fixed=TRUE)
		formula_elem_withlatent_var[-length(formula_elem_withlatent_var)] = paste0(formula_elem_withlatent_var[-length(formula_elem_withlatent_var)], " + ")
		formula_withlatent_var[[i]] = paste0(formula_outcome, " ~ ", paste0(formula_elem_withlatent_var,collapse=""))
		if (sum(grepl("1",formula_elem_withlatent_var, fixed=TRUE))==0) formula_withlatent_var[[i]] = paste0(formula_withlatent_var[[i]], " - 1")
		formula_withlatent_var[[i]] = stats::as.formula(formula_withlatent_var[[i]])

		# Making formula for step i
		latent_var_names = colnames(latent_var[[i]])
		latent_var_names[-length(latent_var[[i]])] = paste0(colnames(latent_var[[i]])[-length(latent_var[[i]])], " + ")
		formula_step[[i]] = paste0(formula_outcome, " ~ ", paste0(latent_var_names,collapse=""))
		formula_step[[i]] = paste0(formula_step[[i]], " offset(R0_",i,") - 1")
		formula_step[[i]] = stats::as.formula(formula_step[[i]])
	}

	for (j in 1:maxiter){
		## Step a : fit main model
		fit_a = stats::glm(formula, data=data, family=family, y=FALSE, model=FALSE)
		conv_latent_var = TRUE
		for (i in 1:k){
			if (NCOL(latent_var[[i]])>1){
				# Reparametrizing variables for step i (estimating i-th latent_var)
				data_expanded_withoutlatent_var = stats::model.matrix(formula_withoutlatent_var[[i]], data=data)
				data[,paste0("R0_",i)] = data_expanded_withoutlatent_var%*%stats::coef(fit_a)[!index_with_latent_var[[i]]]
				data_expanded_withlatent_var = stats::model.matrix(formula_withlatent_var[[i]], data=data)
				R1 = data_expanded_withlatent_var%*%stats::coef(fit_a)[index_with_latent_var[[i]]]
				R1_latent_var = latent_var[[i]]*as.vector(R1)
				data[,colnames(latent_var[[i]])]=R1_latent_var

				## Step i-th : fit model for i-th latent_var
				fit_[[i]] = stats::glm(formula_step[[i]], data=data, family=family, y=FALSE, model=FALSE)
				weights_latent_var_ = stats::coef(fit_[[i]])

				# Updating latent_var estimates and checking convergence
				weights_latent_var_old[[i]] = weights_latent_var[[i]]
				weights_latent_var[[i]] = weights_latent_var_/sum(abs(weights_latent_var_))
				data[,names(latent_var)[i]] = latent_var[[i]]%*%weights_latent_var[[i]]
				if(sqrt(sum((weights_latent_var_old[[i]]-weights_latent_var[[i]])^2)) < eps) conv_latent_var = conv_latent_var & TRUE
				else conv_latent_var = FALSE
			}
			else conv_latent_var = conv_latent_var & TRUE
		}
		if (conv_latent_var) break
	}

	# Rerunning last time and scaling to return as results

	fit_a = stats::glm(formula, data=data, family=family, y=FALSE, model=FALSE)

	warn = FALSE
	total_rank = 0
	for (i in 1:k){
		# Reparametrizing variables for step i (estimating i-th latent_var)
		data_expanded_withoutlatent_var = stats::model.matrix(formula_withoutlatent_var[[i]], data=data)
		data[,paste0("R0_",i)] = data_expanded_withoutlatent_var%*%stats::coef(fit_a)[!index_with_latent_var[[i]]]
		data_expanded_withlatent_var = stats::model.matrix(formula_withlatent_var[[i]], data=data)
		R1 = data_expanded_withlatent_var%*%stats::coef(fit_a)[index_with_latent_var[[i]]]
		R1_latent_var = latent_var[[i]]*as.vector(R1)
		data[,colnames(latent_var[[i]])]=R1_latent_var

		fit_[[i]] = stats::glm(formula_step[[i]], data=data, family=family, y=FALSE, model=FALSE)
		# Changing data temporarily to get model to make sure the model parameters sum to 1
		temp = data[,colnames(latent_var[[i]])]
		data[,colnames(latent_var[[i]])] = data[,colnames(latent_var[[i]])]*sum(abs(stats::coef(fit_[[i]])))
		fit_[[i]] = stats::glm(formula_step[[i]], data=data, family=family, y=FALSE, model=FALSE)
		data[,colnames(latent_var[[i]])] = temp
		data[,names(latent_var)[i]]	= latent_var[[i]] %*% stats::coef(fit_[[i]])
		if (abs(((fit_a$deviance-fit_[[i]]$deviance)/fit_a$deviance))>=.01 && !warn){
			warning("Deviance differs by more than 1% between model parts. Make sure that everything was set up properly and try increasing the number of iterations (maxiter).")
			warn = TRUE
		}
		total_rank = total_rank + fit_[[i]]$rank - 1
	}

	# Make sure data make sense for user (and for plot function)
	fit_a$data = data
	for (i in 1:k){
		fit_[[i]]$data = data
	}

	#Change some arguments so that we get the right AIC, BIC and dispersion for the model
	true_rank = fit_a$rank + total_rank
	# true_rank might be different from df of logLik, this is because glm() assume that variance components should be counted
	# I don't agree but we need to follow the same guideline to make sure it matches with glm()
	logLik_df =  attr(logLik(fit_a), "df") + total_rank
	true_aic = 2*logLik_df - 2*logLik(fit_a)[1]
	true_aicc = true_aic + ((2*logLik_df*(logLik_df+1))/(stats::nobs(fit_a)-logLik_df-1))
	true_bic = log(stats::nobs(fit_a))*logLik_df - 2*logLik(fit_a)[1]
	true_df.residual = stats::nobs(fit_a) - true_rank
	true_null.deviance = fit_a$null.deviance

	# print convergences stuff;
	if (conv_latent_var){
		if (print) cat(paste0("Converged in ",j, " iterations\n"))
	} 
	else{
		warning(paste0("Did not reach convergence in maxiter iterations. Try increasing maxiter or make eps smaller."))
	}
	if (is.null(ylim)) result = list(fit_main = fit_a, fit_latent_var = fit_, true_model_parameters=list(AIC = true_aic, AICc = true_aicc, BIC = true_bic, rank = true_rank, df.residual = true_df.residual, null.deviance=true_null.deviance))
	else result = list(fit_main = fit_a, fit_latent_var = fit_, true_model_parameters=list(AIC = true_aic, AICc = true_aicc, BIC = true_bic, rank = true_rank, df.residual = true_df.residual, null.deviance=true_null.deviance), ylim=ylim)
	class(result) <- "IMLEGIT"
	return(result)
}

elastic_net_var_select = function(data, latent_var, formula, latent_var_searched=NULL, cross_validation=FALSE, alpha=.75, standardize=TRUE, lambda_path=NULL, lambda_mult=1, lambda_min = .0001, n_lambda = 100, start_latent_var=NULL, eps=.001, maxiter=100, family=gaussian, ylim=NULL, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, classification=FALSE, print=TRUE)
{
	if (!is.null(ylim)){
		if (!is.numeric(ylim) || length(ylim) !=2) stop("ylim must either be NULL or a numeric vector of size two")
	}
	# Setting up latent_var and checks
	if (class(latent_var)!="list") stop("latent_var must be a list of datasets")
	k = length(latent_var)
	if (k==0) stop("latent_var cannot be an empty list")
	if (is.null(names(latent_var))){
		if (print) cat("You have not specified names for the latent variables, assigning names to latent_var is highly recommended to prevent confusion. For now, they will be named L1, L2, ...\n")
		names(latent_var) = paste0("L",1:k)
	}
	for (i in 1:k){
		latent_var[[i]] = as.matrix(data.frame(latent_var[[i]],fix.empty.names=FALSE))
		if (sum(colnames(latent_var[[i]])=="") > 0){
			if (print) cat(paste0("You have not specified column names for certain elements in ",names(latent_var)[i], ", elements of this latent variable will be named ",names(latent_var)[i],1,", ",names(latent_var)[i],2," ...\n"))
			colnames(latent_var[[i]]) = paste0(names(latent_var)[i],1:NCOL(latent_var[[i]]))
		}
	}
	# More checks
	if (maxiter <= 0) warning("maxiter must be > 0")
	if (k > 1) for (i in 1:(k-1)) if (NROW(latent_var[[i]]) != NROW(latent_var[[i+1]])) stop("Some datasets in latent_var don't have the same number of observations")
	if(!is.null(start_latent_var)){
		if (class(start_latent_var)!="list") stop("start_latent_var must be a lit of vectors (or NULL)")
		if (k!=length(start_latent_var)) stop("start_latent_var must have the same size as latent_var")
		for (i in 1:k){
			if (!is.null(latent_var[[i]])){
				if (NCOL(latent_var[[i]])!=length(start_latent_var[[i]])) stop("All elements of start_latent_var must either be NULL or have the same length as the number of the elements in its associated latent variable")
			}
		}
	}
	if (class(data) != "data.frame" && class(data) != "matrix") stop("data must be a data.frame")

	# getting right formats
	# Retaining only the needed variables from the dataset (need to set elements in latent_var for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	for (i in 1:k) data[,names(latent_var)[i]] = 0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	formula = stats::as.formula(formula)

	# Error message about factors
	if (sum(apply(data,2,is.numeric)) != NCOL(data)) stop("All variables used must be numeric, factors are not allowed. Please dummy code all categorical variables inside your datasets (data, latent_var[[1]], latent_var[[2]], ...)")
	for (i in 1:k) if (sum(apply(latent_var[[i]],2,is.numeric)) != NCOL(latent_var[[i]])) stop("All variables used must be numeric, factors are not allowed. Please dummy code all categorical variables inside your datasets (data, latent_var[[1]], latent_var[[2]], ...)")

	# remove missing data
	comp = stats::complete.cases(data,latent_var[[1]])
	if (k > 1) for (i in 2:k) comp = comp & stats::complete.cases(latent_var[[i]])
	data = data[comp,, drop=FALSE]
	for (i in 1:k) latent_var[[i]] = latent_var[[i]][comp,, drop=FALSE]
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	#Adding empty variables in main dataset for latent_var
	for (i in 1:k){
		data[,colnames(latent_var[[i]])]=0
		data[,paste0("R0_",i)]=0
	}

	formula_outcome = get.vars(formula)[1]
	if (is.null(latent_var_searched)) latent_var_searched = c(1:k)

	## Standardize variables: (need to use n instead of (n-1) as denominator)
	mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))

	#### Lambda path (taken from https://stackoverflow.com/questions/23686067/default-lambda-sequence-in-glmnet-for-cross-validation)
	if (is.null(lambda_path)){

		# We find lambda max for all latent_var and use the maximum
		n = dim(data)[1]
		lambda_max = -Inf
		for (i in latent_var_searched){
			sx = scale(latent_var[[i]], scale = apply(latent_var[[i]], 2, mysd))
			sx = as.matrix(sx, ncol = NCOL(latent_var[[i]]), nrow = n)
			sy = as.vector(scale(data[,formula_outcome], scale = mysd(data[,formula_outcome])))
			lambda_max_current = max(lambda_mult*abs(colSums(sx*sy)))/n # Multiplying by k because otherwise its too small
			if (lambda_max < lambda_max_current) lambda_max = lambda_max_current
		}

		## Calculate lambda path (first get lambda_max):
		lambda_path = round(exp(seq(log(lambda_max), log(lambda_max*lambda_min), length.out = n_lambda)), digits = 10)
	}

	if (standardize){
		for (i in 1:k){
			sx = scale(latent_var[[i]], scale = apply(latent_var[[i]], 2, mysd))
			latent_var[[i]] = as.matrix(sx, ncol = NCOL(latent_var[[i]]), nrow = dim(data)[1])
		}
	}

	if (cross_validation){
		if (classification) results = matrix(NA,length(lambda_path),7)
		else results = matrix(NA,length(lambda_path),6)
	}
	else results = matrix(NA,length(lambda_path),3)
	if (cross_validation & classification) colnames(results) = c("AIC","AICc","BIC","cv_R2","cv_Huber","cv_L1","cv_AUC")
	else if (cross_validation & !classification) colnames(results) = c("AIC","AICc","BIC","cv_R2","cv_Huber","cv_L1")
	else colnames(results) = c("AIC","AICc","BIC")
	fit = vector("list", length(lambda_path))
	glmnet_coef = vector("list", length(lambda_path))
	n_var_zero_prev = -Inf
	for (i in 1:length(lambda_path)){
		result = IMLEGIT_net(data=data, latent_var=latent_var, formula=formula, latent_var_searched=latent_var_searched, cross_validation = FALSE, alpha=alpha, lambda=lambda_path[i], start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, family_string=as.character(substitute(family)), ylim=ylim, print=FALSE, warn=FALSE)

		# Only do cross-validation, if worth it (i.e., if a variable has been now been added to the list of the ones used). This reduces computation.
		n_var_zero = sum(ceiling(abs(result$glmnet_coef))==0)
		if (cross_validation & n_var_zero_prev != n_var_zero & !is.null(result$fit)) result = IMLEGIT_net(data=data, latent_var=latent_var, formula=formula, latent_var_searched=latent_var_searched, cross_validation = TRUE, alpha=alpha, lambda=lambda_path[i], start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, family_string=as.character(substitute(family)), ylim=ylim, print=FALSE, warn=FALSE)
		n_var_zero_prev = n_var_zero

		if (!is.null(result$fit)) fit[[i]] = result$fit
		glmnet_coef[[i]] = result$glmnet_coef
		if (is.null(fit[[i]])){
			if (cross_validation & classification) results[i,] = rep(NA, 7)
			else if (cross_validation & !classification) results[i,] = rep(NA, 6)
			else results[i,] = rep(NA, 3)
		}
		else{
			if (cross_validation){
				if (is.null(result$fit_cv)) R2_cv = results[max(i-1,1),"cv_R2"]
				else R2_cv = mean(result$fit_cv$R2_cv)
				if (is.null(result$fit_cv)) R2_h = results[max(i-1,1),"cv_Huber"]
				else R2_h = mean(result$fit_cv$Huber_cv)
				if (is.null(result$fit_cv)) R2_l1 = results[max(i-1,1),"cv_L1"]
				else R2_l1 = mean(result$fit_cv$L1_cv)
				if (classification){
					if (is.null(result$fit_cv)) AUC = results[max(i-1,1),"cv_AUC"]
					else AUC = mean(result$fit_cv$AUC)
				}
			}
			if (cross_validation & classification) results[i,] = c(fit[[i]]$true_model_parameters$AIC, fit[[i]]$true_model_parameters$AICc, fit[[i]]$true_model_parameters$BIC, R2_cv, R2_h, R2_l1, AUC)
			else if (cross_validation & !classification) results[i,] = c(fit[[i]]$true_model_parameters$AIC, fit[[i]]$true_model_parameters$AICc, fit[[i]]$true_model_parameters$BIC, R2_cv,R2_h, R2_l1)
			else results[i,] = c(fit[[i]]$true_model_parameters$AIC, fit[[i]]$true_model_parameters$AICc, fit[[i]]$true_model_parameters$BIC)
		}
	}
	out = list(results=results, fit=fit, glmnet_coef=glmnet_coef, lambda_path=lambda_path)
	class(out) = "elastic_net_var_select"
	return(out)
}

IMLEGIT_net = function(data, latent_var, formula, latent_var_searched=NULL, cross_validation=FALSE, alpha=1, lambda=.0001, start_latent_var=NULL, eps=.001, maxiter=100, family=gaussian, ylim=NULL, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, classification=FALSE, print=TRUE, warn=TRUE, family_string=NULL)
{
	if (!is.null(ylim)){
		if (!is.numeric(ylim) || length(ylim) !=2) stop("ylim must either be NULL or a numeric vector of size two")
	}
	# Setting up latent_var and checks
	if (class(latent_var)!="list") stop("latent_var must be a list of datasets")
	k = length(latent_var)
	if (k==0) stop("latent_var cannot be an empty list")
	if (is.null(names(latent_var))){
		if (print) cat("You have not specified names for the latent variables, assigning names to latent_var is highly recommended to prevent confusion. For now, they will be named L1, L2, ...\n")
		names(latent_var) = paste0("L",1:k)
	}
	for (i in 1:k){
		latent_var[[i]] = as.matrix(data.frame(latent_var[[i]],fix.empty.names=FALSE))
		if (sum(colnames(latent_var[[i]])=="") > 0){
			if (print) cat(paste0("You have not specified column names for certain elements in ",names(latent_var)[i], ", elements of this latent variable will be named ",names(latent_var)[i],1,", ",names(latent_var)[i],2," ...\n"))
			colnames(latent_var[[i]]) = paste0(names(latent_var)[i],1:NCOL(latent_var[[i]]))
		}
	}

	# More checks
	if (maxiter <= 0) warning("maxiter must be > 0")
	if (k > 1) for (i in 1:(k-1)) if (NROW(latent_var[[i]]) != NROW(latent_var[[i+1]])) stop("Some datasets in latent_var don't have the same number of observations")
	if(!is.null(start_latent_var)){
		if (class(start_latent_var)!="list") stop("start_latent_var must be a lit of vectors (or NULL)")
		if (k!=length(start_latent_var)) stop("start_latent_var must have the same size as latent_var")
		for (i in 1:k){
			if (!is.null(latent_var[[i]])){
				if (NCOL(latent_var[[i]])!=length(start_latent_var[[i]])) stop("All elements of start_latent_var must either be NULL or have the same length as the number of the elements in its associated latent variable")
			}
		}
	}
	if (class(data) != "data.frame" && class(data) != "matrix") stop("data must be a data.frame")

	# getting right formats
	# Retaining only the needed variables from the dataset (need to set elements in latent_var for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	for (i in 1:k) data[,names(latent_var)[i]] = 0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	formula = stats::as.formula(formula)

	# Error message about factors
	if (sum(apply(data,2,is.numeric)) != NCOL(data)) stop("All variables used must be numeric, factors are not allowed. Please dummy code all categorical variables inside your datasets (data, latent_var[[1]], latent_var[[2]], ...)")
	for (i in 1:k) if (sum(apply(latent_var[[i]],2,is.numeric)) != NCOL(latent_var[[i]])) stop("All variables used must be numeric, factors are not allowed. Please dummy code all categorical variables inside your datasets (data, latent_var[[1]], latent_var[[2]], ...)")

	# remove missing data
	comp = stats::complete.cases(data,latent_var[[1]])
	if (k > 1) for (i in 2:k) comp = comp & stats::complete.cases(latent_var[[i]])
	data = data[comp,, drop=FALSE]
	for (i in 1:k) latent_var[[i]] = latent_var[[i]][comp,, drop=FALSE]
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	#Adding empty variables in main dataset for latent_var
	for (i in 1:k){
		data[,colnames(latent_var[[i]])]=0
		data[,paste0("R0_",i)]=0
	}

	# Setting up initial weighted latent_var
	weights_latent_var_old = vector("list", k)
	weights_latent_var = vector("list", k)
	if (is.null(start_latent_var)){
		for (i in 1:k) weights_latent_var[[i]] = rep(1/dim(latent_var[[i]])[2],dim(latent_var[[i]])[2])
	}
	else{
		for (i in 1:k){
			if (sum(abs(start_latent_var[[i]]))==0) weights_latent_var[[i]] = rep(1/dim(latent_var[[i]])[2],dim(latent_var[[i]])[2])
			else weights_latent_var[[i]] = start_latent_var[[i]]/sum(abs(start_latent_var[[i]]))
		}		
	}
	for (i in 1:k) data[,names(latent_var)[i]] = latent_var[[i]]%*%weights_latent_var[[i]]

	# Lists needed for later
	index_with_latent_var = vector("list", k)
	formula_withoutlatent_var  = vector("list", k)
	formula_withlatent_var  = vector("list", k)
	formula_step  = vector("list", k)
	fit_ = vector("list", k)
	names(fit_) = names(latent_var)

	# Deconstructing formula into parts (With latent_var and without latent_var)
	formula_full = stats::terms(formula,simplify=TRUE)
	formula_outcome = get.vars(formula)[1]
	formula_elem_ = attributes(formula_full)$term.labels
	# Adding white spaces before and after to recognize a "E" as opposed to another string like "Elephant"
	formula_elem = paste("", formula_elem_,"")
	for (i in 1:k) index_with_latent_var[[i]] = grepl(paste0(" ",names(latent_var)[i]," "),formula_elem, fixed=TRUE) | grepl(paste0(" ",names(latent_var)[i],":"),formula_elem, fixed=TRUE) | grepl(paste0(":",names(latent_var)[i],":"),formula_elem, fixed=TRUE) | grepl(paste0(":",names(latent_var)[i]," "),formula_elem, fixed=TRUE)
	data_expanded = stats::model.matrix(formula, data=data)
	if (colnames(data_expanded)[1] == "(Intercept)"){
		intercept=TRUE
		formula_elem = c("1",formula_elem)
		for (i in 1:k) index_with_latent_var[[i]] = c(FALSE,index_with_latent_var[[i]])
	}
	else intercept=FALSE

	for (i in 1:k){
		## Formulas for reparametrizations in each steps
		formula_elem_withoutlatent_var = formula_elem[!index_with_latent_var[[i]]]
		formula_elem_withoutlatent_var[-length(formula_elem_withoutlatent_var)] = paste0(formula_elem_withoutlatent_var[-length(formula_elem_withoutlatent_var)], " + ")
		formula_withoutlatent_var[[i]] = paste0(formula_outcome, " ~ ", paste0(formula_elem_withoutlatent_var,collapse=""))
		if (formula_elem[1] != "1") formula_withoutlatent_var[[i]] = paste0(formula_withoutlatent_var[[i]], " - 1")
		formula_withoutlatent_var[[i]] = stats::as.formula(formula_withoutlatent_var[[i]])

		formula_elem_withlatent_var = formula_elem[index_with_latent_var[[i]]]
		# Remove G elements from formula because we want (b1 + b2*E + ...)*G rather than b1*G + b2*E*G + ...
		formula_elem_withlatent_var = gsub(paste0(" ",names(latent_var)[i]," "),"1",formula_elem_withlatent_var, fixed=TRUE)
		formula_elem_withlatent_var = gsub(paste0(" ",names(latent_var)[i],":"),"",formula_elem_withlatent_var, fixed=TRUE)
		formula_elem_withlatent_var = gsub(paste0(":",names(latent_var)[i],":"),":",formula_elem_withlatent_var, fixed=TRUE)
		formula_elem_withlatent_var = gsub(paste0(":",names(latent_var)[i]," "),"",formula_elem_withlatent_var, fixed=TRUE)
		formula_elem_withlatent_var[-length(formula_elem_withlatent_var)] = paste0(formula_elem_withlatent_var[-length(formula_elem_withlatent_var)], " + ")
		formula_withlatent_var[[i]] = paste0(formula_outcome, " ~ ", paste0(formula_elem_withlatent_var,collapse=""))
		if (sum(grepl("1",formula_elem_withlatent_var, fixed=TRUE))==0) formula_withlatent_var[[i]] = paste0(formula_withlatent_var[[i]], " - 1")
		formula_withlatent_var[[i]] = stats::as.formula(formula_withlatent_var[[i]])

		# Making formula for step i
		latent_var_names = colnames(latent_var[[i]])
		latent_var_names[-length(latent_var[[i]])] = paste0(colnames(latent_var[[i]])[-length(latent_var[[i]])], " + ")
		formula_step[[i]] = paste0(formula_outcome, " ~ ", paste0(latent_var_names,collapse=""))
		formula_step[[i]] = paste0(formula_step[[i]], " offset(R0_",i,") - 1")
		formula_step[[i]] = stats::as.formula(formula_step[[i]])
	}

	if (is.null(latent_var_searched)) latent_var_searched = c(1:k)

	done = FALSE
	for (j in 1:maxiter){
		## Step a : fit main model
		fit_a = stats::glm(formula, data=data, family=family, y=FALSE, model=FALSE)
		# L1 standardize the parameters of the main model
		weights_a = stats::coef(fit_a)
		weights_a = weights_a/sum(abs(weights_a))

		conv_latent_var = TRUE
		for (i in 1:k){
			if (NCOL(latent_var[[i]])>1){
				# Reparametrizing variables for step i (estimating i-th latent_var)
				data_expanded_withoutlatent_var = stats::model.matrix(formula_withoutlatent_var[[i]], data=data)
				data[,paste0("R0_",i)] = data_expanded_withoutlatent_var%*%weights_a[!index_with_latent_var[[i]]]
				data_expanded_withlatent_var = stats::model.matrix(formula_withlatent_var[[i]], data=data)
				R1 = data_expanded_withlatent_var%*%weights_a[index_with_latent_var[[i]]]
				R1_latent_var = latent_var[[i]]*as.vector(R1)
				data[,colnames(latent_var[[i]])]=R1_latent_var

				## Step i-th : fit model for i-th latent_var
				if (mean(is.na(R1))){
					if (warn) warning(paste0("Lambda is too big, one latent variable has weight 0 everywhere. Use a smaller lambda."))
					weights_latent_var[[i]] = rep(0, NCOL(latent_var[[i]]))
					return(list(fit=NULL, glmnet_coef=unlist(weights_latent_var)))
				}
				else{

					if (i %in% latent_var_searched){ #glmnet
						if (is.null(family_string)) fit_[[i]] = glmnet::glmnet(as.matrix(latent_var[[i]]), data[,formula_outcome], family=as.character(substitute(family)), alpha=alpha, lambda=lambda, offset=data$R0_1, intercept=FALSE, standardize=FALSE)
						else fit_[[i]] = glmnet::glmnet(as.matrix(latent_var[[i]]), data[,formula_outcome], family=family_string, alpha=alpha, lambda=lambda, offset=data$R0_1, intercept=FALSE)
						weights_latent_var_ = as.numeric(stats::coef(fit_[[i]])[-1]) # removing the 0 intercept
					}
					else{ #glm
						fit_[[i]] = stats::glm(formula_step[[i]], data=data, family=family, y=FALSE, model=FALSE)
						weights_latent_var_ = stats::coef(fit_[[i]])
					}

					# Updating latent_var estimates and checking convergence
					weights_latent_var_old[[i]] = weights_latent_var[[i]]
					if (i %in% latent_var_searched) weights_latent_var[[i]] = weights_latent_var_
					else weights_latent_var[[i]] = weights_latent_var_/sum(abs(weights_latent_var_))
					data[,names(latent_var)[i]] = latent_var[[i]]%*%weights_latent_var[[i]]
					if(sqrt(sum((weights_latent_var_old[[i]]-weights_latent_var[[i]])^2)) < eps) conv_latent_var = conv_latent_var & TRUE
					else conv_latent_var = FALSE
					#print(weights_latent_var_)
				}
			}
			else conv_latent_var = conv_latent_var & TRUE
		}
		if (conv_latent_var) break
	}

	# Rerunning last time and scaling to return as results
	# We simply fit a IMLEGIT with the same starting points (they are automatically standardized by IMLEGIT).
	removed = rep(FALSE,k)
	weights_latent_var_backup = weights_latent_var
	for (i in 1:k){
		names(weights_latent_var_backup[[i]]) = colnames(latent_var[[i]]) # backup names
		latent_var[[i]] = latent_var[[i]][,colnames(latent_var[[i]])[weights_latent_var[[i]]!=0],drop=FALSE] # Returns only the elements of each latent var which has weight non-zero
		weights_latent_var[[i]] = weights_latent_var[[i]][weights_latent_var[[i]]!=0]
	}
	result = IMLEGIT(data=data, latent_var=latent_var, formula=formula, start_latent_var=weights_latent_var, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=print)
	if (cross_validation) result_cv = IMLEGIT_cv(data=data, latent_var=latent_var, formula=formula, start_latent_var=weights_latent_var, eps=eps, maxiter=maxiter, family=family, ylim=ylim, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification)
	else result_cv = NULL
	# print convergences stuff;
	if (conv_latent_var){
		if (print) cat(paste0("Converged in ",j, " iterations\n"))
	} 
	else{
		warning(paste0("Did not reach convergence in maxiter iterations. Try increasing maxiter or make eps smaller."))
	}
	return(list(fit=result, glmnet_coef=unlist(weights_latent_var_backup), fit_cv=result_cv))
}

summary.elastic_net_var_select = function(object, ...){
	# Only grab the index of the unique IMLEGIT models (we only care when a variable is dropped out)
	n_var_zero_prev = -Inf
	indexes_zero = c()

	for (i in 1:length(object$lambda_path)){
		n_var_zero = sum(ceiling(abs(object$glmnet_coef[[i]]))==0)
		if (n_var_zero_prev != n_var_zero){
			indexes_zero = c(indexes_zero, i)
			if (i == 1) coef = rbind(pmin(ceiling(abs(object$glmnet_coef[[i]])),1))
			else coef = rbind(coef,pmin(ceiling(abs(object$glmnet_coef[[i]])),1))
		}
		n_var_zero_prev = n_var_zero
	}
	results = cbind(object$lambda_path[indexes_zero], indexes_zero, object$results[indexes_zero,,drop=FALSE],coef)
	colnames(results)[1] = "Lambda"
	colnames(results)[2] = "Model index"
	rownames(results) = NULL
	return(results)
}

best_model <- function(object, ...) UseMethod("best_model")

best_model.elastic_net_var_select = function(object, criterion, ...){
	if (criterion == "cv_R2" | criterion=="cv_AUC") i = which.max(object$results[,criterion])
	else i = which.min(object$results[,criterion])
	return(list(results=object$results[i,], fit=object$fit[[i]], coef=pmin(ceiling(abs(object$glmnet_coef[[i]])),1), lambda=object$lambda_path[i], index = i))
}

plot.elastic_net_var_select = function(x, lwd=2, start=1, ...){
	object = x
	n = length(object$glmnet_coef)
	k = length(object$fit[[length(object$glmnet_coef)]]$fit_latent_var)
	lty_ = c()
	n_var = vector("list", k)
	n_var_total = 0
	for (i in 1:k){
		if (is.null(object$fit[[length(object$glmnet_coef)]])) stop("No fit to plot. This must be because all your choices of lambda were extremely large.")
		n_var[[i]] = length(coef(object$fit[[length(object$glmnet_coef)]]$fit_latent_var[[i]]))
		lty_ = c(lty_, rep(i,n_var[[i]]))
		n_var_total = n_var_total + n_var[[i]]
	}
	A = matrix(unlist(object$glmnet_coef), ncol = n_var_total, byrow = TRUE)
	path = object$lambda_path[start:n]
	A = A[start:n,]

	matplot(path,A,type="l", col=RColorBrewer::brewer.pal(n_var_total,"Paired"), xlab="log(lambda)", ylab="Coefficients", lty=lty_, lwd=lwd)
	abline(0,0, col="grey", lty=3)
	legend("topright",legend=names(object$glmnet_coef[[n]]), col=RColorBrewer::brewer.pal(n_var_total,"Paired"), lty=lty_, bg="white", lwd=lwd)
}

predict.LEGIT = function(object, data, genes, env, ...){
	data = data.frame(data)
	genes = as.matrix(genes, drop=FALSE)
	env = as.matrix(env, drop=FALSE)
	# missing data
	comp = stats::complete.cases(genes,env)
	data = data[comp,, drop=FALSE]
	genes = genes[comp,, drop=FALSE]
	data$G = genes%*%stats::coef(object[[2]])
	if (!is.null(object$crossover) && !object$crossover_fixed) data$E = cbind(-1,env)%*%stats::coef(object[[3]])
	else if (!is.null(object$crossover) && object$crossover_fixed) data$E = env%*%stats::coef(object[[3]]) - object$crossover
	else data$E = env%*%stats::coef(object[[3]])
	results = stats::predict.glm(object[[1]], newdata=data, ...)
	if (is.null(object$ylim)) return(results)
	else return(pmin(pmax(results,object$ylim[1]),object$ylim[2]))
}

predict.IMLEGIT = function(object, data, latent_var, ...){
	data = data.frame(data)
	k = length(latent_var)
	for (i in 1:k) latent_var[[i]] = as.matrix(data.frame(latent_var[[i]],fix.empty.names=FALSE))
	if (is.null(names(latent_var))){
		cat("You have not specified names for the latent variables, assigning names to latent_var is highly recommended to prevent confusion. For now, they will be named L1, L2, ...\n")
		names(latent_var) = paste0("L",1:k)
	}
	# remove missing data
	comp = stats::complete.cases(latent_var[[1]])
	if (k > 1) for (i in 2:k) comp = comp & stats::complete.cases(latent_var[[i]])
	for (i in 1:k) latent_var[[i]] = latent_var[[i]][comp,, drop=FALSE]
	for (i in 1:k) data[,names(latent_var)[i]] = latent_var[[i]]%*%stats::coef(object[[2]][[i]])
	results = stats::predict.glm(object[[1]], newdata=data, ...)
	if (is.null(object$ylim)) return(results)
	else return(pmin(pmax(results,object$ylim[1]),object$ylim[2]))
}

summary.LEGIT = function(object, ...){
	lapply(object[1:3],function(object_current, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...){
		# Using the right values
		object_current$aic = object$true_model_parameters$AIC
		object_current$rank = object$true_model_parameters$rank
		object_current$df.residual = object$true_model_parameters$df.residual
		object_current$null.deviance = object$true_model_parameters$null.deviance
	    est.disp <- FALSE
	    df.r <- object_current$df.residual
	    if (is.null(dispersion)) 
	        dispersion <- if (object_current$family$family %in% c("poisson", "binomial")) 1
	        else if (df.r > 0) {
	            est.disp <- TRUE
	            if (any(object_current$weights == 0)) 
	                warning("observations with zero weight not used for calculating dispersion")
	            sum((object_current$weights * object_current$residuals^2)[object_current$weights > 
	                0])/df.r
	        }
	        else {
	            est.disp <- TRUE
	            NaN
	        }
	    aliased <- is.na(stats::coef(object_current))
	    p <- object_current$qr$rank
	    if (p > 0) {
	        p1 <- 1L:p
	        coef.p <- object_current$coefficients[object_current$qr$pivot[p1]]
	        covmat.unscaled <- chol2inv(object_current$qr$qr)
	        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
	        covmat <- dispersion * covmat.unscaled
	        var.cf <- diag(covmat)
	        s.err <- sqrt(var.cf)
	        tvalue <- coef.p/s.err
	        dn <- c("Estimate", "Std. Error")
	        if (!est.disp) {
	            pvalue <- 2 * pnorm(-abs(tvalue))
	            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
	            dimnames(coef.table) <- list(names(coef.p), c(dn, 
	                "z value", "Pr(>|z|)"))
	        }
	        else if (df.r > 0) {
	            pvalue <- 2 * stats::pt(-abs(tvalue), df.r)
	            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
	            dimnames(coef.table) <- list(names(coef.p), c(dn, 
	                "t value", "Pr(>|t|)"))
	        }
	        else {
	            coef.table <- cbind(coef.p, NaN, NaN, NaN)
	            dimnames(coef.table) <- list(names(coef.p), c(dn, 
	                "t value", "Pr(>|t|)"))
	        }
	        df.f <- NCOL(object_current$qr$qr)
	    }
	    else {
	        coef.table <- matrix(, 0L, 4L)
	        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
	            "t value", "Pr(>|t|)"))
	        covmat.unscaled <- covmat <- matrix(, 0L, 0L)
	        df.f <- length(aliased)
	    }
	    keep <- match(c("call", "terms", "family", "deviance", "aic", 
	        "contrasts", "df.residual", "null.deviance", "df.null", 
	        "iter", "na.action"), names(object_current), 0L)
	    ans <- c(object_current[keep], list(deviance.resid = stats::residuals(object_current, 
	        type = "deviance"), coefficients = coef.table, aliased = aliased, 
	        dispersion = dispersion, df = c(object_current$rank, df.r, df.f), 
	        cov.unscaled = covmat.unscaled, cov.scaled = covmat))
	    if (correlation && p > 0) {
	        dd <- sqrt(diag(covmat.unscaled))
	        ans$correlation <- covmat.unscaled/outer(dd, dd)
	        ans$symbolic.cor <- symbolic.cor
	    }
	    class(ans) <- "summary.glm"
	    return(ans)
	})
}

summary.IMLEGIT = function(object, ...){
	newobject = list(fit_main=object$fit_main)
	for (i in 1:length(object$fit_latent_var)) newobject[[i+1]] = object$fit_latent_var[[i]]
	names(newobject) = c("fit_main",paste0("fit_",names(object$fit_latent_var)))
	lapply(newobject,function(object_current, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...){
		# Using the right values
		object_current$aic = object$true_model_parameters$AIC
		object_current$rank = object$true_model_parameters$rank
		object_current$df.residual = object$true_model_parameters$df.residual
		object_current$null.deviance = object$true_model_parameters$null.deviance
	    est.disp <- FALSE
	    df.r <- object_current$df.residual
	    if (is.null(dispersion)) 
	        dispersion <- if (object_current$family$family %in% c("poisson", "binomial")) 1
	        else if (df.r > 0) {
	            est.disp <- TRUE
	            if (any(object_current$weights == 0)) 
	                warning("observations with zero weight not used for calculating dispersion")
	            sum((object_current$weights * object_current$residuals^2)[object_current$weights > 
	                0])/df.r
	        }
	        else {
	            est.disp <- TRUE
	            NaN
	        }
	    aliased <- is.na(stats::coef(object_current))
	    p <- object_current$qr$rank
	    if (p > 0) {
	        p1 <- 1L:p
	        coef.p <- object_current$coefficients[object_current$qr$pivot[p1]]
	        covmat.unscaled <- chol2inv(object_current$qr$qr)
	        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
	        covmat <- dispersion * covmat.unscaled
	        var.cf <- diag(covmat)
	        s.err <- sqrt(var.cf)
	        tvalue <- coef.p/s.err
	        dn <- c("Estimate", "Std. Error")
	        if (!est.disp) {
	            pvalue <- 2 * pnorm(-abs(tvalue))
	            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
	            dimnames(coef.table) <- list(names(coef.p), c(dn, 
	                "z value", "Pr(>|z|)"))
	        }
	        else if (df.r > 0) {
	            pvalue <- 2 * stats::pt(-abs(tvalue), df.r)
	            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
	            dimnames(coef.table) <- list(names(coef.p), c(dn, 
	                "t value", "Pr(>|t|)"))
	        }
	        else {
	            coef.table <- cbind(coef.p, NaN, NaN, NaN)
	            dimnames(coef.table) <- list(names(coef.p), c(dn, 
	                "t value", "Pr(>|t|)"))
	        }
	        df.f <- NCOL(object_current$qr$qr)
	    }
	    else {
	        coef.table <- matrix(, 0L, 4L)
	        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
	            "t value", "Pr(>|t|)"))
	        covmat.unscaled <- covmat <- matrix(, 0L, 0L)
	        df.f <- length(aliased)
	    }
	    keep <- match(c("call", "terms", "family", "deviance", "aic", 
	        "contrasts", "df.residual", "null.deviance", "df.null", 
	        "iter", "na.action"), names(object_current), 0L)
	    ans <- c(object_current[keep], list(deviance.resid = stats::residuals(object_current, 
	        type = "deviance"), coefficients = coef.table, aliased = aliased, 
	        dispersion = dispersion, df = c(object_current$rank, df.r, df.f), 
	        cov.unscaled = covmat.unscaled, cov.scaled = covmat))
	    if (correlation && p > 0) {
	        dd <- sqrt(diag(covmat.unscaled))
	        ans$correlation <- covmat.unscaled/outer(dd, dd)
	        ans$symbolic.cor <- symbolic.cor
	    }
	    class(ans) <- "summary.glm"
	    return(ans)
	})
}

LEGIT_cv = function (data, genes, env, formula, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, classification=FALSE, start_genes=NULL, start_env=NULL, eps=.001, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, id=NULL, crossover = NULL, crossover_fixed = FALSE){
	if (!is.null(ylim)){
		if (!is.numeric(ylim) || length(ylim) !=2) stop("ylim must either be NULL or a numeric vector of size two")
	}

	# Renaming it because there is already an id variable
	if (!is.null(id)) obs_id = id
	else obs_id = NULL
	id = NULL

	# getting right formats
	# Retaining only the needed variables from the dataset (need to set G and E variables for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	data$G=0
	data$E=0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	genes = as.matrix(genes, drop=FALSE)
	if (is.null(colnames(genes))){
		if (print) cat("You have not specified column names for genes, they will be named gene1, gene2, ...\n")
		colnames(genes) = paste0("gene",1:NCOL(genes))
	}
	env = as.matrix(env, drop=FALSE)
	if (is.null(colnames(env))){
		if (print) cat("You have not specified column names for env, they will be named env1, env2, ...\n")
		colnames(env) = paste0("env",1:NCOL(env))
	}
	formula = stats::as.formula(formula)

	# Error message about factors
	if (sum(apply(data,2,is.numeric)) != NCOL(data) || sum(apply(genes,2,is.numeric)) != NCOL(genes) || sum(apply(env,2,is.numeric)) != NCOL(env)) stop("All variables used must be numeric, factors are not allowed. Please dummy code all categorical variables inside your datasets (data, gene, env)")

	# remove missing data
	comp = stats::complete.cases(data,genes,env)
	data = data[comp,, drop=FALSE]
	genes = genes[comp,, drop=FALSE]
	env = env[comp,, drop=FALSE]
	if (!is.null(obs_id)){
		if (!is.null(dim(obs_id))) obs_id = obs_id[comp,,drop=FALSE]
		else obs_id = obs_id[comp]
	}
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	formula_outcome = get.vars(formula)[1]
	R2_cv = c()
	Huber_cv = c()
	L1_cv = c()
	AUC = c()
	best_threshold = c()
	roc_curve = list()
	residuals = rep(0,dim(data)[1])
	pearson_residuals = rep(0,dim(data)[1])

	if (!is.null(folds)) cv_iter = length(folds)	

	for (j in 1:cv_iter){
		if (!is.null(seed)) set.seed(seed*j)
		# Folds
		if (is.null(folds)){
			s = sample(NROW(data))
			data_n = data[s,, drop=FALSE]
			genes_n = genes[s,, drop=FALSE]
			env_n = env[s,, drop=FALSE]
			id = cut(seq(1,NROW(data_n)),breaks=cv_folds,labels=FALSE)
			list = 1:cv_folds
		}
		else{
			s = 1:NROW(data)
			data_n = data
			genes_n = genes
			env_n = env
			id = folds[[j]]
			list = unique(id)
		}
		pred=c()
		y_test=c()

		for (i in list){
			# Train and test datasets
			data_train = subset(data_n, id != i, drop = FALSE)
			genes_train = subset(genes_n, id != i, drop = FALSE)
	 		env_train = subset(env_n, id != i, drop = FALSE)
	 		data_test = subset(data_n, id == i, drop = FALSE)
	 		genes_test = subset(genes_n, id == i, drop = FALSE)
	 		env_test = subset(env_n, id == i, drop = FALSE)
	 		y_test_new = data_test[,formula_outcome]

	 		# Fit model and add predictions
	 		fit_train = LEGIT(data=data_train, genes=genes_train, env=env_train, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE, crossover = crossover, crossover_fixed = crossover_fixed)
			pred_new = predict(fit_train, data=data_test,genes=genes_test,env=env_test,type="response", crossover = crossover, crossover_fixed = crossover_fixed)
			pred = c(pred,pred_new)
			y_test = c(y_test, y_test_new)
		}

		# Cross-validated R2
		ssres = sum((pred-y_test)^2)
		sstotal = sum((y_test-mean(y_test))^2)
		R2_cv = c(R2_cv, 1 - ssres/sstotal)
		# Outlier-resistant cross-validation criterion
		L1_cv = c(L1_cv, sum(abs(pred-y_test))/length(pred))
		Huber_index = abs(pred-y_test) > Huber_p
		Huber_cv_err = (((pred-y_test)^2)/2)
		Huber_cv_err[Huber_index] = (Huber_p*abs(pred-y_test)-(Huber_p^2)/2)[Huber_index]
		Huber_cv = c(Huber_cv, sum(Huber_cv_err)/length(pred))

		#Cross-validated confusion matrix and ROC curve
		if (classification){
			roc_curve_n = pROC::roc(y_test,pred)
			roc_curve = append(roc_curve, list(roc_curve_n))
			AUC = c(AUC, pROC::auc(roc_curve_n))
			best_threshold =  rbind(pROC::coords(roc_curve_n, "best"),best_threshold)
		}

		#Residuals (To detect outliers)
		residuals = residuals + scale(pred-y_test)[s]
		if(class(family)=="function") pearson_residuals = pearson_residuals + scale((pred-y_test)/sqrt(family()$variance(pred)))[s]
		else pearson_residuals = pearson_residuals + scale((pred-y_test)/sqrt(family$variance(pred)))[s]
	}
	residuals = residuals/cv_iter
	pearson_residuals = pearson_residuals/cv_iter

	possible_outliers = abs(residuals)>2.5 | abs(pearson_residuals)>2.5
	if (is.null(obs_id)) possible_outliers_data = cbind(rownames(data)[possible_outliers],residuals[possible_outliers],pearson_residuals[possible_outliers])
	else{
		if (!is.null(dim(obs_id))) obs_id = obs_id[possible_outliers,,drop=FALSE]
		else obs_id = obs_id[possible_outliers]
		possible_outliers_data = cbind(obs_id,residuals[possible_outliers],pearson_residuals[possible_outliers])
	}
	if (NCOL(possible_outliers_data)==3){
		if (is.null(obs_id)) colnames(possible_outliers_data) = c("Observation","Standardized_residual","Standardized_pearson_residual")
		else colnames(possible_outliers_data) = c("ID","Standardized_residual","Standardized_pearson_residual")
	}
	else{
		if (!is.null(colnames(obs_id))) colnames(possible_outliers_data) = c(colnames(obs_id),"Standardized_residual","Standardized_pearson_residual")
		else colnames(possible_outliers_data) = c(rep("ID",NCOL(obs_id)),"Standardized_residual","Standardized_pearson_residual")
	}

	if (classification) return(list(R2_cv = R2_cv, Huber_cv = Huber_cv, L1_cv=L1_cv, AUC=AUC, best_threshold=best_threshold, roc_curve = roc_curve, possible_outliers = possible_outliers_data))
	return(list(R2_cv = R2_cv, Huber_cv = Huber_cv, L1_cv=L1_cv, possible_outliers = possible_outliers_data))
}

IMLEGIT_cv = function (data, latent_var, formula, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, classification=FALSE, start_latent_var=NULL, eps=.001, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, id=NULL){
	if (!is.null(ylim)){
		if (!is.numeric(ylim) || length(ylim) !=2) stop("ylim must either be NULL or a numeric vector of size two")
	}

	# Renaming it because there is already an id variable
	if (!is.null(id)) obs_id = id
	else obs_id = NULL
	id = NULL

	# Setting up latent_var and checks
	if (class(latent_var)!="list") stop("latent_var must be a list of datasets")
	k = length(latent_var)
	if (k==0) stop("latent_var cannot be an empty list")
	if (is.null(names(latent_var))){
		cat("You have not specified names for the latent variables, assigning names to latent_var is highly recommended to prevent confusion. For now, they will be named L1, L2, ...\n")
		names(latent_var) = paste0("L",1:k)
	}
	for (i in 1:k){
		latent_var[[i]] = as.matrix(data.frame(latent_var[[i]],fix.empty.names=FALSE))
		if (sum(colnames(latent_var[[i]])=="") > 0){
			if (print) cat(paste0("You have not specified column names for certain elements in ",names(latent_var)[i], ", elements of this latent variable will be named ",names(latent_var)[i],1,", ",names(latent_var)[i],2," ...\n"))
			colnames(latent_var[[i]]) = paste0(names(latent_var)[i],1:NCOL(latent_var[[i]]))
		}
	}

	# More checks
	if (maxiter <= 0) warning("maxiter must be > 0")
	if (k > 1) for (i in 1:(k-1)) if (NROW(latent_var[[i]]) != NROW(latent_var[[i+1]])) stop("Some datasets in latent_var don't have the same number of observations")
	if(!is.null(start_latent_var)){
		if (class(start_latent_var)!="list") stop("start_latent_var must be a lit of vectors (or NULL)")
		if (k!=length(start_latent_var)) stop("start_latent_var must have the same size as latent_var")
		for (i in 1:k){
			if (!is.null(latent_var[[i]])){
				if (NCOL(latent_var[[i]])!=length(start_latent_var[[i]])) stop("All elements of start_latent_var must either be NULL or have the same length as the number of the elements in its associated latent variable")
			}
		}
	}
	if (class(data) != "data.frame" && class(data) != "matrix") stop("data must be a data.frame")

	# getting right formats
	# Retaining only the needed variables from the dataset (need to set elements in latent_var for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	for (i in 1:k) data[,names(latent_var)[i]] = 0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	formula = stats::as.formula(formula)

	# Error message about factors
	if (sum(apply(data,2,is.numeric)) != NCOL(data)) stop("All variables used must be numeric, factors are not allowed. Please dummy code all categorical variables inside your datasets (data, latent_var[[1]], latent_var[[2]], ...)")
	for (i in 1:k) if (sum(apply(latent_var[[i]],2,is.numeric)) != NCOL(latent_var[[i]])) stop("All variables used must be numeric, factors are not allowed. Please dummy code all categorical variables inside your datasets (data, latent_var[[1]], latent_var[[2]], ...)")

	# remove missing data
	comp = stats::complete.cases(data,latent_var[[1]])
	if (k > 1) for (i in 2:k) comp = comp & stats::complete.cases(latent_var[[i]])
	data = data[comp,, drop=FALSE]
	for (i in 1:k) latent_var[[i]] = latent_var[[i]][comp,, drop=FALSE]
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	formula_outcome = get.vars(formula)[1]
	R2_cv = c()
	Huber_cv = c()
	L1_cv = c()
	AUC = c()
	best_threshold = c()
	roc_curve = list()
	residuals = rep(0,dim(data)[1])
	pearson_residuals = rep(0,dim(data)[1])

	if (!is.null(folds)) cv_iter = length(folds)	

	for (j in 1:cv_iter){
		if (!is.null(seed)) set.seed(seed*j)
		# Folds
		if (is.null(folds)){
			s = sample(NROW(data))
			data_n = data[s,, drop=FALSE]
			latent_var_new = latent_var
			for (l in 1:k) latent_var_new[[l]] = latent_var_new[[l]][s,, drop=FALSE]
			id = cut(seq(1,NROW(data_n)),breaks=cv_folds,labels=FALSE)
			list = 1:cv_folds
		}
		else{
			s = 1:NROW(data)
			data_n = data
			latent_var_new = latent_var
			id = folds[[j]]
			list = unique(id)
		}
		pred=c()
		y_test=c()

		for (i in list){
			# Train and test datasets
			data_train = subset(data_n, id != i, drop = FALSE)
			latent_var_train = latent_var_new
			for (l in 1:k) latent_var_train[[l]] = subset(latent_var_new[[l]], id != i, drop = FALSE)
	 		data_test = subset(data_n, id == i, drop = FALSE)
	 		latent_var_test = latent_var_new
	 		for (l in 1:k) latent_var_test[[l]] = subset(latent_var_new[[l]], id == i, drop = FALSE)
	 		y_test_new = data_test[,formula_outcome]

	 		# Fit model and add predictions
	 		fit_train = IMLEGIT(data=data_train, latent_var=latent_var_train, formula=formula, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
			pred_new = predict(fit_train, data=data_test,latent_var=latent_var_test,type="response")
			pred = c(pred,pred_new)
			y_test = c(y_test, y_test_new)
		}

		# Cross-validated R2
		ssres = sum((pred-y_test)^2)
		sstotal = sum((y_test-mean(y_test))^2)
		R2_cv = c(R2_cv, 1 - ssres/sstotal)
		# Outlier-resistant cross-validation criterion
		L1_cv = c(L1_cv, sum(abs(pred-y_test))/length(pred))
		Huber_index = abs(pred-y_test) > Huber_p
		Huber_cv_err = (((pred-y_test)^2)/2)
		Huber_cv_err[Huber_index] = (Huber_p*abs(pred-y_test)-(Huber_p^2)/2)[Huber_index]
		Huber_cv = c(Huber_cv, sum(Huber_cv_err)/length(pred))

		#Cross-validated confusion matrix and ROC curve
		if (classification){
			roc_curve_n = pROC::roc(y_test,pred)
			roc_curve = append(roc_curve, list(roc_curve_n))
			AUC = c(AUC, pROC::auc(roc_curve_n))
			best_threshold =  rbind(pROC::coords(roc_curve_n, "best"),best_threshold)
		}

		#Residuals (To detect outliers)
		residuals = residuals + scale(pred-y_test)[s]
		if(class(family)=="function") pearson_residuals = pearson_residuals + scale((pred-y_test)/sqrt(family()$variance(pred)))[s]
		else pearson_residuals = pearson_residuals + scale((pred-y_test)/sqrt(family$variance(pred)))[s]
	}
	residuals = residuals/cv_iter
	pearson_residuals = pearson_residuals/cv_iter

	possible_outliers = abs(residuals)>2.5 | abs(pearson_residuals)>2.5
	if (is.null(obs_id)) possible_outliers_data = cbind(rownames(data)[possible_outliers],residuals[possible_outliers],pearson_residuals[possible_outliers])
	else{
		if (!is.null(dim(obs_id))) obs_id = obs_id[possible_outliers,,drop=FALSE]
		else obs_id = obs_id[possible_outliers]
		possible_outliers_data = cbind(obs_id,residuals[possible_outliers],pearson_residuals[possible_outliers])
	}
	if (NCOL(possible_outliers_data)==3){
		if (is.null(obs_id)) colnames(possible_outliers_data) = c("Observation","Standardized_residual","Standardized_pearson_residual")
		else colnames(possible_outliers_data) = c("ID","Standardized_residual","Standardized_pearson_residual")
	}
	else{
		if (!is.null(colnames(obs_id))) colnames(possible_outliers_data) = c(colnames(obs_id),"Standardized_residual","Standardized_pearson_residual")
		else colnames(possible_outliers_data) = c(rep("ID",NCOL(obs_id)),"Standardized_residual","Standardized_pearson_residual")
	}

	if (classification) return(list(R2_cv = R2_cv, Huber_cv = Huber_cv, L1_cv=L1_cv, AUC=AUC, best_threshold=best_threshold, roc_curve = roc_curve, possible_outliers = possible_outliers_data))
	return(list(R2_cv = R2_cv, Huber_cv = Huber_cv, L1_cv=L1_cv, possible_outliers = possible_outliers_data))
}

forward_step = function(empty_start_dataset, fit, data, formula, interactive_mode=FALSE, genes_current=NULL, env_current=NULL, genes_toadd=NULL, env_toadd=NULL, search="genes", search_criterion="AIC", p_threshold = .20, exclude_worse_AIC=TRUE, max_steps = 100, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, classification=FALSE, start_genes=NULL, start_env=NULL, eps=.01, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, print=TRUE){
	# How much genes or env to add
	if (search=="genes") elements_N = NCOL(genes_toadd)
	if (search=="env") elements_N = NCOL(env_toadd)
	if (elements_N == 0){
		if (search=="genes" && print) cat("No gene added\n")
		if (search=="env" && print) cat("No environment added\n")
		return(NULL)
	}
	# Vector which says which extra variables are "good" (worth exploring)
	good = rep(TRUE, elements_N)
	# Vector which says how much the criterion changed from including the variable
	criterion_before = rep(NA, elements_N)
	criterion_after = rep(NA, elements_N)
	criterion_diff = rep(NA, elements_N)
	# In interactive model, we must keep track of every AIC, BIC, p-value, Cross-validated R2 and AUC to show the user at every iteration
	if (interactive_mode){
		if (search=="genes") interactive = data.frame(variable=colnames(genes_toadd), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N))
		if (search=="env") interactive = data.frame(variable=colnames(env_toadd), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N))
		if (search=="genes" && (search_criterion == "cv" || search_criterion == "cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1")) interactive = data.frame(variable=colnames(genes_toadd), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N),cv_R2_old=rep(NA, elements_N),cv_R2_new=rep(NA, elements_N),cv_AUC_old=rep(NA, elements_N),cv_AUC_new=rep(NA, elements_N),cv_Huber_old=rep(NA, elements_N),cv_Huber_new=rep(NA, elements_N),cv_L1_old=rep(NA, elements_N),cv_L1_new=rep(NA, elements_N))
		if (search=="env" && (search_criterion == "cv" || search_criterion == "cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1")) interactive = data.frame(variable=colnames(env_toadd), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N),cv_R2_old=rep(NA, elements_N),cv_R2_new=rep(NA, elements_N),cv_AUC_old=rep(NA, elements_N),cv_AUC_new=rep(NA, elements_N),cv_Huber_old=rep(NA, elements_N),cv_Huber_new=rep(NA, elements_N),cv_L1_old=rep(NA, elements_N),cv_L1_new=rep(NA, elements_N))
	}
	# Complete dataset
	if (NCOL(genes_current)==0) genes_current_nomiss = NULL
	else genes_current_nomiss = genes_current
	if (NCOL(env_current)==0) env_current_nomiss = NULL
	else env_current_nomiss = env_current
	comp_without = stats::complete.cases(data,genes_current_nomiss,env_current_nomiss)
	# Non-cross-validated models
	for (j in 1:elements_N){
		if (search=="genes") fit_with = LEGIT(data=data, genes=cbind(genes_current,genes_toadd[,j,drop=FALSE]), env=env_current, formula=formula, start_genes=c(start_genes,0), start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
		if (search=="env") fit_with = LEGIT(data=data, genes=genes_current, env=cbind(env_current,env_toadd[,j,drop=FALSE]), formula=formula, start_genes=start_genes, start_env=c(start_env,0), eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
		if (search=="genes"){
			p_value = stats::coef(summary(fit_with)$fit_genes)[,4]
			good[j] = p_value[length(p_value)] <= p_threshold
		}
		if (search=="env"){
			p_value = stats::coef(summary(fit_with)$fit_env)[,4]
			good[j] = p_value[length(p_value)] <= p_threshold
		}
		if (empty_start_dataset) fit_without = NULL
		else if (fit$fit_main$df.null != fit_with$fit_main$df.null){
			if (search=="genes") comp = stats::complete.cases(data,genes_current,genes_toadd[,j,drop=FALSE],env_current)
			if (search=="env") comp = stats::complete.cases(data,genes_current,env_toadd[,j,drop=FALSE],env_current)
			fit_without = LEGIT(data=data[comp,,drop=FALSE], genes=genes_current[comp,,drop=FALSE], env=env_current[comp,,drop=FALSE], formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
		}
		else fit_without = fit

		if (exclude_worse_AIC && !empty_start_dataset){
			if (fit_with$true_model_parameters$AIC <= fit_without$true_model_parameters$AIC) good[j]=good[j] && TRUE
		}
		if (search_criterion=="AIC"){
			if (empty_start_dataset) criterion_before[j] = Inf
			else criterion_before[j] = fit_without$true_model_parameters$AIC
			criterion_after[j] = fit_with$true_model_parameters$AIC
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		}
		if (search_criterion=="AICc"){
			if (empty_start_dataset) criterion_before[j] = Inf
			else criterion_before[j] = fit_without$true_model_parameters$AICc
			criterion_after[j] = fit_with$true_model_parameters$AICc
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		}
		if (search_criterion=="BIC"){
			if (empty_start_dataset) criterion_before[j] = Inf
			else criterion_before[j] = fit_without$true_model_parameters$BIC
			criterion_after[j] = fit_with$true_model_parameters$BIC
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		} 
		# Keep p-value, AIC and BIC if in interactive model
		if (interactive_mode){
			if (empty_start_dataset) interactive$N_old[j] = NA
			else interactive$N_old[j] = sum(comp_without)
			interactive$N_new[j] = fit_with$true_model_parameters$df.residual + fit_with$true_model_parameters$rank

			interactive$p_value[j] = round(p_value[length(p_value)],6)

			if (empty_start_dataset) interactive$AIC_old[j] = Inf
			else interactive$AIC_old[j] = fit_without$true_model_parameters$AIC
			interactive$AIC_new[j] = fit_with$true_model_parameters$AIC

			if (empty_start_dataset) interactive$AICc_old[j] = Inf
			else interactive$AICc_old[j] = fit_without$true_model_parameters$AICc
			interactive$AICc_new[j] = fit_with$true_model_parameters$AICc

			if (empty_start_dataset) interactive$BIC_old[j] = Inf
			else interactive$BIC_old[j] = fit_without$true_model_parameters$BIC
			interactive$BIC_new[j] = fit_with$true_model_parameters$BIC
		}
	}
	if (sum(good)==0){
		if (search=="genes" && print) cat("No gene added\n")
		if (search=="env" && print) cat("No environment added\n")
		return(NULL)
	}
	# Only do this if NOT in interactive mode, otherwise at the end of the algorithm, we show the the data.frame interactive and let the user choose
	if (!interactive_mode){
		if (search_criterion=="AIC" || search_criterion=="AICc" || search_criterion=="BIC"){
			if (min(criterion_diff,na.rm=TRUE) > 0){
				if (search=="genes" && print) cat("No gene added\n")
				if (search=="env" && print) cat("No environment added\n")
				return(NULL)
			}
			if (empty_start_dataset) best_var = which.min(criterion_after)
			else best_var = which.min(criterion_diff)
			if (search=="genes"){
				if (print) cat(paste0("Adding gene: ",colnames(genes_toadd)[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after = ",round(criterion_after[best_var],5),")\n"))
				genes_current = cbind(genes_current, genes_toadd[,best_var, drop=FALSE])
				genes_toadd = genes_toadd[,-best_var, drop=FALSE]
			}
			if (search=="env"){
				if (print) cat(paste0("Adding environment: ",colnames(env_toadd)[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after = ",round(criterion_after[best_var],5),")\n"))
				env_current = cbind(env_current, env_toadd[,best_var, drop=FALSE])
				env_toadd = env_toadd[,-best_var, drop=FALSE]						
			}
		}
	}
	# Cross-validated models
	if (search_criterion == "cv" || search_criterion == "cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1"){
		# Dropping variables with p < threshold and worse AIC
		if (interactive_mode) interactive = interactive[good,]
		if (search=="genes") genes_toadd = genes_toadd[,good, drop=FALSE]
		if (search=="env") env_toadd = env_toadd[,good, drop=FALSE]
		if (search=="genes") elements_N = NCOL(genes_toadd)
		if (search=="env") elements_N = NCOL(env_toadd)
		# Vector which says how much the criterion changed from including the variable
		criterion_before = rep(NA, elements_N)
		criterion_after = rep(NA, elements_N)
		criterion_diff = rep(NA, elements_N)
		# Set seed
		if (!is.null(seed)) current_seed = seed
		else current_seed = NULL
		if (!empty_start_dataset) fit_cv = LEGIT_cv(data=data, genes=genes_current, env=env_current, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
		for (j in 1:elements_N){
			if (search=="genes") fit_cv_with = LEGIT_cv(data=data, genes=cbind(genes_current,genes_toadd[,j,drop=FALSE]), env=env_current, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=c(start_genes,0), start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
			if (search=="env") fit_cv_with = LEGIT_cv(data=data, genes=genes_current, env=cbind(env_current,env_toadd[,j,drop=FALSE]), formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=c(start_env,0), eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
			if (empty_start_dataset) fit_cv_without = NULL
			else{
				if (search=="genes") comp_with = stats::complete.cases(data,genes_current,genes_toadd[,j,drop=FALSE],env_current)
				if (search=="env") comp_with = stats::complete.cases(data,genes_current,env_toadd[,j,drop=FALSE],env_current)
				if (sum(comp_without) != sum(comp_with)) fit_cv_without = LEGIT_cv(data=data[comp_with,,drop=FALSE], genes=genes_current[comp_with,,drop=FALSE], env=env_current[comp_with,,drop=FALSE], formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
				else fit_cv_without = fit_cv
			}
			if (search_criterion=="cv"){
				if (empty_start_dataset) criterion_before[j] = 0
				else criterion_before[j] = mean(fit_cv_without$R2_cv)
				criterion_after[j] = mean(fit_cv_with$R2_cv)
				criterion_diff[j] = criterion_after[j] - criterion_before[j]
			} 
			if (search_criterion=="cv_Huber"){
				if (empty_start_dataset) criterion_before[j] = Inf
				else criterion_before[j] = mean(fit_cv_without$Huber_cv)
				criterion_after[j] = mean(fit_cv_with$Huber_cv)
				criterion_diff[j] = criterion_after[j] - criterion_before[j]
			} 
			if (search_criterion=="cv_L1"){
				if (empty_start_dataset) criterion_before[j] = Inf
				else criterion_before[j] = mean(fit_cv_without$L1_cv)
				criterion_after[j] = mean(fit_cv_with$L1_cv)
				criterion_diff[j] = criterion_after[j] - criterion_before[j]
			} 
			if (search_criterion=="cv_AUC"){
				if (empty_start_dataset) criterion_before[j] = 0
				else criterion_before[j] = mean(fit_cv_without$AUC)
				criterion_after[j] = mean(fit_cv_with$AUC)
				criterion_diff[j] = criterion_after[j] - criterion_before[j]
			}
			# Keep cross-validation R2 and AUC in interactive model
			if (interactive_mode){
				if (empty_start_dataset) interactive$cv_R2_old[j] = 0
				else interactive$cv_R2_old[j] = mean(fit_cv_without$R2_cv)
				interactive$cv_R2_new[j] = mean(fit_cv_with$R2_cv)

				if (empty_start_dataset) interactive$cv_L1_old[j] = Inf
				else interactive$cv_L1_old[j] = mean(fit_cv_without$L1_cv)
				interactive$cv_L1_new[j] = mean(fit_cv_with$L1_cv)

				if (empty_start_dataset) interactive$cv_Huber_old[j] = Inf
				else interactive$cv_Huber_old[j] = mean(fit_cv_without$Huber_cv)
				interactive$cv_Huber_new[j] = mean(fit_cv_with$Huber_cv)

				if(classification){
					if (empty_start_dataset) interactive$cv_AUC_old[j] = 0
					else interactive$cv_AUC_old[j] = mean(fit_cv_without$AUC)
					interactive$cv_AUC_new[j] = mean(fit_cv_with$AUC)
				}
			}
		}
		# Only do this if NOT in interactive mode, otherwise at the end of the algorithm, we show the the data.frame interactive and let the user choose
		if (!interactive_mode){
			if ((max(criterion_diff,na.rm=TRUE) < 0 && !(search_criterion=="cv_Huber" || search_criterion=="cv_L1")) || (min(criterion_diff,na.rm=TRUE) > 0 && (search_criterion=="cv_Huber" || search_criterion=="cv_L1"))){
				if (search=="genes" && print) cat("No gene added\n")
				if (search=="env" && print) cat("No environment added\n")
				return(NULL)
			}
			if (search_criterion=="cv_Huber" || search_criterion=="cv_L1") best_var = which.min(criterion_diff)
			else best_var = which.max(criterion_diff)
			if (search=="genes"){
				if (print) cat(paste0("Adding gene: ",colnames(genes_toadd)[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after = ",round(criterion_after[best_var],5),")\n"))
				genes_current = cbind(genes_current, genes_toadd[,best_var, drop=FALSE])
				genes_toadd = genes_toadd[,-best_var, drop=FALSE]
			}
			if (search=="env"){
				if (print) cat(paste0("Adding environment: ",colnames(env_toadd)[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after = ",round(criterion_after[best_var],5),")\n"))
				env_current = cbind(env_current, env_toadd[,best_var, drop=FALSE])
				env_toadd = env_toadd[,-best_var, drop=FALSE]				
			}
		}
	}
	if (interactive_mode){
		interactive_n = min(5,elements_N)
		if (search_criterion=="AIC" && empty_start_dataset) neworder = order(interactive$AIC_new, decreasing=FALSE)
		if (search_criterion=="AIC" && !empty_start_dataset) neworder = order(interactive$AIC_new - interactive$AIC_old, decreasing=FALSE)
		if (search_criterion=="AICc" && empty_start_dataset) neworder = order(interactive$AICc_new, decreasing=FALSE)
		if (search_criterion=="AICc" && !empty_start_dataset) neworder = order(interactive$AICc_new - interactive$AICc_old, decreasing=FALSE)
		if (search_criterion=="BIC" && empty_start_dataset) neworder = order(interactive$BIC_new, decreasing=FALSE)
		if (search_criterion=="BIC" && !empty_start_dataset) neworder = order(interactive$BIC_new - interactive$BIC_old, decreasing=FALSE)
		if (search_criterion=="cv" && empty_start_dataset) neworder = order(interactive$cv_R2_new, decreasing=TRUE)
		if (search_criterion=="cv" && !empty_start_dataset) neworder = order(interactive$cv_R2_new - interactive$cv_R2_old, decreasing=TRUE)
		if (search_criterion=="cv_Huber" && empty_start_dataset) neworder = order(interactive$cv_Huber_new, decreasing=FALSE)
		if (search_criterion=="cv_Huber" && !empty_start_dataset) neworder = order(interactive$cv_Huber_new - interactive$cv_Huber_old, decreasing=FALSE)
		if (search_criterion=="cv_L1" && empty_start_dataset) neworder = order(interactive$cv_L1_new, decreasing=FALSE)
		if (search_criterion=="cv_L1" && !empty_start_dataset) neworder = order(interactive$cv_L1_new - interactive$cv_L1_old, decreasing=FALSE)
		if (search_criterion=="cv_AUC" && empty_start_dataset) neworder = order(interactive$cv_AUC_new, decreasing=TRUE)
		if (search_criterion=="cv_AUC" && !empty_start_dataset) neworder = order(interactive$cv_AUC_new - interactive$cv_AUC_old, decreasing=TRUE)
		interactive = interactive[neworder[1:interactive_n],]
		rownames(interactive)=1:interactive_n
		interactive = format(interactive,scientific=FALSE)
		print(interactive)
		if (search=="genes") input_user = readline(prompt="Enter the index of the gene to be added: ")
		if (search=="env") input_user = readline(prompt="Enter the index of the environment to be added: ")
		if (sum(input_user == rownames(interactive))==0){
			if (search=="genes" && print) cat("No gene added\n")
			if (search=="env" && print) cat("No environment added\n")
			return(NULL)
		}
		best_var = neworder[1:interactive_n][input_user == rownames(interactive)]
		if (search=="genes"){
			if (print) cat(paste0("Adding gene: ",colnames(genes_toadd)[best_var],"\n"))
			genes_current = cbind(genes_current, genes_toadd[,best_var, drop=FALSE])
			genes_toadd = genes_toadd[,-best_var, drop=FALSE]
		}
		if (search=="env"){
			if (print) cat(paste0("Adding environment: ",colnames(env_toadd)[best_var],"\n"))
			env_current = cbind(env_current, env_toadd[,best_var, drop=FALSE])
			env_toadd = env_toadd[,-best_var, drop=FALSE]				
		}
	}
	# Updated model and coefficients
	if (search=="genes") start_genes=c(start_genes,0)
	if (search=="env") start_env=c(start_env,0)
	fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
	start_genes = stats::coef(fit$fit_genes)
	start_env = stats::coef(fit$fit_env)
	if (search=="genes") return(list(fit=fit, start_genes=start_genes,start_env=start_env,genes_current=genes_current,genes_toadd=genes_toadd))
	if (search=="env") return(list(fit=fit, start_genes=start_genes,start_env=start_env,env_current=env_current,env_toadd=env_toadd))
}

forward_step_IM = function(empty_start_dataset, fit, data, formula, interactive_mode=FALSE, latent_var_current=NULL, latent_var_toadd=NULL, search=NULL, search_criterion="AIC", p_threshold = .20, exclude_worse_AIC=TRUE, max_steps = 100, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, classification=FALSE, start_latent_var=start_latent_var, eps=.01, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, print=TRUE){
	k = length(latent_var_current)
	# How much genes or env to add
	elements_N = NCOL(latent_var_toadd[[search]])
	if (elements_N == 0){
		if (print) cat(paste0("No element from ", names(latent_var_current)[search]," was added\n"))
		return(NULL)
	}
	# Vector which says which extra variables are "good" (worth exploring)
	good = rep(TRUE, elements_N)
	# Vector which says how much the criterion changed from including the variable
	criterion_before = rep(NA, elements_N)
	criterion_after = rep(NA, elements_N)
	criterion_diff = rep(NA, elements_N)
	# In interactive model, we must keep track of every AIC, BIC, p-value, Cross-validated R2 and AUC to show the user at every iteration
	if (interactive_mode){
		if (search_criterion == "cv" || search_criterion == "cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1") interactive = data.frame(variable=colnames(latent_var_toadd[[search]]), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N),cv_R2_old=rep(NA, elements_N),cv_R2_new=rep(NA, elements_N),cv_AUC_old=rep(NA, elements_N),cv_AUC_new=rep(NA, elements_N),cv_Huber_old=rep(NA, elements_N),cv_Huber_new=rep(NA, elements_N),cv_L1_old=rep(NA, elements_N),cv_L1_new=rep(NA, elements_N))
		else interactive = data.frame(variable=colnames(latent_var_toadd[[search]]), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N))
	}	
	# Complete dataset
	if (NCOL(latent_var_current[[1]])==0) latent_var_current_nomiss = NULL
	else latent_var_current_nomiss = latent_var_current[[1]]
	comp_without = stats::complete.cases(data,latent_var_current_nomiss)
	if (k > 1) for (i in 2:k){
		if (NCOL(latent_var_current[[i]])==0) latent_var_current_nomiss = NULL
		else latent_var_current_nomiss = latent_var_current[[i]]
		comp_without = comp_without & stats::complete.cases(latent_var_current_nomiss)
	}

	# Non-cross-validated models
	for (j in 1:elements_N){
		# Running IMLEGIT with new variable
		latent_var_new = latent_var_current
		if (is.null(latent_var_current[[search]])) latent_var_new[[search]] = latent_var_toadd[[search]][,j,drop=FALSE]
		else latent_var_new[[search]] = cbind(latent_var_current[[search]],latent_var_toadd[[search]][,j,drop=FALSE])
		start_latent_var_new = start_latent_var
		start_latent_var_new[[search]] = c(start_latent_var[[search]],0)
		fit_with = IMLEGIT(data=data, latent_var=latent_var_new, formula=formula, start_latent_var=start_latent_var_new, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
		# p-values
		p_value = stats::coef(summary(fit_with)[[search+1]])[,4]
		good[j] = p_value[length(p_value)] <= p_threshold
		if (empty_start_dataset) fit_without = NULL
		else if (fit$fit_main$df.null != fit_with$fit_main$df.null){
			# Removing observations missing the new variable and rerunnning model without the variable (Note : with and without is confusing here)
			comp_with = comp_without & stats::complete.cases(latent_var_toadd[[search]][,j,drop=FALSE])
			data_with = data[comp_with,, drop=FALSE]
			latent_var_with = latent_var_current
			for (i in 1:k) latent_var_with[[i]] = latent_var_current[[i]][comp_with,, drop=FALSE]
			fit_without = IMLEGIT(data=data_with, latent_var=latent_var_with, formula=formula, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
		}
		else fit_without = fit

		if (exclude_worse_AIC && !empty_start_dataset){
			if (fit_with$true_model_parameters$AIC <= fit_without$true_model_parameters$AIC) good[j]=good[j] && TRUE
		}
		if (search_criterion=="AIC"){
			if (empty_start_dataset) criterion_before[j] = Inf
			else criterion_before[j] = fit_without$true_model_parameters$AIC
			criterion_after[j] = fit_with$true_model_parameters$AIC
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		}
		if (search_criterion=="AICc"){
			if (empty_start_dataset) criterion_before[j] = Inf
			else criterion_before[j] = fit_without$true_model_parameters$AICc
			criterion_after[j] = fit_with$true_model_parameters$AICc
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		}
		if (search_criterion=="BIC"){
			if (empty_start_dataset) criterion_before[j] = Inf
			else criterion_before[j] = fit_without$true_model_parameters$BIC
			criterion_after[j] = fit_with$true_model_parameters$BIC
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		} 
		# Keep p-value, AIC and BIC if in interactive model
		if (interactive_mode){
			if (empty_start_dataset) interactive$N_old[j] = NA
			else interactive$N_old[j] = sum(comp_without)
			interactive$N_new[j] = fit_with$true_model_parameters$df.residual + fit_with$true_model_parameters$rank

			interactive$p_value[j] = round(p_value[length(p_value)],6)

			if (empty_start_dataset) interactive$AIC_old[j] = Inf
			else interactive$AIC_old[j] = fit_without$true_model_parameters$AIC
			interactive$AIC_new[j] = fit_with$true_model_parameters$AIC

			if (empty_start_dataset) interactive$AICc_old[j] = Inf
			else interactive$AICc_old[j] = fit_without$true_model_parameters$AICc
			interactive$AICc_new[j] = fit_with$true_model_parameters$AICc

			if (empty_start_dataset) interactive$BIC_old[j] = Inf
			else interactive$BIC_old[j] = fit_without$true_model_parameters$BIC
			interactive$BIC_new[j] = fit_with$true_model_parameters$BIC
		}
	}
	if (sum(good)==0){
		if (print) cat(paste0("No element from ", names(latent_var_current)[search]," was added\n"))
		return(NULL)
	}
	# Only do this if NOT in interactive mode, otherwise at the end of the algorithm, we show the the data.frame interactive and let the user choose
	if (!interactive_mode){
		if (search_criterion=="AIC" || search_criterion=="AICc" || search_criterion=="BIC"){
			if (min(criterion_diff,na.rm=TRUE) > 0){
				if (print) cat(paste0("No element from ", names(latent_var_current)[search]," was added\n"))
				return(NULL)
			}
			if (empty_start_dataset) best_var = which.min(criterion_after)
			else best_var = which.min(criterion_diff)
			if (print) cat(paste0("Adding element from ", names(latent_var_current)[search], ": ",colnames(latent_var_toadd[[search]])[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after = ",round(criterion_after[best_var],5),")\n"))
			latent_var_current[[search]] = cbind(latent_var_current[[search]], latent_var_toadd[[search]][,best_var, drop=FALSE])
			latent_var_toadd[[search]] = latent_var_toadd[[search]][,-best_var, drop=FALSE]						
		}
	}
	# Cross-validated models
	if (search_criterion == "cv" || search_criterion == "cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1"){
		# Dropping variables with p < threshold and worse AIC
		if (interactive_mode) interactive = interactive[good,]
		latent_var_toadd[[search]] = latent_var_toadd[[search]][,good, drop=FALSE]
		elements_N = NCOL(latent_var_toadd[[search]])
		# Vector which says how much the criterion changed from including the variable
		criterion_before = rep(NA, elements_N)
		criterion_after = rep(NA, elements_N)
		criterion_diff = rep(NA, elements_N)
		# Set seed
		if (!is.null(seed)) current_seed = seed
		else current_seed = NULL
		if (!empty_start_dataset) fit_cv = IMLEGIT_cv(data=data, latent_var=latent_var_current, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
		for (j in 1:elements_N){
			# Running IMLEGIT with new variable
			latent_var_new = latent_var_current
			if (is.null(latent_var_current[[search]])) latent_var_new[[search]] = latent_var_toadd[[search]][,j,drop=FALSE]
			latent_var_new[[search]] = cbind(latent_var_current[[search]],latent_var_toadd[[search]][,j,drop=FALSE])
			start_latent_var_new = start_latent_var
			start_latent_var_new[[search]] = c(start_latent_var[[search]],0)
			fit_cv_with = IMLEGIT_cv(data=data, latent_var=latent_var_new, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var_new, eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
			if (empty_start_dataset) fit_cv_without = NULL
			else{
				# If new variable has new missing data: Remove observations missing the new variable and rerun model without the variable (Note : with and without is confusing here)
				comp_with = comp_without & stats::complete.cases(latent_var_toadd[[search]][,j,drop=FALSE])
				if (sum(comp_without) != sum(comp_with)){
					data_with = data[comp_with,, drop=FALSE]
					latent_var_with = latent_var_current
					for (i in 1:k) latent_var_with[[i]] = latent_var_current[[i]][comp_with,, drop=FALSE]
					fit_cv_without = IMLEGIT_cv(data=data_with, latent_var=latent_var_with, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
				}
				else fit_cv_without = fit_cv
			}
			if (search_criterion=="cv"){
				if (empty_start_dataset) criterion_before[j] = 0
				else criterion_before[j] = mean(fit_cv_without$R2_cv)
				criterion_after[j] = mean(fit_cv_with$R2_cv)
				criterion_diff[j] = criterion_after[j] - criterion_before[j]
			} 
			if (search_criterion=="cv_Huber"){
				if (empty_start_dataset) criterion_before[j] = Inf
				else criterion_before[j] = mean(fit_cv_without$Huber_cv)
				criterion_after[j] = mean(fit_cv_with$Huber_cv)
				criterion_diff[j] = criterion_after[j] - criterion_before[j]
			} 
			if (search_criterion=="cv_L1"){
				if (empty_start_dataset) criterion_before[j] = Inf
				else criterion_before[j] = mean(fit_cv_without$L1_cv)
				criterion_after[j] = mean(fit_cv_with$L1_cv)
				criterion_diff[j] = criterion_after[j] - criterion_before[j]
			} 
			if (search_criterion=="cv_AUC"){
				if (empty_start_dataset) criterion_before[j] = 0
				else criterion_before[j] = mean(fit_cv_without$AUC)
				criterion_after[j] = mean(fit_cv_with$AUC)
				criterion_diff[j] = criterion_after[j] - criterion_before[j]
			}
			# Keep cross-validation R2 and AUC in interactive model
			if (interactive_mode){
				if (empty_start_dataset) interactive$cv_R2_old[j] = 0
				else interactive$cv_R2_old[j] = mean(fit_cv_without$R2_cv)
				interactive$cv_R2_new[j] = mean(fit_cv_with$R2_cv)

				if (empty_start_dataset) interactive$cv_L1_old[j] = Inf
				else interactive$cv_L1_old[j] = mean(fit_cv_without$L1_cv)
				interactive$cv_L1_new[j] = mean(fit_cv_with$L1_cv)

				if (empty_start_dataset) interactive$cv_Huber_old[j] = Inf
				else interactive$cv_Huber_old[j] = mean(fit_cv_without$Huber_cv)
				interactive$cv_Huber_new[j] = mean(fit_cv_with$Huber_cv)

				if(classification){
					if (empty_start_dataset) interactive$cv_AUC_old[j] = 0
					else interactive$cv_AUC_old[j] = mean(fit_cv_without$AUC)
					interactive$cv_AUC_new[j] = mean(fit_cv_with$AUC)
				}
			}
		}
		# Only do this if NOT in interactive mode, otherwise at the end of the algorithm, we show the the data.frame interactive and let the user choose
		if (!interactive_mode){
			if ((max(criterion_diff,na.rm=TRUE) < 0 && !(search_criterion=="cv_Huber" || search_criterion=="cv_L1")) || (min(criterion_diff,na.rm=TRUE) > 0 && (search_criterion=="cv_Huber" || search_criterion=="cv_L1"))){
				if (print) cat(paste0("No element from ", names(latent_var_current)[search]," was added\n"))
				return(NULL)
			}
			if (search_criterion=="cv_Huber" || search_criterion=="cv_L1") best_var = which.min(criterion_diff)
			else best_var = which.max(criterion_diff)
			if (print) cat(paste0("Adding element from ", names(latent_var_current)[search], ": ",colnames(latent_var_toadd[[search]])[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after = ",round(criterion_after[best_var],5),")\n"))
			latent_var_current[[search]] = cbind(latent_var_current[[search]], latent_var_toadd[[search]][,best_var, drop=FALSE])
			latent_var_toadd[[search]] = latent_var_toadd[[search]][,-best_var, drop=FALSE]							
		}
	}
	if (interactive_mode){
		interactive_n = min(5,elements_N)
		if (search_criterion=="AIC" && empty_start_dataset) neworder = order(interactive$AIC_new, decreasing=FALSE)
		if (search_criterion=="AIC" && !empty_start_dataset) neworder = order(interactive$AIC_new - interactive$AIC_old, decreasing=FALSE)
		if (search_criterion=="AICc" && empty_start_dataset) neworder = order(interactive$AICc_new, decreasing=FALSE)
		if (search_criterion=="AICc" && !empty_start_dataset) neworder = order(interactive$AICc_new - interactive$AICc_old, decreasing=FALSE)
		if (search_criterion=="BIC" && empty_start_dataset) neworder = order(interactive$BIC_new, decreasing=FALSE)
		if (search_criterion=="BIC" && !empty_start_dataset) neworder = order(interactive$BIC_new - interactive$BIC_old, decreasing=FALSE)
		if (search_criterion=="cv" && empty_start_dataset) neworder = order(interactive$cv_R2_new, decreasing=TRUE)
		if (search_criterion=="cv" && !empty_start_dataset) neworder = order(interactive$cv_R2_new - interactive$cv_R2_old, decreasing=TRUE)
		if (search_criterion=="cv_Huber" && empty_start_dataset) neworder = order(interactive$cv_Huber_new, decreasing=FALSE)
		if (search_criterion=="cv_Huber" && !empty_start_dataset) neworder = order(interactive$cv_Huber_new - interactive$cv_Huber_old, decreasing=FALSE)
		if (search_criterion=="cv_L1" && empty_start_dataset) neworder = order(interactive$cv_L1_new, decreasing=FALSE)
		if (search_criterion=="cv_L1" && !empty_start_dataset) neworder = order(interactive$cv_L1_new - interactive$cv_L1_old, decreasing=FALSE)
		if (search_criterion=="cv_AUC" && empty_start_dataset) neworder = order(interactive$cv_AUC_new, decreasing=TRUE)
		if (search_criterion=="cv_AUC" && !empty_start_dataset) neworder = order(interactive$cv_AUC_new - interactive$cv_AUC_old, decreasing=TRUE)
		interactive = interactive[neworder[1:interactive_n],]
		rownames(interactive)=1:interactive_n
		interactive = format(interactive,scientific=FALSE)
		print(interactive)
		input_user = readline(prompt=paste0("Enter the index of the element from ", names(latent_var_current)[search], " to be added: "))
		if (sum(input_user == rownames(interactive))==0){
			if (print) cat(paste0("No element from ", names(latent_var_current)[search]," was added\n"))
			return(NULL)
		}
		best_var = neworder[1:interactive_n][input_user == rownames(interactive)]
		if (print) cat(paste0("Adding element from ", names(latent_var_current)[search], ": ",colnames(latent_var_toadd[[search]])[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after = ",round(criterion_after[best_var],5),")\n"))
		latent_var_current[[search]] = cbind(latent_var_current[[search]], latent_var_toadd[[search]][,best_var, drop=FALSE])
		latent_var_toadd[[search]] = latent_var_toadd[[search]][,-best_var, drop=FALSE]						

	}
	# Updated model and coefficients
	start_latent_var[[search]] = c(start_latent_var[[search]],0)
	fit = IMLEGIT(data=data, latent_var=latent_var_current, formula=formula, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
	for (i in 1:k) start_latent_var[[i]] = stats::coef(fit$fit_latent_var[[i]])
	return(list(fit=fit, start_latent_var=start_latent_var, latent_var_current=latent_var_current,latent_var_toadd=latent_var_toadd))
}


backward_step = function(fit, data, formula, interactive_mode=FALSE, genes_current=NULL, env_current=NULL, genes_dropped=NULL, env_dropped=NULL, search="genes", search_criterion="AIC", p_threshold = .20, exclude_worse_AIC=TRUE, max_steps = 100, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, classification=FALSE, start_genes=NULL, start_env=NULL, eps=.01, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, print=TRUE){
	# How much genes or env to add
	if (search=="genes") elements_N = NCOL(genes_current)
	if (search=="env") elements_N = NCOL(env_current)
	if (elements_N <= 1){
		if (search=="genes" && print) cat("No gene removed\n")
		if (search=="env" && print) cat("No environment removed\n")
		return(NULL)
	}
	# Vector which says which variables are "good" (worth keeping)
	good = rep(FALSE, elements_N)
	# Vector which says how much the criterion changed from excluding the variable
	criterion_before = rep(NA, elements_N)
	criterion_after = rep(NA, elements_N)
	criterion_diff = rep(NA, elements_N)
	# In interactive model, we must keep track of every AIC, BIC, p-value, Cross-validated R2 and AUC to show the user at every iteration
	if (interactive_mode){
		if (search=="genes") interactive = data.frame(variable=colnames(genes_current), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N))
		if (search=="env") interactive = data.frame(variable=colnames(env_current), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N))
		if (search=="genes" && (search_criterion == "cv" || search_criterion == "cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1")) interactive = data.frame(variable=colnames(genes_current), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N),cv_R2_old=rep(NA, elements_N),cv_R2_new=rep(NA, elements_N),cv_AUC_old=rep(NA, elements_N),cv_AUC_new=rep(NA, elements_N),cv_Huber_old=rep(NA, elements_N),cv_Huber_new=rep(NA, elements_N),cv_L1_old=rep(NA, elements_N),cv_L1_new=rep(NA, elements_N))
		if (search=="env" && (search_criterion == "cv" || search_criterion == "cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1")) interactive = data.frame(variable=colnames(env_current), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N),cv_R2_old=rep(NA, elements_N),cv_R2_new=rep(NA, elements_N),cv_AUC_old=rep(NA, elements_N),cv_AUC_new=rep(NA, elements_N),cv_Huber_old=rep(NA, elements_N),cv_Huber_new=rep(NA, elements_N),cv_L1_old=rep(NA, elements_N),cv_L1_new=rep(NA, elements_N))
	}

	# we need to always use this sample as it could be different when removing a variable
	comp_with = stats::complete.cases(data,genes_current,env_current)
	fit_with = fit
	if (search=="genes") p_value = stats::coef(summary(fit_with)$fit_genes)[,4]
	if (search=="env") p_value = stats::coef(summary(fit_with)$fit_env)[,4]
	# Non-cross-validated models
	for (j in 1:elements_N){
		if (search=="genes"){
			comp_without = stats::complete.cases(data,genes_current[,-j,drop=FALSE],env_current)
			fit_without = LEGIT(data=data[comp_with,,drop=FALSE], genes=genes_current[comp_with,-j,drop=FALSE], env=env_current[comp_with,,drop=FALSE], formula=formula, start_genes=start_genes[-j], start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
		}
		if (search=="env"){
			comp_without = stats::complete.cases(data,genes_current,env_current[,-j,drop=FALSE])
			fit_without = LEGIT(data=data[comp_with,,drop=FALSE], genes=genes_current[comp_with,,drop=FALSE], env=env_current[comp_with,-j,drop=FALSE], formula=formula, start_genes=start_genes, start_env=start_env[-j], eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
		}
		good[j] = p_value[j] <= p_threshold

		if (exclude_worse_AIC){
			if (fit_with$true_model_parameters$AIC <= fit_without$true_model_parameters$AIC) good[j]= good[j] || TRUE
		}
		if (search_criterion=="AIC"){
			criterion_before[j] = fit_with$true_model_parameters$AIC
			criterion_after[j] = fit_without$true_model_parameters$AIC
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		}
		if (search_criterion=="AICc"){
			criterion_before[j] = fit_with$true_model_parameters$AICc
			criterion_after[j] = fit_without$true_model_parameters$AICc
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		}
		if (search_criterion=="BIC"){
			criterion_before[j] = fit_with$true_model_parameters$BIC
			criterion_after[j] = fit_without$true_model_parameters$BIC
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		}
		# Keep p-value, AIC and BIC if in interactive model
		if (interactive_mode){
			interactive$N_old[j] = sum(comp_with)
			interactive$N_new[j] = sum(comp_without)

			interactive$p_value[j] = round(p_value[j],6)

			interactive$AIC_old[j] = fit_with$true_model_parameters$AIC
			interactive$AIC_new[j] = fit_without$true_model_parameters$AIC

			interactive$AICc_old[j] = fit_with$true_model_parameters$AICc
			interactive$AICc_new[j] = fit_without$true_model_parameters$AICc

			interactive$BIC_old[j] = fit_with$true_model_parameters$BIC
			interactive$BIC_new[j] = fit_without$true_model_parameters$BIC
		}
	}
	if (sum(!good)==0){
		if (search=="genes" && print) cat("No gene removed\n")
		if (search=="env" && print) cat("No environment removed\n")
		return(NULL)
	}
	# Only do this if NOT in interactive mode, otherwise at the end of the algorithm, we show the the data.frame interactive and let the user choose
	if (!interactive_mode){
		if (search_criterion=="AIC" || search_criterion=="AICc" || search_criterion=="BIC"){
			if (min(criterion_diff,na.rm=TRUE) > 0){
				if (search=="genes" && print) cat("No gene removed\n")
				if (search=="env" && print) cat("No environment removed\n")
				return(NULL)
			}
			worst_var = which.min(criterion_diff)
			if (search=="genes"){
				if (print) cat(paste0("Removing gene: ",colnames(genes_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after = ",round(criterion_after[worst_var],5),")\n"))
				genes_dropped = cbind(genes_dropped, genes_current[,worst_var, drop=FALSE])
				genes_current = genes_current[,-worst_var, drop=FALSE]
			}
			if (search=="env"){
				if (print) cat(paste0("Removing environment: ",colnames(env_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after = ",round(criterion_after[worst_var],5),")\n"))
				env_dropped = cbind(env_dropped, env_current[,worst_var, drop=FALSE])
				env_current = env_current[,-worst_var, drop=FALSE]
			}
		}
	}
	# Cross-validated models
	if (search_criterion == "cv" || search_criterion == "cv_AUC" || search_criterion =="cv_Huber" || search_criterion =="cv_L1"){
		# Not looking at variables that labelled as good
		elements_N_cv = sum(!good)
		# Vector which says how much the criterion changed from removing the variable
		criterion_before = rep(NA, elements_N)
		criterion_after = rep(NA, elements_N)
		criterion_diff = rep(NA, elements_N)
		# Set seed
		if (!is.null(seed)) current_seed = seed
		else current_seed = NULL
		fit_cv_with = LEGIT_cv(data=data, genes=genes_current, env=env_current, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
		for (j in 1:elements_N){
			# Only do this if not labelled as good
			if (!good[j]){
				if (search=="genes") comp_without = stats::complete.cases(data,genes_current[,-j,drop=FALSE],env_current)
				if (search=="env") comp_without = stats::complete.cases(data,genes_current,env_current[,-j,drop=FALSE])
				if (search=="genes") fit_cv_without = LEGIT_cv(data=data[comp_with,,drop=FALSE], genes=genes_current[comp_with,-j,drop=FALSE], env=env_current[comp_with,,drop=FALSE], formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes[-j], start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
				if (search=="env") fit_cv_without = LEGIT_cv(data=data[comp_with,,drop=FALSE], genes=genes_current[comp_with,,drop=FALSE], env=env_current[comp_with,-j,drop=FALSE], formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env[-j], eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
				if (search_criterion=="cv"){
					criterion_before[j] = mean(fit_cv_with$R2_cv)
					criterion_after[j] = mean(fit_cv_without$R2_cv)
					criterion_diff[j] = criterion_after[j] - criterion_before[j]
				} 
				if (search_criterion=="cv_Huber"){
					criterion_before[j] = mean(fit_cv_with$Huber_cv)
					criterion_after[j] = mean(fit_cv_without$Huber_cv)
					criterion_diff[j] = criterion_after[j] - criterion_before[j]
				}
				if (search_criterion=="cv_L1"){
					criterion_before[j] = mean(fit_cv_with$L1_cv)
					criterion_after[j] = mean(fit_cv_without$L1_cv)
					criterion_diff[j] = criterion_after[j] - criterion_before[j]
				}
				if (search_criterion=="cv_AUC"){
					criterion_before[j] = mean(fit_cv_with$AUC)
					criterion_after[j] = mean(fit_cv_without$AUC)
					criterion_diff[j] = criterion_after[j] - criterion_before[j]
				}
				# Keep cross-validation R2 and AUC in interactive model
				if (interactive_mode){
					interactive$cv_R2_old[j] = mean(fit_cv_with$R2_cv)
					interactive$cv_R2_new[j] = mean(fit_cv_without$R2_cv)

					interactive$cv_Huber_old[j] = mean(fit_cv_with$Huber_cv)
					interactive$cv_Huber_new[j] = mean(fit_cv_without$Huber_cv)

					interactive$cv_L1_old[j] = mean(fit_cv_with$L1_cv)
					interactive$cv_L1_new[j] = mean(fit_cv_without$L1_cv)

					if(classification){
						interactive$cv_AUC_old[j] = mean(fit_cv_with$AUC)
						interactive$cv_AUC_new[j] = mean(fit_cv_without$AUC)
					}
				}
			}
		}
		# Only do this if NOT in interactive mode, otherwise at the end of the algorithm, we show the the data.frame interactive and let the user choose
		if (!interactive_mode){
			if ((max(criterion_diff,na.rm=TRUE) < 0 && !(search_criterion=="cv_Huber" || search_criterion=="cv_L1")) || (min(criterion_diff,na.rm=TRUE) > 0 && (search_criterion=="cv_Huber" || search_criterion=="cv_L1"))){
				if (search=="genes" && print) cat("No gene removed\n")
				if (search=="env" && print) cat("No environment removed\n")
				return(NULL)
			}
			if (search_criterion=="cv_Huber" || search_criterion=="cv_L1") worst_var = which.min(criterion_diff)
			else worst_var = which.max(criterion_diff)
			if (search=="genes"){
				if (print) cat(paste0("Removing gene: ",colnames(genes_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after = ",round(criterion_after[worst_var],5),")\n"))
				genes_dropped = cbind(genes_dropped, genes_current[,worst_var, drop=FALSE])
				genes_current = genes_current[,-worst_var, drop=FALSE]
			}
			if (search=="env"){
				if (print) cat(paste0("Removing environment: ",colnames(env_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after = ",round(criterion_after[worst_var],5),")\n"))
				env_dropped = cbind(env_dropped, env_current[,worst_var, drop=FALSE])
				env_current = env_current[,-worst_var, drop=FALSE]
			}
		}
	}
	if (interactive_mode){
		if (search_criterion=="cv") interactive_n = min(5,elements_N_cv)
		else interactive_n = min(5,elements_N)
		if (search_criterion=="AIC") neworder = order(interactive$AIC_new - interactive$AIC_old, decreasing=FALSE)
		if (search_criterion=="AICc") neworder = order(interactive$AICc_new - interactive$AICc_old, decreasing=FALSE)
		if (search_criterion=="BIC") neworder = order(interactive$BIC_new - interactive$BIC_old, decreasing=FALSE)
		if (search_criterion=="cv") neworder = order(interactive$cv_R2_new - interactive$cv_R2_old, decreasing=TRUE)
		if (search_criterion=="cv_Huber") neworder = order(interactive$cv_Huber_new - interactive$cv_Huber_old, decreasing=FALSE)
		if (search_criterion=="cv_L1") neworder = order(interactive$cv_L1_new - interactive$cv_L1_old, decreasing=FALSE)
		if (search_criterion=="cv_AUC") neworder = order(interactive$cv_AUC_new - interactive$cv_AUC_old, decreasing=TRUE)
		interactive = interactive[neworder[1:interactive_n],]
		rownames(interactive)=1:interactive_n
		interactive = format(interactive,scientific=FALSE)
		print(interactive)
		if (search=="genes") input_user = readline(prompt="Enter the index of the gene to be removed: ")
		if (search=="env") input_user = readline(prompt="Enter the index of the environment to be removed: ")
		if (sum(input_user == rownames(interactive))==0){
			if (search=="genes" && print) cat("No gene removed\n")
			if (search=="env" && print) cat("No environment removed\n")
			return(NULL)
		}
		worst_var = neworder[1:interactive_n][input_user == rownames(interactive)]
		if (search=="genes"){
			if (print) cat(paste0("Removing gene: ",colnames(genes_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after = ",round(criterion_after[worst_var],5),")\n"))
			genes_dropped = cbind(genes_dropped, genes_current[,worst_var, drop=FALSE])
			genes_current = genes_current[,-worst_var, drop=FALSE]
		}
		if (search=="env"){
			if (print) cat(paste0("Removing environment: ",colnames(env_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after = ",round(criterion_after[worst_var],5),")\n"))
			env_dropped = cbind(env_dropped, env_current[,worst_var, drop=FALSE])
			env_current = env_current[,-worst_var, drop=FALSE]
		}
	}
	# Updated model and coefficients
	if (search=="genes") start_genes=start_genes[-worst_var]
	if (search=="env") start_env=start_env[-worst_var]
	fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
	start_genes = stats::coef(fit$fit_genes)
	start_env = stats::coef(fit$fit_env)
	if (search=="genes") return(list(fit=fit, start_genes=start_genes,start_env=start_env,genes_current=genes_current,genes_dropped=genes_dropped))
	if (search=="env") return(list(fit=fit, start_genes=start_genes,start_env=start_env,env_current=env_current,env_dropped=env_dropped))
}

backward_step_IM = function(fit, data, formula, interactive_mode=FALSE, latent_var_current=NULL, latent_var_dropped=NULL, search=NULL, search_criterion="AIC", p_threshold = .20, exclude_worse_AIC=TRUE, max_steps = 100, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, classification=FALSE, start_latent_var=start_latent_var, eps=.01, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, print=TRUE){
	k = length(latent_var_current)
	# How much genes or env to add
	elements_N = NCOL(latent_var_current[[search]])
	if (elements_N <= 1){
		if (print) cat(paste0("No element from ", names(latent_var_current)[search]," was removed\n"))
		return(NULL)
	}
	# Vector which says which variables are "good" (worth keeping)
	good = rep(FALSE, elements_N)
	# Vector which says how much the criterion changed from excluding the variable
	criterion_before = rep(NA, elements_N)
	criterion_after = rep(NA, elements_N)
	criterion_diff = rep(NA, elements_N)
	# In interactive model, we must keep track of every AIC, BIC, p-value, Cross-validated R2 and AUC to show the user at every iteration
	if (interactive_mode){
		if (search_criterion == "cv" || search_criterion == "cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1") interactive = data.frame(variable=colnames(latent_var_current[[search]]), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N),cv_R2_old=rep(NA, elements_N),cv_R2_new=rep(NA, elements_N),cv_AUC_old=rep(NA, elements_N),cv_AUC_new=rep(NA, elements_N),cv_Huber_old=rep(NA, elements_N),cv_Huber_new=rep(NA, elements_N),cv_L1_old=rep(NA, elements_N),cv_L1_new=rep(NA, elements_N))
		else interactive = data.frame(variable=colnames(latent_var_current[[search]]), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),AICc_old=rep(NA, elements_N),AICc_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N))
	}	
	# we need to always use this sample as it could be different when removing a variable
	comp_with = stats::complete.cases(data,latent_var_current[[1]])
	if (k > 1) for (i in 2:k){
		comp_with = comp_with & stats::complete.cases(latent_var_current[[i]])
	}
	fit_with = fit
	p_value = stats::coef(summary(fit_with)[[search+1]])[,4]
	# Non-cross-validated models
	for (j in 1:elements_N){
		# Compelte dataset without variable
		if (search == 1) comp_without = stats::complete.cases(data,latent_var_current[[1]][,-j,drop=FALSE])
		else comp_without = stats::complete.cases(data,latent_var_current[[1]])
		if (k > 1) for (i in 2:k){
			if (search == i) comp_without = comp_without & stats::complete.cases(latent_var_current[[i]][,-j,drop=FALSE])
			else comp_without = comp_without & stats::complete.cases(latent_var_current[[i]])
		}
		latent_var_without = latent_var_current
		for (i in 1:k){
			if (i == search) latent_var_without[[i]] = latent_var_current[[i]][comp_with,-j, drop=FALSE]
			else latent_var_without[[i]] = latent_var_current[[i]][comp_with,, drop=FALSE]
		}
		start_latent_var_without = start_latent_var
		start_latent_var_without[[search]] = start_latent_var[[search]][-j]
		fit_without = IMLEGIT(data=data[comp_with,,drop=FALSE], latent_var=latent_var_without, formula=formula, start_latent_var=start_latent_var_without, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
		good[j] = p_value[j] <= p_threshold

		if (exclude_worse_AIC){
			if (fit_with$true_model_parameters$AIC <= fit_without$true_model_parameters$AIC) good[j]= good[j] || TRUE
		}
		if (search_criterion=="AIC"){
			criterion_before[j] = fit_with$true_model_parameters$AIC
			criterion_after[j] = fit_without$true_model_parameters$AIC
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		}
		if (search_criterion=="AICc"){
			criterion_before[j] = fit_with$true_model_parameters$AICc
			criterion_after[j] = fit_without$true_model_parameters$AICc
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		}
		if (search_criterion=="BIC"){
			criterion_before[j] = fit_with$true_model_parameters$BIC
			criterion_after[j] = fit_without$true_model_parameters$BIC
			criterion_diff[j] = criterion_after[j] - criterion_before[j]
		}
		# Keep p-value, AIC and BIC if in interactive model
		if (interactive_mode){
			interactive$N_old[j] = sum(comp_with)
			interactive$N_new[j] = sum(comp_without)

			interactive$p_value[j] = round(p_value[j],6)

			interactive$AIC_old[j] = fit_with$true_model_parameters$AIC
			interactive$AIC_new[j] = fit_without$true_model_parameters$AIC

			interactive$AICc_old[j] = fit_with$true_model_parameters$AICc
			interactive$AICc_new[j] = fit_without$true_model_parameters$AICc

			interactive$BIC_old[j] = fit_with$true_model_parameters$BIC
			interactive$BIC_new[j] = fit_without$true_model_parameters$BIC
		}
	}
	if (sum(!good)==0){
		if (print) cat(paste0("No element from ", names(latent_var_current)[search]," was removed\n"))
		return(NULL)
	}
	# Only do this if NOT in interactive mode, otherwise at the end of the algorithm, we show the the data.frame interactive and let the user choose
	if (!interactive_mode){
		if (search_criterion=="AIC" || search_criterion=="AICc" || search_criterion=="BIC"){
			if (min(criterion_diff,na.rm=TRUE) > 0){
				if (print) cat(paste0("No element from ", names(latent_var_current)[search]," was removed\n"))
				return(NULL)
			}
			worst_var = which.min(criterion_diff)
			if (print) cat(paste0("Removing element from ", names(latent_var_current)[search], ": ",colnames(latent_var_current[[search]])[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after = ",round(criterion_after[worst_var],5),")\n"))
			latent_var_dropped[[search]] = cbind(latent_var_dropped[[search]], latent_var_current[[search]][,worst_var, drop=FALSE])
			latent_var_current[[search]] = latent_var_current[[search]][,-worst_var, drop=FALSE]
		}
	}
	# Cross-validated models
	if (search_criterion == "cv" || search_criterion == "cv_AUC" || search_criterion =="cv_Huber" || search_criterion =="cv_L1"){
		# Not looking at variables that labelled as good
		elements_N_cv = sum(!good)
		# Vector which says how much the criterion changed from removing the variable
		criterion_before = rep(NA, elements_N)
		criterion_after = rep(NA, elements_N)
		criterion_diff = rep(NA, elements_N)
		# Set seed
		if (!is.null(seed)) current_seed = seed
		else current_seed = NULL
		fit_cv_with = IMLEGIT_cv(data=data, latent_var=latent_var_current, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
		for (j in 1:elements_N){
			# Only do this if not labelled as good
			if (!good[j]){
				# Compelte dataset without variable
				if (search == 1) comp_without = stats::complete.cases(data,latent_var_current[[1]][,-j,drop=FALSE])
				else comp_without = stats::complete.cases(data,latent_var_current[[1]])
				if (k > 1) for (i in 2:k){
					if (search == i) comp_without = comp_without & stats::complete.cases(latent_var_current[[i]][,-j,drop=FALSE])
					else comp_without = comp_without & stats::complete.cases(latent_var_current[[i]])
				}
				latent_var_without = latent_var_current
				for (i in 1:k){
					if (i == search) latent_var_without[[i]] = latent_var_current[[i]][comp_with,-j, drop=FALSE]
					else latent_var_without[[i]] = latent_var_current[[i]][comp_with,, drop=FALSE]
				}
				start_latent_var_without = start_latent_var
				start_latent_var_without[[search]] = start_latent_var[[search]][-j]

				fit_cv_without = IMLEGIT_cv(data=data[comp_with,,drop=FALSE], latent_var=latent_var_without, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var_without, eps=eps, maxiter=maxiter, family=family, seed=current_seed, ylim=ylim)
				if (search_criterion=="cv"){
					criterion_before[j] = mean(fit_cv_with$R2_cv)
					criterion_after[j] = mean(fit_cv_without$R2_cv)
					criterion_diff[j] = criterion_after[j] - criterion_before[j]
				} 
				if (search_criterion=="cv_Huber"){
					criterion_before[j] = mean(fit_cv_with$Huber_cv)
					criterion_after[j] = mean(fit_cv_without$Huber_cv)
					criterion_diff[j] = criterion_after[j] - criterion_before[j]
				}
				if (search_criterion=="cv_L1"){
					criterion_before[j] = mean(fit_cv_with$L1_cv)
					criterion_after[j] = mean(fit_cv_without$L1_cv)
					criterion_diff[j] = criterion_after[j] - criterion_before[j]
				}
				if (search_criterion=="cv_AUC"){
					criterion_before[j] = mean(fit_cv_with$AUC)
					criterion_after[j] = mean(fit_cv_without$AUC)
					criterion_diff[j] = criterion_after[j] - criterion_before[j]
				}
				# Keep cross-validation R2 and AUC in interactive model
				if (interactive_mode){
					interactive$cv_R2_old[j] = mean(fit_cv_with$R2_cv)
					interactive$cv_R2_new[j] = mean(fit_cv_without$R2_cv)

					interactive$cv_Huber_old[j] = mean(fit_cv_with$Huber_cv)
					interactive$cv_Huber_new[j] = mean(fit_cv_without$Huber_cv)

					interactive$cv_L1_old[j] = mean(fit_cv_with$L1_cv)
					interactive$cv_L1_new[j] = mean(fit_cv_without$L1_cv)

					if(classification){
						interactive$cv_AUC_old[j] = mean(fit_cv_with$AUC)
						interactive$cv_AUC_new[j] = mean(fit_cv_without$AUC)
					}
				}
			}
		}
		# Only do this if NOT in interactive mode, otherwise at the end of the algorithm, we show the the data.frame interactive and let the user choose
		if (!interactive_mode){
			if ((max(criterion_diff,na.rm=TRUE) < 0 && !(search_criterion=="cv_Huber" || search_criterion=="cv_L1")) || (min(criterion_diff,na.rm=TRUE) > 0 && (search_criterion=="cv_Huber" || search_criterion=="cv_L1"))){
				if (print) cat(paste0("No element from ", names(latent_var_current)[search]," was removed\n"))
				return(NULL)
			}
			if (search_criterion=="cv_Huber" || search_criterion=="cv_L1") worst_var = which.min(criterion_diff)
			else worst_var = which.max(criterion_diff)
			if (print) cat(paste0("Removing element from ", names(latent_var_current)[search], ": ",colnames(latent_var_current[[search]])[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after = ",round(criterion_after[worst_var],5),")\n"))
			latent_var_dropped[[search]] = cbind(latent_var_dropped[[search]], latent_var_current[[search]][,worst_var, drop=FALSE])
			latent_var_current[[search]] = latent_var_current[[search]][,-worst_var, drop=FALSE]
		}
	}
	if (interactive_mode){
		if (search_criterion=="cv") interactive_n = min(5,elements_N_cv)
		else interactive_n = min(5,elements_N)
		if (search_criterion=="AIC") neworder = order(interactive$AIC_new - interactive$AIC_old, decreasing=FALSE)
		if (search_criterion=="AICc") neworder = order(interactive$AICc_new - interactive$AICc_old, decreasing=FALSE)
		if (search_criterion=="BIC") neworder = order(interactive$BIC_new - interactive$BIC_old, decreasing=FALSE)
		if (search_criterion=="cv") neworder = order(interactive$cv_R2_new - interactive$cv_R2_old, decreasing=TRUE)
		if (search_criterion=="cv_Huber") neworder = order(interactive$cv_Huber_new - interactive$cv_Huber_old, decreasing=FALSE)
		if (search_criterion=="cv_L1") neworder = order(interactive$cv_L1_new - interactive$cv_L1_old, decreasing=FALSE)
		if (search_criterion=="cv_AUC") neworder = order(interactive$cv_AUC_new - interactive$cv_AUC_old, decreasing=TRUE)
		interactive = interactive[neworder[1:interactive_n],]
		rownames(interactive)=1:interactive_n
		interactive = format(interactive,scientific=FALSE)
		print(interactive)
		input_user = readline(prompt=paste0("Enter the index of the element from ", names(latent_var_current)[search], " to be removed: "))
		if (sum(input_user == rownames(interactive))==0){
			if (print) cat(paste0("No element from ", names(latent_var_current)[search]," was removed\n"))
			return(NULL)
		}
		worst_var = neworder[1:interactive_n][input_user == rownames(interactive)]
		if (print) cat(paste0("Removing element from ", names(latent_var_current)[search], ": ",colnames(latent_var_current[[search]])[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after = ",round(criterion_after[worst_var],5),")\n"))
		latent_var_dropped[[search]] = cbind(latent_var_dropped[[search]], latent_var_current[[search]][,worst_var, drop=FALSE])
		latent_var_current[[search]] = latent_var_current[[search]][,-worst_var, drop=FALSE]

	}
	# Updated model and coefficients
	start_latent_var[[search]] = start_latent_var[[search]][-worst_var]
	fit = IMLEGIT(data=data, latent_var=latent_var_current, formula=formula, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
	for (i in 1:k) start_latent_var[[i]] = stats::coef(fit$fit_latent_var[[i]])
	return(list(fit=fit, start_latent_var=start_latent_var,latent_var_current=latent_var_current,latent_var_dropped=latent_var_dropped))
}


stepwise_search = function(data, formula, interactive_mode=FALSE, genes_original=NULL, env_original=NULL, genes_extra=NULL, env_extra=NULL, search_type="bidirectional-forward", search="both", search_criterion="AIC", forward_exclude_p_bigger = .20, backward_exclude_p_smaller = .01, exclude_worse_AIC=TRUE, max_steps = 100, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, classification=FALSE, start_genes=NULL, start_env=NULL, eps=.01, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, print=TRUE, remove_miss=FALSE){
	if (forward_exclude_p_bigger > 1 || forward_exclude_p_bigger <= 0) stop("forward_exclude_p_bigger must be between 0 and 1 (Set to 1 to ignore p-values in forward step)")
	if (backward_exclude_p_smaller >= 1 || backward_exclude_p_smaller < 0) stop("backward_exclude_p_smaller must be between 0 and 1 (Set to 0 to ignore p-values in backward step)")
	if (search_criterion=="AIC") string_choice="lowest AIC"
	else if (search_criterion=="AICc") string_choice="lowest AICc"
	else if (search_criterion=="BIC") string_choice="lowest BIC"
	else if (search_criterion=="cv") string_choice="lowest cross-validation error"
	else if (search_criterion=="cv_Huber") string_choice="lowest cross-validation Huber error"
	else if (search_criterion=="cv_L1") string_choice="lowest cross-validation L1-norm error"
	else if (search_criterion=="cv_AUC") string_choice="biggest cross-validated area under the curve"
	else stop("Not a valid search_criterion, use: AIC, AICc, BIC, cv, cv_Huber, cv_L1 or cv_AUC")

	# Retaining only the needed variables from the dataset (need to set G and E variables for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	data$G=0
	data$E=0
	data = stats::model.frame(formula, data=data, na.action=na.pass)

	comp = stats::complete.cases(data, genes_original, env_original, genes_extra, env_extra)
	if (remove_miss){
		data = data[comp,, drop=FALSE]
		if (!is.null(genes_original)) genes_original = genes_original[comp,, drop=FALSE]
		if (!is.null(genes_extra)) genes_extra = genes_extra[comp,, drop=FALSE]
		if (!is.null(env_original)) env_original = env_original[comp,, drop=FALSE]
		if (!is.null(env_extra)) env_extra = env_extra[comp,, drop=FALSE]
		if (dim(data)[1] <= 0) stop("no valid observation without missing values")
	}
	else if (sum(comp)!=length(comp) && print) cat("Note: Missing data was found. This will increase computing time in forward steps \n and could lead to an incorrect subset of variable if the sample size change too much. \n You can use remove_miss=TRUE to force the removal of missing data. \n")

	if (interactive_mode && print) cat("<<~ Interative mode enabled ~>>\n")

	# If true, then we start with no genes or envand must find the best one first
	if ((is.null(genes_original) && search=="genes") || (is.null(env_original) && search=="env")) empty_start_dataset = TRUE
	else empty_start_dataset = FALSE
	if ((is.null(genes_original) || is.null(env_original)) && search=="both") stop("To do a stepwise search on both G and E, you must start with at least a single variable in G and in E. Please set genes_original and env_original.")

	# Setting up initial weighted scores
	if (is.null(genes_original) && search=="genes") start_genes = c()
	else if (is.null(start_genes)) start_genes = rep(1/dim(genes_original)[2],dim(genes_original)[2])

	if (is.null(env_original) && search=="env") start_env = c()
	else if (is.null(start_env)) start_env = rep(1/dim(env_original)[2],dim(env_original)[2])

	if (is.null(start_genes) && search=="both") start_genes = rep(1/dim(genes_original)[2],dim(genes_original)[2])
	if (is.null(start_env) && search=="both") start_env = rep(1/dim(env_original)[2],dim(env_original)[2])

	fit = NULL

	if (search_type == "forward"){
		if (print){
			if (forward_exclude_p_bigger < 1 && (exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Keeping only variables with p-values smaller than ", forward_exclude_p_bigger, " and which inclusion decrease the AIC\n"))
			else if (forward_exclude_p_bigger < 1 && !(exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Keeping only variables with p-values smaller than ", forward_exclude_p_bigger,"\n"))
			else if (exclude_worse_AIC || search_criterion=="AIC") cat(paste0("Keeping only variables which inclusion decrease the AIC\n"))
			if ((search_criterion=="cv" || search_criterion=="cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1") && (!exclude_worse_AIC || forward_exclude_p_bigger > .2)) cat("Note : We recommend using exclude_worse_AIC=TRUE and forward_exclude_p_bigger <= .20 to reduce the amount of cross-validations needed and to make the algorithm much faster.\n")
			if (search_criterion=="cv_AUC") classification=TRUE
			if (search=="genes") cat(paste0("Forward search of the genes to find the model with the ", string_choice,"\n"))
			if (search=="env") cat(paste0("Forward search of the environments to find the model with the ", string_choice,"\n"))
			if (search=="both") cat(paste0("Forward search of the genes and environments to find the model with the ", string_choice,"\n"))
		}
		genes_current = genes_original
		env_current = env_original
		if (!empty_start_dataset){
			# Original model
			fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
			start_genes = stats::coef(fit$fit_genes)
			start_env = stats::coef(fit$fit_env)
		}
		else{
			# If dataset is NULL, then make the current dataset be an empty data frame of same size as the dataset with extra variables 
			if (is.null(genes_original)) genes_current = genes_extra[,-c(1:NCOL(genes_extra))]
			if (is.null(env_original)) env_current = env_extra[,-c(1:NCOL(env_extra))]
		}
		genes_toadd = NULL
		env_toadd = NULL
		if (search=="genes" || search=="both") genes_toadd = genes_extra
		if (search=="env" || search=="both") env_toadd = env_extra
		if (print) cat("\n")

		if (search=="both"){
			G_conv = FALSE
			E_conv = FALSE
			# Starts looking at genes when search="both"
			search_current = "genes"
		}
		else search_current = search
		# Overall convergence (equal to G_conv && E_conv if search=both)
		conv = FALSE

		for (i in 1:max_steps){
			if (print) cat(paste0("[Iteration: ",i,"]\n"))
			results = forward_step(empty_start_dataset=empty_start_dataset, fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, genes_current=genes_current, env_current=env_current, genes_toadd=genes_toadd, env_toadd=env_toadd, search=search_current, search_criterion=search_criterion, p_threshold = forward_exclude_p_bigger, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
			if (is.null(results)){
				if (search=="both"){
					if (search_current == "genes"){
						search_current = "env"
						G_conv = TRUE
						if (E_conv) conv = TRUE
					}
					else{
						search_current = "genes"
						E_conv = TRUE
						if (G_conv) conv = TRUE
					}
				}
				else conv = TRUE
				if (conv) break
			}
			else{
				if (search_current == "genes") E_conv = FALSE
				if (search_current == "env") G_conv = FALSE
				# Resetting parameters based on iteration results
				empty_start_dataset = FALSE
				fit = results$fit
				start_genes=results$start_genes
				start_env=results$start_env
				if (search_current=="genes"){
					genes_current=results$genes_current
					genes_toadd=results$genes_toadd
				}
				if (search_current=="env"){
					env_current=results$env_current
					env_toadd=results$env_toadd
				}
			}
		}
	}
	if (search_type == "backward"){
		if (print){
			if (backward_exclude_p_smaller < 1 && (exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Dropping only variables with p-values bigger than ", backward_exclude_p_smaller, " and which removal decrease the AIC\n"))
			else if (backward_exclude_p_smaller < 1 && !(exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Dropping only variables with p-values bigger than ", backward_exclude_p_smaller,"\n"))
			else if (exclude_worse_AIC || search_criterion=="AIC") cat(paste0("Dropping only variables which removal decrease the AIC\n"))
			if ((search_criterion=="cv" || search_criterion=="cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1") && (!exclude_worse_AIC || backward_exclude_p_smaller < .01)) cat("Note : We recommend using exclude_worse_AIC=TRUE and backward_exclude_p_smaller >= .01 to reduce the amount of cross-validations needed and to make the algorithm much faster.\n")
			if (search_criterion=="cv_AUC") classification=TRUE
			if (search=="genes") cat(paste0("Backward search of the genes to find the model with the ", string_choice,"\n"))
			if (search=="env") cat(paste0("Backward search of the environments to find the model with the ", string_choice,"\n"))
			if (search=="both") cat(paste0("Backward search of the genes and environments to find the model with the ", string_choice,"\n"))
		}
		genes_current = genes_original
		env_current = env_original
		# Original model
		fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
		start_genes = stats::coef(fit$fit_genes)
		start_env = stats::coef(fit$fit_env)
		# Creating empty data frames of same size as the other data frames
		genes_dropped = genes_original[,-c(1:NCOL(genes_original))]
		env_dropped = env_original[,-c(1:NCOL(env_original))]
		if (print) cat("\n")

		if (search=="both"){
			G_conv = FALSE
			E_conv = FALSE
			# Starts looking at genes when search="both"
			search_current = "genes"
		}
		else search_current = search
		# Overall convergence (equal to G_conv && E_conv if search=both)
		conv = FALSE

		for (i in 1:max_steps){
			if (print) cat(paste0("[Iteration: ",i,"]\n"))
			results = backward_step(fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, genes_current=genes_current, genes_dropped=genes_dropped, env_current=env_current, env_dropped=env_dropped, search=search_current, search_criterion=search_criterion, p_threshold = backward_exclude_p_smaller, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
			if (is.null(results)){
				if (search=="both"){
					if (search_current == "genes"){
						search_current = "env"
						G_conv = TRUE
						if (E_conv) conv = TRUE
					}
					else{
						search_current = "genes"
						E_conv = TRUE
						if (G_conv) conv = TRUE
					}
				}
				else conv = TRUE
				if (conv) break
			}
			else{
				if (search_current == "genes") E_conv = FALSE
				if (search_current == "env") G_conv = FALSE
				# Resetting parameters based on iteration results
				fit = results$fit
				start_genes=results$start_genes
				start_env=results$start_env
				if (search_current=="genes"){
					genes_current=results$genes_current
					# Only used in bidirectionnal
					genes_dropped=results$genes_dropped
					if (NCOL(genes_current)==1){
						if(search=="both"){
							search_current = "env"
							G_conv = TRUE
							if (E_conv) conv = TRUE
						}
						else conv = TRUE
						if (conv) break
					}
				}
				if (search_current=="env"){
					env_current=results$env_current
					# Only used in bidirectionnal
					env_dropped=results$env_dropped
					if (NCOL(env_current)==1){
						if(search=="both"){
							search_current = "genes"
							E_conv = TRUE
							if (G_conv) conv = TRUE
						}
						else conv = TRUE
						if (conv) break
					}
				}
			}
		}
	}
	if (search_type == "bidirectional-forward" || search_type == "bidirectional-backward"){
		if (print){
			if (forward_exclude_p_bigger < 1 && (exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Keeping only variables with p-values smaller than ", forward_exclude_p_bigger, " and which inclusion decrease the AIC\n"))
			else if (forward_exclude_p_bigger < 1 && !(exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Keeping only variables with p-values smaller than ", forward_exclude_p_bigger,"\n"))
			else if (exclude_worse_AIC || search_criterion=="AIC") cat(paste0("Keeping only variables which inclusion decrease the AIC\n"))
			if (backward_exclude_p_smaller < 1 && (exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Dropping only variables with p-values bigger than ", backward_exclude_p_smaller, " and which removal decrease the AIC\n"))
			else if (backward_exclude_p_smaller < 1 && !(exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Dropping only variables with p-values bigger than ", backward_exclude_p_smaller,"\n"))
			else if (exclude_worse_AIC || search_criterion=="AIC") cat(paste0("Dropping only variables which removal decrease the AIC\n"))
			if ((search_criterion=="cv" || search_criterion=="cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1") && (!exclude_worse_AIC || forward_exclude_p_bigger > .2)) cat("Note : We recommend using exclude_worse_AIC=TRUE and forward_exclude_p_bigger <= .20 to reduce the amount of cross-validations needed and to make the algorithm much faster.\n")
			if (search_criterion=="cv_AUC") classification=TRUE
			if (search=="genes") cat(paste0("Bidirectional search of the genes to find the model with the ", string_choice,"\n"))
			if (search=="env") cat(paste0("Bidirectional search of the environments to find the model with the ", string_choice,"\n"))
			if (search=="both") cat(paste0("Bidirectional search of the genes and environments to find the model with the ", string_choice,"\n"))
		}
		if (search_type == "bidirectional-forward"){
			genes_current = genes_original
			env_current = env_original
			if (!empty_start_dataset){
				# Original model
				fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
				start_genes = stats::coef(fit$fit_genes)
				start_env = stats::coef(fit$fit_env)
			}
			else{
				# If dataset is NULL, then make the current dataset be an empty data frame of same size as the dataset with extra variables 
				if (is.null(genes_original)) genes_current = genes_extra[,-c(1:NCOL(genes_extra))]
				if (is.null(env_original)) env_current = env_extra[,-c(1:NCOL(env_extra))]
			}
			genes_toadd_ordrop = NULL
			env_toadd_ordrop = NULL
			if (search=="genes" || search=="both") genes_toadd_ordrop = genes_extra
			if (search=="env" || search=="both") env_toadd_ordrop = env_extra
			direction="forward"
			forward_failed=FALSE
			# Count as failed because we can't backward if forward doesn't work
			backward_failed=TRUE
		}
		if (search_type == "bidirectional-backward"){
			genes_current = genes_original
			env_current = env_original
			# Original model
			fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
			start_genes = stats::coef(fit$fit_genes)
			start_env = stats::coef(fit$fit_env)
			# Creating empty data frames of same size as the other data frames
			genes_toadd_ordrop = genes_original[,-c(1:NCOL(genes_original))]
			env_toadd_ordrop = env_original[,-c(1:NCOL(env_original))]
			direction="backward"
			# Count as failed because we can't backward if forward doesn't work
			forward_failed=TRUE
			backward_failed=FALSE
		}
		if (print) cat("\n")

		if (search=="both"){
			G_conv = FALSE
			E_conv = FALSE
			# Starts looking at genes when search="both"
			search_current = "genes"
		}
		else search_current = search
		# Overall convergence (equal to G_conv && E_conv if search=both)
		conv = FALSE

		for (i in 1:max_steps){
			if (print) cat(paste0("[Iteration: ",i,"]\n"))
			if (direction=="forward"){
				results = forward_step(empty_start_dataset=empty_start_dataset, fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, genes_current=genes_current, env_current=env_current, genes_toadd=genes_toadd_ordrop, env_toadd=env_toadd_ordrop, search=search_current, search_criterion=search_criterion, p_threshold = forward_exclude_p_bigger, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
				if (is.null(results)) forward_failed=TRUE
				else{
					if (search_current == "genes") E_conv = FALSE
					if (search_current == "env") G_conv = FALSE
					forward_failed=FALSE
					# Resetting parameters based on iteration results
					empty_start_dataset = FALSE
					fit = results$fit
					start_genes=results$start_genes
					start_env=results$start_env
					if (search_current=="genes"){
						genes_current=results$genes_current
						genes_toadd_ordrop=results$genes_toadd
						if (NCOL(genes_toadd_ordrop)==0){
							# Count as a fail
							forward_failed=TRUE
						}
					}
					if (search_current=="env"){
						env_current=results$env_current
						env_toadd_ordrop=results$env_toadd
						if (NCOL(env_toadd_ordrop)==0){
							# Count as a fail
							forward_failed=TRUE
						}
					}
				}
				direction="backward"
			}
			if (forward_failed && backward_failed){
				if (search=="both"){
					backward_failed = FALSE
					forward_failed = FALSE
					if (search_current == "genes"){
						search_current = "env"
						# Count as a fail
						if (NCOL(env_current)==1) backward_failed=TRUE
						if (NCOL(env_toadd_ordrop)==0) forward_failed=TRUE
						G_conv = TRUE
						if (E_conv) conv = TRUE
					}
					else{
						search_current = "genes"
						# Count as a fail
						if (NCOL(genes_current)==1) backward_failed=TRUE
						if (NCOL(genes_toadd_ordrop)==0) forward_failed=TRUE
						E_conv = TRUE
						if (G_conv) conv = TRUE
					}
					if (search_type == "bidirectional-forward") direction="forward"
					if (search_type == "bidirectional-backward") direction="backward"
				}
				else conv = TRUE
				if (conv) break
			}
			if (direction=="backward"){
				# Can't backward if only one remaining variable
				if (!((NCOL(genes_current)<=1 && search_current=="genes") || (NCOL(env_current)<=1 && search_current=="env"))){
					results = backward_step(fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, genes_current=genes_current, genes_dropped=genes_toadd_ordrop, env_current=env_current, env_dropped=env_toadd_ordrop, search=search_current, search_criterion=search_criterion, p_threshold = backward_exclude_p_smaller, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
					if (is.null(results)) backward_failed=TRUE
					else{
						if (search_current == "genes") E_conv = FALSE
						if (search_current == "env") G_conv = FALSE
						backward_failed=FALSE
						# Resetting parameters based on iteration results
						fit = results$fit
						start_genes=results$start_genes
						start_env=results$start_env
						one_remain=FALSE
						if (search_current=="genes"){
							genes_current=results$genes_current
							if (NCOL(genes_current)==1){
								# Count as a fail
								backward_failed=TRUE
							}
							genes_toadd_ordrop=results$genes_dropped
						}
						if (search_current=="env"){
							env_current=results$env_current
							if (NCOL(env_current)==1){
								# Count as a fail
								backward_failed=TRUE
							}
							env_toadd_ordrop=results$env_dropped
						}
					}
				}
				else{
					backward_failed=TRUE
					if (search_current=="genes" && print) cat("No gene removed\n")
					if (search_current=="env" && print) cat("No environment removed\n")
				}
				direction="forward"
			}
			if (forward_failed && backward_failed){
				if (search=="both"){
					backward_failed = FALSE
					forward_failed = FALSE
					if (search_current == "genes"){
						search_current = "env"
						# Count as a fail
						if (NCOL(env_current)==1) backward_failed=TRUE
						if (NCOL(env_toadd_ordrop)==0) forward_failed=TRUE
						G_conv = TRUE
						if (E_conv) conv = TRUE
					}
					else{
						search_current = "genes"
						# Count as a fail
						if (NCOL(genes_current)==1) backward_failed=TRUE
						if (NCOL(genes_toadd_ordrop)==0) forward_failed=TRUE
						E_conv = TRUE
						if (G_conv) conv = TRUE
					}
					if (search_type == "bidirectional-forward") direction="forward"
					if (search_type == "bidirectional-backward") direction="backward"
				}
				else conv = TRUE
				if (conv) break
			}
		}
	}
	if (i >= max_steps) warning("Stepwise search did not reach convergence in max_steps steps. Try increasing max_steps.")
	return(fit)
}


stepwise_search_IM = function(data, formula, interactive_mode=FALSE, latent_var_original=NULL, latent_var_extra=NULL, search_type="bidirectional-forward", search=0, search_criterion="AIC", forward_exclude_p_bigger = .20, backward_exclude_p_smaller = .01, exclude_worse_AIC=TRUE, max_steps = 100, cv_iter=5, cv_folds=10, folds=NULL, Huber_p=1.345, classification=FALSE, start_latent_var=NULL, eps=.01, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, print=TRUE, remove_miss=FALSE){
	k = max(length(latent_var_original),length(latent_var_extra))
	if (forward_exclude_p_bigger > 1 || forward_exclude_p_bigger <= 0) stop("forward_exclude_p_bigger must be between 0 and 1 (Set to 1 to ignore p-values in forward step)")
	if (backward_exclude_p_smaller >= 1 || backward_exclude_p_smaller < 0) stop("backward_exclude_p_smaller must be between 0 and 1 (Set to 0 to ignore p-values in backward step)")
	if (search_criterion=="AIC") string_choice="lowest AIC"
	else if (search_criterion=="AICc") string_choice="lowest AICc"
	else if (search_criterion=="BIC") string_choice="lowest BIC"
	else if (search_criterion=="cv") string_choice="lowest cross-validation error"
	else if (search_criterion=="cv_Huber") string_choice="lowest cross-validation Huber error"
	else if (search_criterion=="cv_L1") string_choice="lowest cross-validation L1-norm error"
	else if (search_criterion=="cv_AUC") string_choice="biggest cross-validated area under the curve"
	else stop("Not a valid search_criterion, use: AIC, AICc, BIC, cv, cv_Huber, cv_L1 or cv_AUC")

	# Retaining only the needed variables from the dataset (need to set G and E variables for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	for (i in 1:k) data[,names(latent_var_original)[i]] = 0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	#  Check for empty latent variable (Note: Probably shouldn't be named empty_start_dataset but empty_start_latent_var althought that would be even more confusing considering start_latent_var)
	if (is.null(latent_var_original)) stop("latent_var_original cannot be null. However, if search=i, then you could set latent_var_original[[i]]=NULL.")
	# Make it to have NULL elements but not be NULL
	comp = rep(TRUE, NROW(data))
	for (i in 1:k){
		if (!is.null(latent_var_original[[i]])) comp = comp & stats::complete.cases(data, latent_var_original[[i]])
		if (!is.null(latent_var_extra) && !is.null(latent_var_extra[[i]])) comp = comp & stats::complete.cases(data, latent_var_extra[[i]])
	}
	if (remove_miss){
		data = data[comp,, drop=FALSE]
		for (i in 1:k){
			if (!is.null(latent_var_original[[i]])) latent_var_original[[i]] = latent_var_original[[i]][comp,, drop=FALSE]
			if (!is.null(latent_var_extra)){
				if (!is.null(latent_var_extra[[i]])) latent_var_extra[[i]] = latent_var_extra[[i]][comp,, drop=FALSE]
			}
		}
		if (dim(data)[1] <= 0) stop("no valid observation without missing values")
	}
	else if (sum(comp)!=length(comp) && print) cat("Note: Missing data was found. This will increase computing time in forward steps \n and could lead to an incorrect subset of variable if the sample size change too much. \n You can use remove_miss=TRUE to force the removal of missing data. \n")

	if (interactive_mode && print) cat("<<~ Interative mode enabled ~>>\n")

	# Making all elements NULL instead
	if (is.null(latent_var_extra)) latent_var_extra = vector("list", k)

	empty_start_dataset = FALSE
	for (i in 1:k){
		if (is.null(latent_var_original[[i]]) || NCOL(latent_var_original[[i]])==0){
			if (i != search) stop("If search=i, you can set latent_var_original[[i]]=NULL but not latent_var_original[[j]]=NULL. When search=0, you must start with at least an element in every latent variable.")
			else if (search_type=="bidirectional-backward" || search_type=="backward") stop ("If search=i, you can normally set latent_var_original[[i]]=NULL but not with backward search because you must start from somewhere before dropping variables.")
			else empty_start_dataset = TRUE
			# Make an empty dataframe instead of NULL
			latent_var_original[[i]] = latent_var_extra[[i]][,-c(1:NCOL(latent_var_extra[[i]])), drop=FALSE]
		}
		if (is.null(latent_var_extra[[i]]) || NCOL(latent_var_extra[[i]])==0){
			# Make an empty dataframe instead of NULL
			latent_var_extra[[i]] = latent_var_original[[i]][,-c(1:NCOL(latent_var_original[[i]])), drop=FALSE]
		}
	}

	# Setting up initial start
	if (is.null(start_latent_var)){
		start_latent_var = vector("list", k)
		for (i in 1:k){
			if (empty_start_dataset && search==i) start_latent_var[[i]]=c()
			else start_latent_var[[i]] = rep(1/NCOL(latent_var_original[[i]]),NCOL(latent_var_original[[i]]))
		}
	}
	fit = NULL

	if (search_type == "forward"){
		if (print){
			if (forward_exclude_p_bigger < 1 && (exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Keeping only variables with p-values smaller than ", forward_exclude_p_bigger, " and which inclusion decrease the AIC\n"))
			else if (forward_exclude_p_bigger < 1 && !(exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Keeping only variables with p-values smaller than ", forward_exclude_p_bigger,"\n"))
			else if (exclude_worse_AIC || search_criterion=="AIC") cat(paste0("Keeping only variables which inclusion decrease the AIC\n"))
			if ((search_criterion=="cv" || search_criterion=="cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1") && (!exclude_worse_AIC || forward_exclude_p_bigger > .2)) cat("Note : We recommend using exclude_worse_AIC=TRUE and forward_exclude_p_bigger <= .20 to reduce the amount of cross-validations needed and to make the algorithm much faster.\n")
			if (search_criterion=="cv_AUC") classification=TRUE
			if (search==0) cat(paste0("Forward search of the elements from all latent variables to find the model with the ", string_choice,"\n"))
			else cat(paste0("Forward search of the elements from ", names(latent_var_original)[search]," to find the model with the ", string_choice,"\n"))
		}
		latent_var_current = latent_var_original
		if (!empty_start_dataset){
			# Original model
			fit = IMLEGIT(data=data, latent_var=latent_var_current, formula=formula, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
			for (i in 1:k) start_latent_var[[i]] = stats::coef(fit$fit_latent_var[[i]])
		}
		else{
			# If dataset is NULL, then make the current dataset be an empty data frame of same size as the dataset with extra variables 
			if (is.null(latent_var_original[[search]])) latent_var_original[[search]] = latent_var_extra[[search]][,-c(1:NCOL(latent_var_extra[[search]]))]
		}
		latent_var_toadd = latent_var_extra
		if (print) cat("\n")

		if (search == 0){
			latent_var_conv = rep(FALSE, k)
			search_current = 1
		}
		else search_current = search
		# Overall convergence
		conv = FALSE

		for (i in 1:max_steps){
			if (print) cat(paste0("[Iteration: ",i,"]\n"))
			results = forward_step_IM(empty_start_dataset=empty_start_dataset, fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, latent_var_current=latent_var_current, latent_var_toadd=latent_var_toadd, search=search_current, search_criterion=search_criterion, p_threshold = forward_exclude_p_bigger, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
			if (is.null(results)){
				if (search == 0){
					latent_var_conv[search_current] = TRUE
					if (sum(latent_var_conv)==k) conv = TRUE
					if (search_current != k) search_current = search_current + 1
					else search_current = 1
				}
				else conv = TRUE
				if (conv) break
			}
			else{
				# If this one didn't converge then all others have to be retried again since things can change
				if (search == 0) latent_var_conv = rep(FALSE, k)
				# Resetting parameters based on iteration results
				empty_start_dataset = FALSE
				fit = results$fit
				for (i in 1:k) start_latent_var[[i]] = stats::coef(fit$fit_latent_var[[i]])
				latent_var_current=results$latent_var_current
				latent_var_toadd=results$latent_var_toadd
			}
		}
	}
	if (search_type == "backward"){
		if (print){
			if (backward_exclude_p_smaller < 1 && (exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Dropping only variables with p-values bigger than ", backward_exclude_p_smaller, " and which removal decrease the AIC\n"))
			else if (backward_exclude_p_smaller < 1 && !(exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Dropping only variables with p-values bigger than ", backward_exclude_p_smaller,"\n"))
			else if (exclude_worse_AIC || search_criterion=="AIC") cat(paste0("Dropping only variables which removal decrease the AIC\n"))
			if ((search_criterion=="cv" || search_criterion=="cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1") && (!exclude_worse_AIC || backward_exclude_p_smaller < .01)) cat("Note : We recommend using exclude_worse_AIC=TRUE and backward_exclude_p_smaller >= .01 to reduce the amount of cross-validations needed and to make the algorithm much faster.\n")
			if (search_criterion=="cv_AUC") classification=TRUE
			if (search==0) cat(paste0("Backward search of the elements from all latent variables to find the model with the ", string_choice,"\n"))
			else cat(paste0("Backward search of the elements from ", names(latent_var_original)[search]," to find the model with the ", string_choice,"\n"))
		}
		latent_var_current = latent_var_original
		# Original model
		fit = IMLEGIT(data=data, latent_var=latent_var_current, formula=formula, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
		for (i in 1:k) start_latent_var[[i]] = stats::coef(fit$fit_latent_var[[i]])
		# Creating empty data frames of same size as the other data frames
		latent_var_dropped = latent_var_original
		for (i in 1:k) latent_var_dropped[[i]] = latent_var_original[[i]][,-c(1:NCOL(latent_var_original[[i]]))]

		if (print) cat("\n")

		if (search == 0){
			latent_var_conv = rep(FALSE, k)
			search_current = 1
		}
		else search_current = search
		# Overall convergence
		conv = FALSE

		for (i in 1:max_steps){
			if (print) cat(paste0("[Iteration: ",i,"]\n"))
			results = backward_step_IM(fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, latent_var_current=latent_var_current, latent_var_dropped=latent_var_dropped, search=search_current, search_criterion=search_criterion, p_threshold = backward_exclude_p_smaller, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
			if (is.null(results)){
				if (search==0){
					latent_var_conv[search_current] = TRUE
					if (sum(latent_var_conv)==k) conv = TRUE
					if (search_current != k) search_current = search_current + 1
					else search_current = 1
				}
				else conv = TRUE
				if (conv) break
			}
			else{
				# If this one didn't converge then all others have to be retried again since things can change
				if (search == 0) latent_var_conv = rep(FALSE, k)
				# Resetting parameters based on iteration results
				fit = results$fit
				for (i in 1:k) start_latent_var[[i]] = stats::coef(fit$fit_latent_var[[i]])
				latent_var_current=results$latent_var_current
				latent_var_dropped=results$latent_var_dropped
				if (NCOL(latent_var_current[[search_current]])==1){
				if (search == 0){
					latent_var_conv[search_current] = TRUE
					if (sum(latent_var_conv) == k) conv = TRUE	
					if (search_current != k) search_current = search_current + 1
					else search_current = 1
					}
				else conv = TRUE
				if (conv) break
			}
			}
		}
	}
	if (search_type == "bidirectional-forward" || search_type == "bidirectional-backward"){
		if (print){
			if (forward_exclude_p_bigger < 1 && (exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Keeping only variables with p-values smaller than ", forward_exclude_p_bigger, " and which inclusion decrease the AIC\n"))
			else if (forward_exclude_p_bigger < 1 && !(exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Keeping only variables with p-values smaller than ", forward_exclude_p_bigger,"\n"))
			else if (exclude_worse_AIC || search_criterion=="AIC") cat(paste0("Keeping only variables which inclusion decrease the AIC\n"))
			if (backward_exclude_p_smaller < 1 && (exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Dropping only variables with p-values bigger than ", backward_exclude_p_smaller, " and which removal decrease the AIC\n"))
			else if (backward_exclude_p_smaller < 1 && !(exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Dropping only variables with p-values bigger than ", backward_exclude_p_smaller,"\n"))
			else if (exclude_worse_AIC || search_criterion=="AIC") cat(paste0("Dropping only variables which removal decrease the AIC\n"))
			if ((search_criterion=="cv" || search_criterion=="cv_AUC" || search_criterion=="cv_Huber" || search_criterion=="cv_L1") && (!exclude_worse_AIC || forward_exclude_p_bigger > .2)) cat("Note : We recommend using exclude_worse_AIC=TRUE and forward_exclude_p_bigger <= .20 to reduce the amount of cross-validations needed and to make the algorithm much faster.\n")
			if (search_criterion=="cv_AUC") classification=TRUE
			if (search==0) cat(paste0("Bidirectional search of the elements from all latent variables to find the model with the ", string_choice,"\n"))
			else cat(paste0("Bidirectional search of the elements from ", names(latent_var_original)[search]," to find the model with the ", string_choice,"\n"))

		}
		if (search_type == "bidirectional-forward"){
			latent_var_current = latent_var_original
			if (!empty_start_dataset){
				# Original model
				fit = IMLEGIT(data=data, latent_var=latent_var_current, formula=formula, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
				for (i in 1:k) start_latent_var[[i]] = stats::coef(fit$fit_latent_var[[i]])
			}
			else{
				# If dataset is NULL, then make the current dataset be an empty data frame of same size as the dataset with extra variables 
				if (is.null(latent_var_original[[search]])) latent_var_original[[search]] = latent_var_extra[[search]][,-c(1:NCOL(latent_var_extra[[search]]))]
			}
			latent_var_toadd_ordrop = latent_var_extra

			direction="forward"
			forward_failed=FALSE
			# Count as failed because we can't backward if forward doesn't work
			backward_failed=TRUE
		}
		if (search_type == "bidirectional-backward"){
			latent_var_current = latent_var_original
			# Original model
			fit = IMLEGIT(data=data, latent_var=latent_var_current, formula=formula, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
			for (i in 1:k) start_latent_var[[i]] = stats::coef(fit$fit_latent_var[[i]])
			# Creating empty data frames of same size as the other data frames
			latent_var_toadd_ordrop = latent_var_original
			for (i in 1:k) latent_var_toadd_ordrop[[i]] = latent_var_original[[i]][,-c(1:NCOL(latent_var_original[[i]]))]
			direction="backward"
			# Count as failed because we can't backward if forward doesn't work
			forward_failed=TRUE
			backward_failed=FALSE
		}
		if (print) cat("\n")

		if (search == 0){
			latent_var_conv = rep(FALSE, k)
			search_current = 1
		}
		else search_current = search
		# Overall convergence
		conv = FALSE

		for (i in 1:max_steps){
			if (print) cat(paste0("[Iteration: ",i,"]\n"))
			if (direction=="forward"){
				results = forward_step_IM(empty_start_dataset=empty_start_dataset, fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, latent_var_current=latent_var_current, latent_var_toadd=latent_var_toadd_ordrop, search=search_current, search_criterion=search_criterion, p_threshold = forward_exclude_p_bigger, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
				if (is.null(results)) forward_failed=TRUE
				else{
					if (search == 0) latent_var_conv = rep(FALSE, k)
					forward_failed=FALSE
					# Resetting parameters based on iteration results
					empty_start_dataset = FALSE
					fit = results$fit
					for (i in 1:k) start_latent_var[[i]] = stats::coef(fit$fit_latent_var[[i]])
					latent_var_current=results$latent_var_current
					latent_var_toadd_ordrop=results$latent_var_toadd
					if (NCOL(latent_var_toadd_ordrop[[search_current]])==0){
						# Count as a fail
						forward_failed=TRUE
					}
				}
				direction="backward"
			}
			if (forward_failed && backward_failed){
				if (search == 0){
					backward_failed = FALSE
					forward_failed = FALSE
					latent_var_conv[search_current] = TRUE
					if (sum(latent_var_conv)==k) conv = TRUE
					if (search_current != k) search_current = search_current + 1
					else search_current = 1
					# Count as a fail
					if (NCOL(latent_var_current)==1) backward_failed=TRUE
					if (NCOL(latent_var_toadd_ordrop)==0) forward_failed=TRUE
					if (search_type == "bidirectional-forward") direction="forward"
					if (search_type == "bidirectional-backward") direction="backward"
				}
				else conv = TRUE
				if (conv) break
			}
			if (direction=="backward"){
				# Can't backward if only one remaining variable
				if (!(NCOL(latent_var_current[[search_current]])<=1)){
					results = backward_step_IM(fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, latent_var_current=latent_var_current, latent_var_dropped=latent_var_toadd_ordrop, search=search_current, search_criterion=search_criterion, p_threshold = backward_exclude_p_smaller, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
					if (is.null(results)) backward_failed=TRUE
					else{
						# If this one didn't converge then all others have to be retried again since things can change
						if (search == 0) latent_var_conv = rep(FALSE, k)
						backward_failed=FALSE
						# Resetting parameters based on iteration results
						fit = results$fit
						for (i in 1:k) start_latent_var[[i]] = stats::coef(fit$fit_latent_var[[i]])
						one_remain=FALSE
						latent_var_current=results$latent_var_current
						latent_var_toadd_ordrop=results$latent_var_dropped
						if (NCOL(latent_var_current[[search_current]])==1){
							# Count as a fail
							backward_failed=TRUE
						}
					}
				}
				else{
					backward_failed=TRUE
					if (print) cat(paste0("No element from ", names(latent_var_original)[search_current]," was removed\n"))
				}
				direction="forward"
			}
			if (forward_failed && backward_failed){
				if (search==0){
					backward_failed = FALSE
					forward_failed = FALSE

					if (search == 0) latent_var_conv[search_current] = TRUE
					if (sum(latent_var_conv) == k) conv = TRUE	
					if (search_current != k) search_current = search_current + 1
					else search_current = 1
					# Count as a fail
					if (NCOL(latent_var_current[[search_current]])==1) backward_failed=TRUE
					if (NCOL(latent_var_toadd_ordrop[[search_current]])==0) forward_failed=TRUE

					if (search_type == "bidirectional-forward") direction="forward"
					if (search_type == "bidirectional-backward") direction="backward"
				}
				else conv = TRUE
				if (conv) break
			}
		}
	}
	if (i >= max_steps) warning("Stepwise search did not reach convergence in max_steps steps. Try increasing max_steps.")
	return(fit)
}


bootstrap_var_select = function(data, formula, boot_iter=1000, boot_size=NULL, boot_group=NULL, latent_var_original=NULL, latent_var_extra=NULL, search_type="bidirectional-forward", search=0, search_criterion="AIC", forward_exclude_p_bigger = .20, backward_exclude_p_smaller = .01, exclude_worse_AIC=TRUE, max_steps = 100, start_latent_var=NULL, eps=.01, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, progress=TRUE, n_cluster = 1, best_subsets=5){
	k = max(length(latent_var_original),length(latent_var_extra))
	## Removing missing data and checks
	# Retaining only the needed variables from the dataset (need to set G and E variables for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	for (i in 1:k) data[,names(latent_var_original)[i]] = 0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	#  Check for empty latent variable (Note: Probably shouldn't be named empty_start_dataset but empty_start_latent_var althought that would be even more confusing considering start_latent_var)
	if (is.null(latent_var_original)) stop("latent_var_original cannot be null. However, if search=i, then you could set latent_var_original[[i]]=NULL.")
	# Make it to have NULL elements but not be NULL
	comp = rep(TRUE, NROW(data))
	for (i in 1:k){
		if (!is.null(latent_var_original[[i]])) comp = comp & stats::complete.cases(data, latent_var_original[[i]])
		if (!is.null(latent_var_extra) && !is.null(latent_var_extra[[i]])) comp = comp & stats::complete.cases(data, latent_var_extra[[i]])
	}
	if (!is.null(boot_group)) boot_group = boot_group[comp]
	data = data[comp,, drop=FALSE]
	for (i in 1:k){
		if (!is.null(latent_var_original[[i]])) latent_var_original[[i]] = latent_var_original[[i]][comp,, drop=FALSE]
		if (!is.null(latent_var_extra)){
			if (!is.null(latent_var_extra[[i]])) latent_var_extra[[i]] = latent_var_extra[[i]][comp,, drop=FALSE]
		}
	}
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	# Making all elements NULL instead
	if (is.null(latent_var_extra)) latent_var_extra = vector("list", k)

	## Bootstrapping
	if (is.null(boot_size)) boot_size = NROW(data)
	# Create list of variable counts
	var_select = vector("list", k)
	for (i in 1:k){
		if (is.null(latent_var_extra[[i]])) latent_var_all = latent_var_original[[i]]
		else if (is.null(latent_var_original[[i]])) latent_var_all = latent_var_extra[[i]]
		else latent_var_all = cbind(latent_var_original[[i]],latent_var_extra[[i]])
		var_select[[i]] = rep(0, NCOL(latent_var_all))
		names(var_select[[i]]) = colnames(latent_var_all)
	}
	# Number of times a subset of variable has appeared
	subset_count = c()

	# Setting up parallel
	cl <- snow::makeCluster(n_cluster)
	doSNOW::registerDoSNOW(cl)
	if (progress){
		if (runif(1) < .50) char="(^._.^) "
		else char="(=^o^=) "
		pb = utils::txtProgressBar(max = boot_iter, style = 3, char=char)
		progress <- function(n) utils::setTxtProgressBar(pb, n)
		opts <- list(progress = progress)
	}
	else opts <- list()
	# Need to use this "with(c(),CODEHERE)" to prevent R check from returning a "no visible binding for global variable"
	with(c(),{
		results <- foreach::foreach(b = 1:boot_iter, .options.snow = opts) %dopar% {
			if (!is.null(seed)) set.seed(seed+b)
			if (is.null(boot_group)) boot = sample(1:NROW(data),size=boot_size,replace=TRUE)
			else{
				boot = c()
				N = NROW(data)
				repeat{
					id_sampled = sample(unique(boot_group),size=1)
					toadd = which(boot_group == id_sampled)
					current_size = length(boot)
					if (abs(N-current_size) > abs(N - current_size - length(toadd))) boot = c(boot, toadd)
					else break
				}
			}
			data_b = data[boot,, drop=FALSE]
			latent_var_original_b = latent_var_original
			latent_var_extra_b = latent_var_extra
			for (i in 1:k){
				if (!is.null(latent_var_original[[i]])) latent_var_original_b[[i]] = latent_var_original_b[[i]][boot,, drop=FALSE]
				if (!is.null(latent_var_extra)){
					if (!is.null(latent_var_extra[[i]])) latent_var_extra_b[[i]] = latent_var_extra_b[[i]][boot,, drop=FALSE]
				}
			}
			results = stepwise_search_IM(data_b, formula, interactive_mode=FALSE, latent_var_original=latent_var_original_b, latent_var_extra=latent_var_extra_b, search_type=search_type, search=search, search_criterion=search_criterion, forward_exclude_p_bigger = forward_exclude_p_bigger, backward_exclude_p_smaller = backward_exclude_p_smaller, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var, eps=eps, maxiter=maxiter, family=family, seed=NULL, print=FALSE, remove_miss=TRUE)
			var_select_current = var_select
			for (i in 1:k){
				kept = names(stats::coef(results$fit_latent_var[[i]]))
				var_select_current[[i]] = names(var_select_current[[i]]) %in%  kept
			}
			return(var_select_current)
		}
		close(pb)
		snow::stopCluster(cl)
	})
	# Counting number of times variables appear
	var_select = results[[1]]
	# Keeping list of unique subset of variables and their count
	unique_subset_string = vector("character",boot_iter)
	unique_subset = vector("list", boot_iter)
	unique_subset_count = vector("numeric", boot_iter)

	unique_subset_string[1] = toString(results[[1]])
	unique_subset[[1]] = results[[1]]
	names(unique_subset[[1]]) = names(latent_var_original)
	unique_subset_count[1] = 1
	max_sub = 1

	for (i in 2:boot_iter){
		var_select = Map("+", var_select, results[[i]])
		index_subset = unique_subset_string %in% toString(results[[i]])
		if (sum(index_subset)==1) unique_subset_count[index_subset] = unique_subset_count[index_subset] + 1
		else{
			max_sub = max_sub + 1
			unique_subset_string[max_sub] = toString(results[[i]])
			unique_subset[[max_sub]] = results[[i]]
			unique_subset_count[max_sub] = 1
			names(unique_subset[[max_sub]]) = names(latent_var_original)
		}
	}
	new_index = order(unique_subset_count,decreasing=TRUE)
	unique_subset_count = unique_subset_count[new_index]
	unique_subset = unique_subset[new_index]
	names(var_select) = names(latent_var_original)
	for (i in 1:k){
		if (is.null(latent_var_extra[[i]])) latent_var_all = latent_var_original[[i]]
		else if (is.null(latent_var_original[[i]])) latent_var_all = latent_var_extra[[i]]
		else latent_var_all = cbind(latent_var_original[[i]],latent_var_extra[[i]])
		names(var_select[[i]]) = colnames(latent_var_all)
		new_index_current = order(var_select[[i]], decreasing = TRUE)
		var_select[[i]] = var_select[[i]][new_index_current]
		for (j in 1:max_sub){
			names(unique_subset[[j]][[i]]) = colnames(latent_var_all)
			unique_subset[[j]][[i]] = unique_subset[[j]][[i]][new_index_current]
		}
		var_select[[i]] = var_select[[i]]/boot_iter
	}
	unique_subset = unique_subset[1:min(max_sub,best_subsets)]
	unique_subset_count=unique_subset_count[1:min(max_sub,best_subsets)]
	unique_subset_count = unique_subset_count/boot_iter
	names(unique_subset) = paste0("Subset", seq(1,length(unique_subset)))
	cat("\nDone\n")
	return(list(var_select=var_select, subset_select=unique_subset_count, best_subset=unique_subset))
}


genetic_var_select = function(data, formula, parallel_iter=10, entropy_threshold=.10, popsize=25, mutation_prob=.50, first_pop=NULL, latent_var=NULL, search_criterion="AIC", maxgen = 100, eps=.01, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, progress=TRUE, n_cluster = 1, best_subsets=5, cv_iter=5, cv_folds=5, folds=NULL, Huber_p=1.345, classification=FALSE){
	k = length(latent_var)
	## Removing missing data and checks
	# Retaining only the needed variables from the dataset (need to set G and E variables for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	for (i in 1:k) data[,names(latent_var)[i]] = 0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	#  Check for empty latent variable (Note: Probably shouldn't be named empty_start_dataset but empty_start_latent_var althought that would be even more confusing considering start_latent_var)
	if (is.null(latent_var)) stop("latent_var cannot be null. However, if search=i, then you could set latent_var[[i]]=NULL.")
	# Make it to have NULL elements but not be NULL
	comp = rep(TRUE, NROW(data))
	for (i in 1:k){
		if (!is.null(latent_var[[i]])) comp = comp & stats::complete.cases(data, latent_var[[i]])
	}
	data = data[comp,, drop=FALSE]
	N_var = 0
	var = c()
	var_k = c()
	for (i in 1:k){
		if (!is.null(latent_var[[i]])) latent_var[[i]] = latent_var[[i]][comp,, drop=FALSE]
		N_var = NCOL(latent_var[[i]]) + N_var
		var = c(var, names(latent_var[[i]]))
		var_k = c(var_k, c(rep(i,NCOL(latent_var[[i]]))))
	}
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	# Setting up parallel
	cl <- snow::makeCluster(n_cluster)
	doSNOW::registerDoSNOW(cl)
	if (progress){
		if (runif(1) < .50) char="(^._.^) "
		else char="(=^o^=) "
		pb = utils::txtProgressBar(max = parallel_iter, style = 3, char=char)
		progress <- function(n) utils::setTxtProgressBar(pb, n)
		opts <- list(progress = progress)
	}
	else opts <- list()
	# Need to use this "with(c(),CODEHERE)" to prevent R check from returning a "no visible binding for global variable"
	with(c(),{
		results <- foreach::foreach(b = 1:parallel_iter, .options.snow = opts) %dopar% {
			if (!is.null(seed)) set.seed(seed+b)
			# The names of the variables selected is going to be the names of start_latent_var_pop elements
			start_latent_var_pop  = vector("list", popsize)
			# r if the frequency of appearance
			r = rep(0, N_var)
			indexes_pop = vector("list", popsize)
			# cirterion of the individual
			crit = rep(0, popsize)
			# Generate population
			for (p in 1:popsize){
				start_latent_var_pop[[p]] = vector("list", k)
				names(start_latent_var_pop[[p]]) = names(latent_var)
				# current latent_var (not to be retained, temporary)
				latent_var_pop = vector("list", k)
				names(latent_var_pop) = names(latent_var)
				for (i in 1:k){
					if (is.null(first_pop)){
						indexes = sample(1:NCOL(latent_var[[i]]),size=sample(1:NCOL(latent_var[[i]]),1))
						latent_var_pop[[i]] = latent_var[[i]][,indexes,drop=FALSE]
					}
					else{
						# As 0 0 0 1 0 0 1 coding
						indexes_ = colnames(latent_var[[i]]) %in% colnames(first_pop[[i]])
						# Mutate a lot
						mutated = FALSE
						while (!mutated){
							for (j in 1:length(indexes_)){
								if (runif(1) < .33) indexes_[j] = !indexes_[j]
							}
							# Revert back and restart if this remove all variables
							if (sum(indexes_)==0) indexes_ = colnames(latent_var[[i]]) %in% colnames(first_pop[[i]])
							else mutated = TRUE
						}
						# As 0 0 0 1 0 0 1, but we want c(4,7) instead
						indexes = indexes_*1:length(indexes_)
						indexes = indexes[indexes!=0]
						latent_var_pop[[i]] = latent_var[[i]][,indexes,drop=FALSE]
					}
					r = r + var %in% colnames(latent_var_pop[[i]])
					# Keeping indexes 
					indexes_pop[[p]] = c(indexes_pop[[p]], 1:NCOL(latent_var[[i]]) %in% indexes)
					# Starting point generation  (Dirichlet distribution)
					start_latent_var_pop[[p]][[i]] = rgamma(length(indexes),1)*(1-2*rbinom(length(indexes),1,.5))
					start_latent_var_pop[[p]][[i]] = start_latent_var_pop[[p]][[i]]/sum(abs(start_latent_var_pop[[p]][[i]]))
					names(start_latent_var_pop[[p]][[i]]) = colnames(latent_var_pop[[i]])
				}
				# Fit model
				fit = IMLEGIT(data=data, formula=formula, latent_var = latent_var_pop, start_latent_var=start_latent_var_pop[[p]], eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
				start_latent_var_pop[[p]] = lapply(fit$fit_latent_var,coef)
				if (search_criterion=="AIC") crit[p] = fit$true_model_parameters$AIC
				else if (search_criterion=="AICc") crit[p] = fit$true_model_parameters$AICc
				else if (search_criterion=="BIC") crit[p] = fit$true_model_parameters$BIC
				else{
					fit_cv = IMLEGIT_cv(data=data, formula=formula, latent_var = latent_var_pop, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var_pop[[p]], eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed)
					if (search_criterion=="cv") crit[p] = mean(fit_cv$R2_cv)	
					else if (search_criterion=="cv_Huber") crit[p] = mean(fit_cv$Huber_cv)
					else if (search_criterion=="cv_L1") crit[p] = mean(fit_cv$L1_cv)
					else if (search_criterion=="cv_AUC") crit[p] = mean(fit_cv$AUC)				
				}
			}
			# print(lapply(latent_var_pop,function(x) lapply(x,dim)))
			r_ = r/popsize
			entropy = - (1/N_var)*sum(r_*log2(r_) + (1-r_)*log2(1-r_), na.rm=TRUE)
			step = 0
			conv = TRUE
			while(entropy > entropy_threshold){
				step =  step + 1
				# Selection of the survivors
				if(search_criterion=="cv" || search_criterion=="cv_AUC") ordered_crit = order(crit, decreasing=TRUE)
				else ordered_crit = order(crit, decreasing=FALSE)
				survivors = ordered_crit[1:(popsize/2)]
				deaths = ordered_crit[((popsize/2)+1):popsize]
				# Reproduction, make babies
				indexes_pop_old = indexes_pop
				for (p in deaths){
					parents = sample(survivors,2)
					# Mutation index, do not mutate yet until children is born
					Mutation_index = sample(1:N_var,1)
					Mutation_yesorno = runif(1) < (mutation_prob)
					latent_var_pop = vector("list", k)
					names(latent_var_pop) = names(latent_var)
					for (i in 1:k){
						# crossover for latent variable i
						crossover = sample(which(var_k == i),length(which(var_k==i))/2)
						indexes_pop[[p]][var_k == i] = indexes_pop_old[[parents[1]]][var_k == i]
						indexes_pop[[p]][crossover] = indexes_pop_old[[parents[2]]][crossover]
						if (sum(indexes_pop[[p]][var_k == i])==0){
							# We have no variable in child, we need to reverse order of parenting to fix this
								indexes_pop[[p]][var_k == i] = indexes_pop_old[[parents[2]]][var_k == i]
								indexes_pop[[p]][crossover] = indexes_pop_old[[parents[1]]][crossover]
						}
						# Mutation
						if (var_k[Mutation_index]==i && Mutation_yesorno){
							indexes_pop[[p]][Mutation_index] = !indexes_pop[[p]][Mutation_index]
							# Revert back if this remove all variables
							if (sum(indexes_pop[[p]][var_k == i])==0) indexes_pop[[p]][Mutation_index] = !indexes_pop[[p]][Mutation_index]
						}
						# Removing current one from r
						r = r - (var %in% names(start_latent_var_pop[[p]][[i]]))
						# New latent var
						latent_var_pop[[i]] = latent_var[[i]][,indexes_pop[[p]][var_k == i],drop=FALSE]
						# Added new one to r
						r = r + var %in% colnames(latent_var_pop[[i]])
						# Index of the variables the child has
						b_index = which(indexes_pop[[p]][var_k==i])
						# Index of the variables the parents have
						p1_index = which(indexes_pop[[parents[1]]][var_k==i])
						p2_index = which(indexes_pop[[parents[2]]][var_k==i])
						start_latent_var_pop[[p]][[i]] = rep(NA, length(b_index))
						# Average starting point of parents (Not sure if this ideal? Could also use the splitting thing from birthing)
						for (j in b_index){
							if ((j %in% p1_index) && (j %in% p2_index)) start_latent_var_pop[[p]][[i]][b_index == j] = (start_latent_var_pop[[parents[1]]][[i]][p1_index == j] + start_latent_var_pop[[parents[2]]][[i]][p2_index == j])/2
							else if ((j %in% p1_index) && !(j %in% p2_index)) start_latent_var_pop[[p]][[i]][b_index == j] = start_latent_var_pop[[parents[1]]][[i]][p1_index == j]
							else if (!(j %in% p1_index) && (j %in% p2_index)) start_latent_var_pop[[p]][[i]][b_index == j] = start_latent_var_pop[[parents[2]]][[i]][p2_index == j]
							else start_latent_var_pop[[p]][[i]][b_index == j] = 0
						}
						start_latent_var_pop[[p]][[i]] = start_latent_var_pop[[p]][[i]]/sum(start_latent_var_pop[[p]][[i]])
					}
					# Fit model
					fit = IMLEGIT(data=data, formula=formula, latent_var = latent_var_pop, start_latent_var=start_latent_var_pop[[p]], eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
					start_latent_var_pop[[p]] = lapply(fit$fit_latent_var,coef)
					if (search_criterion=="AIC") crit[p] = fit$true_model_parameters$AIC
					else if (search_criterion=="AICc") crit[p] = fit$true_model_parameters$AICc
					else if (search_criterion=="BIC") crit[p] = fit$true_model_parameters$BIC
					else{
						fit_cv = IMLEGIT_cv(data=data, formula=formula, latent_var = latent_var_pop, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var_pop[[p]], eps=eps, maxiter=maxiter, family=family, ylim=ylim, seed=seed)
						if (search_criterion=="cv") crit[p] = mean(fit_cv$R2_cv)
						else if (search_criterion=="cv_Huber") crit[p] = mean(fit_cv$Huber_cv)
						else if (search_criterion=="cv_L1") crit[p] = mean(fit_cv$L1_cv)
						else if (search_criterion=="cv_AUC") crit[p] = mean(fit_cv$AUC)				
					}
				}
				r_ = r/popsize
				entropy = - (1/N_var)*sum(r_*log2(r_) + (1-r_)*log2(1-r_), na.rm=TRUE)
				if (step >= maxgen){
					conv = FALSE
					break
				}
			}
			# We won't return the whole population but just the % of times variables are selected(r) and top K models
			if (search_criterion=="cv" || search_criterion=="cv_AUC") ordered_crit = order(crit, decreasing=TRUE)
			else ordered_crit = order(crit, decreasing=FALSE)
			list_best_latent_var_crit = crit[ordered_crit[1:best_subsets]]
			list_best_latent_var_start = start_latent_var_pop[ordered_crit[1:best_subsets]]
			return(list(list_best_latent_var_crit=list_best_latent_var_crit, list_best_latent_var_start=list_best_latent_var_start,r=r, conv=conv, entropy=entropy, n_step=step))
		}
		close(pb)
		snow::stopCluster(cl)
	})
	list_best_latent_var_crit = results[[1]]$list_best_latent_var_crit
	list_best_latent_var_start = results[[1]]$list_best_latent_var_start
	list_entropy = results[[1]]$entropy
	list_n_step = results[[1]]$n_step
	# Keeping the best ones
	if (parallel_iter > 1){
		for (i in 2:parallel_iter){
			best_latent_crit = c(list_best_latent_var_crit, results[[i]]$list_best_latent_var_crit)
			if (search_criterion=="cv" || search_criterion=="cv_AUC") ordered_crit = order(best_latent_crit, decreasing=TRUE)
			else ordered_crit = order(best_latent_crit, decreasing=FALSE)
			list_best_latent_var_crit = c(list_best_latent_var_crit, results[[i]]$list_best_latent_var_crit)[ordered_crit[1:best_subsets]]
			list_best_latent_var_start = c(list_best_latent_var_start, results[[i]]$list_best_latent_var_start)[ordered_crit[1:best_subsets]]
			list_entropy = c(list_entropy, results[[i]]$entropy)
			list_n_step = c(list_n_step, results[[i]]$n_step)
		} 
	}
	# Getting percentages of variables used
	var_select = Reduce("+",lapply(results,function(x) x$r/popsize))/parallel_iter
	names(var_select) = var
	return(list(var_select=var_select, best_subsets_crit=list_best_latent_var_crit, best_subsets_start=list_best_latent_var_start, entropy=list_entropy, n_step=list_n_step))
}

nes_var_select = function(data, formula, parallel_iter=3, alpha=c(1,5,10), entropy_threshold=.05, popsize=25, lr = .2, prop_ignored=.50, latent_var=NULL, search_criterion="AICc", n_cluster=3, eps=.01, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, progress=TRUE, cv_iter=5, cv_folds=5, folds=NULL, Huber_p=1.345, classification=FALSE, print=FALSE){
	k = length(latent_var)
	## Removing missing data and checks
	# Retaining only the needed variables from the dataset (need to set G and E variables for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	for (i in 1:k) data[,names(latent_var)[i]] = 0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	#  Check for empty latent variable (Note: Probably shouldn't be named empty_start_dataset but empty_start_latent_var althought that would be even more confusing considering start_latent_var)
	if (is.null(latent_var)) stop("latent_var cannot be null. However, if search=i, then you could set latent_var[[i]]=NULL.")
	# Make it to have NULL elements but not be NULL
	comp = rep(TRUE, NROW(data))
	for (i in 1:k){
		if (!is.null(latent_var[[i]])) comp = comp & stats::complete.cases(data, latent_var[[i]])
	}
	data = data[comp,, drop=FALSE]
	N_var = 0
	var = c()
	var_k = c()
	for (i in 1:k){
		if (!is.null(latent_var[[i]])) latent_var[[i]] = latent_var[[i]][comp,, drop=FALSE]
		N_var = NCOL(latent_var[[i]]) + N_var
		var = c(var, names(latent_var[[i]]))
		var_k = c(var_k, c(rep(i,NCOL(latent_var[[i]]))))
	}
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	# The names of the variables selected is going to be the names of start_latent_var_pop elements
	start_latent_var_pop  = vector("list", popsize)
	# theta is the vector of parameters of the bernouilli distributions
	theta = vector("list", k)
	names(theta) = names(latent_var)
	for (i in 1:k){
		theta[[i]] = rep(.50, NCOL(latent_var[[i]]))
		names(theta[[i]]) = colnames(latent_var[[i]])
	}
	# Can't keep all X matrix because so big, so keeping the binary yes/no representing if a variable is kept in each subset of variable
	indexes_pop = vector("list", popsize)
	for (i in 1:popsize){
		indexes_pop[[i]] = vector("list", k)
		for (j in 1:k){
			indexes_pop[[i]][[j]] = rep(FALSE, NCOL(latent_var[[j]]))
		}
	}
	# cirterion of the individual
	crit = rep(0, popsize)
	crit_ranked = rep(0, popsize)
	theta_diff = Inf
	# Generate population
	conv = FALSE
	iter = 0
	# Setting up parallel
	cl <- snow::makeCluster(n_cluster, outfile="")
	doSNOW::registerDoSNOW(cl)
	if (progress){
		if (runif(1) < .50) char="(^._.^) "
		else char="(=^o^=) "
		pb = utils::txtProgressBar(max = parallel_iter, style = 3, char=char)
		progress <- function(n) utils::setTxtProgressBar(pb, n)
		opts <- list(progress = progress)
	}
	else opts <- list()
	# Need to use this "with(c(),CODEHERE)" to prevent R check from returning a "no visible binding for global variable"
	with(c(),{
		results <- foreach::foreach(b = 1:parallel_iter, .options.snow = opts) %dopar% {
			if (!is.null(seed)) set.seed(seed+b)
			entropy = Inf
			while(entropy > entropy_threshold){
				# r if the frequency of appearance
				r_entropy = rep(0, N_var)
				iter = iter + 1
				if (print){
					toprint = paste0("i=", iter)
					for (i in 1:k) toprint = paste0(toprint, " theta[",i,"]: ", paste0(round(theta[[i]],2), collapse=" "))
					print(toprint)
					flush.console()
				} 
				# grad is the gradient of the log-likelihood of the bernouilli distributions (has to be reinitialized every iteration)
				grad = vector("list", k)
				for (curr in 1:N_var){
					grad[[curr]] = rep(0, popsize)
				}
				for (p in 1:popsize){
					start_latent_var_pop[[p]] = vector("list", k)
					names(start_latent_var_pop[[p]]) = names(latent_var)
					# current latent_var (not to be retained, temporary)
					latent_var_pop = vector("list", k)
					names(latent_var_pop) = names(latent_var)
					curr = 0
					for (i in 1:k){
						# Can't have no variables so looping until we get a correct subset
						curr_backup = curr
						indexes = rep(FALSE, NCOL(latent_var[[i]]))
						while (sum(indexes)==0){
							indexes = rep(FALSE, NCOL(latent_var[[i]]))
							curr = curr_backup
							for (j in 1:NCOL(latent_var[[i]])){
								curr = curr + 1
								indexes[j] = rbinom(1, 1, prob = theta[[i]][j])==1
								if (indexes[j]) grad[[curr]][p] = 1/theta[[i]][j]
								else grad[[curr]][p] = -1/(1 - theta[[i]][j])
							}
						}
						latent_var_pop[[i]] = latent_var[[i]][,indexes,drop=FALSE]

						# Keeping indexes 
						indexes_pop[[p]][[i]] = indexes
						r_entropy = r_entropy + var %in% colnames(latent_var_pop[[i]])
						# Starting point (Dirichlet distribution)
						alpha_index = b %% length(b)
						if (alpha_index == 0) alpha_index = length(b)
						start_latent_var_pop[[p]][[i]] = rgamma(NCOL(latent_var_pop[[i]]),alpha[alpha_index])*(1-2*rbinom(NCOL(latent_var_pop[[i]]),1,.5))
						start_latent_var_pop[[p]][[i]] = start_latent_var_pop[[p]][[i]]/sum(abs(start_latent_var_pop[[p]][[i]]))
						names(start_latent_var_pop[[p]][[i]]) = colnames(latent_var_pop[[i]])
					}
					# Fit model
					fit = IMLEGIT(data=data, formula=formula, latent_var = latent_var_pop, start_latent_var=start_latent_var_pop[[p]], eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
					#fit = IMLEGIT(data=data, formula=formula, latent_var = latent_var_pop, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
					start_latent_var_pop[[p]] = lapply(fit$fit_latent_var,coef)
					if (search_criterion=="AIC") crit[p] = -fit$true_model_parameters$AIC
					else if (search_criterion=="AICc") crit[p] = -fit$true_model_parameters$AICc
					else if (search_criterion=="BIC") crit[p] = -fit$true_model_parameters$BIC
					else{
						fit_cv = IMLEGIT_cv(data=data, formula=formula, latent_var = latent_var_pop, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var_pop[[p]], eps=eps, maxiter=maxiter, family=family, ylim=ylim)
						if (search_criterion=="cv") crit[p] = mean(fit_cv$R2_cv)
						else if (search_criterion=="cv_Huber") crit[p] = -mean(fit_cv$Huber_cv)
						else if (search_criterion=="cv_L1") crit[p] = -mean(fit_cv$L1_cv)
						else if (search_criterion=="cv_AUC") crit[p] = mean(fit_cv$AUC)
					}
				}
				# Rank the criterion and apply fitness shaping function 2i-1 where i is rank in [0,1]
				# Make all low ranks to have the same value so that bad examples are not too penalized but good examples are more prioritized
				crit = (rank(crit)-1)/(length(crit)-1)*2-1
				crit[crit < prop_ignored*2-1] = prop_ignored*2-1
				# Gradient descent
				curr = 0
				for (i in 1:k){
					for (j in 1:NCOL(latent_var[[i]])){
						curr = curr + 1
						grad_theta = mean(crit*grad[[curr]])
						theta[[i]][j] = max(min(theta[[i]][j] + lr*grad_theta,1),0)

					}
				}
				r_ = r_entropy/popsize
				entropy = - (1/N_var)*sum(r_*log2(r_) + (1-r_)*log2(1-r_), na.rm=TRUE)
				if (print) print(entropy)
			}
			p = order(crit, decreasing=TRUE)[1]
			latent_var_best = vector("list", k)
			names(latent_var_best) = names(latent_var)
			for (i in 1:k) latent_var_best[[i]] = latent_var[[i]][,indexes_pop[[p]][[i]],drop=FALSE]
			start_latent_var_best = start_latent_var_pop[[p]]
			return(list(latent_var_best=latent_var_best,start_latent_var_best=start_latent_var_best, crit_best=max(crit)))
		}
		close(pb)
		snow::stopCluster(cl)
	})
	latent_var_best = results[[1]]$latent_var_best
	start_latent_var_best = results[[1]]$start_latent_var_best
	crit = results[[1]]$crit_best
	# Keeping the best ones
	if (parallel_iter > 1){
		for (i in 2:parallel_iter){
			if (results[[i]]$crit_best > crit){
				latent_var_best = results[[i]]$latent_var_best
				start_latent_var_best = results[[i]]$start_latent_var_best
				crit = results[[i]]$crit_best
			}
		}
	}
	best_fit = IMLEGIT(data=data, formula=formula, latent_var = latent_var_best, start_latent_var=start_latent_var_best, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
	best_fit_cv = IMLEGIT_cv(data=data, formula=formula, latent_var = latent_var_best, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var_best, eps=eps, maxiter=maxiter, family=family, ylim=ylim)
	return(list(fit=best_fit, fit_cv=best_fit_cv, data=data, latent_var=latent_var_best, start_latent_var=start_latent_var_best))
}

r1nes_var_select = function(data, formula, parallel_iter=3, alpha=c(1,5,10), entropy_threshold=.05, popsize=25, lr = .2, prop_ignored=.50, latent_var=NULL, search_criterion="AICc", n_cluster=3, eps=.01, maxiter=100, family=gaussian, ylim=NULL, seed=NULL, progress=TRUE, cv_iter=5, cv_folds=5, folds=NULL, Huber_p=1.345, classification=FALSE, print=FALSE){
	k = length(latent_var)
	## Removing missing data and checks
	# Retaining only the needed variables from the dataset (need to set G and E variables for this to work, they will be replaced with their proper values later)
	data=data.frame(data)
	for (i in 1:k) data[,names(latent_var)[i]] = 0
	data = stats::model.frame(formula, data=data, na.action=na.pass)
	#  Check for empty latent variable (Note: Probably shouldn't be named empty_start_dataset but empty_start_latent_var althought that would be even more confusing considering start_latent_var)
	if (is.null(latent_var)) stop("latent_var cannot be null. However, if search=i, then you could set latent_var[[i]]=NULL.")
	# Make it to have NULL elements but not be NULL
	comp = rep(TRUE, NROW(data))
	for (i in 1:k){
		if (!is.null(latent_var[[i]])) comp = comp & stats::complete.cases(data, latent_var[[i]])
	}
	data = data[comp,, drop=FALSE]
	N_var = 0
	index_var = c() # Number determine in which latent variable the observed variable is
	var = c()
	for (i in 1:k){
		if (!is.null(latent_var[[i]])) latent_var[[i]] = latent_var[[i]][comp,, drop=FALSE]
		N_var = NCOL(latent_var[[i]]) + N_var
		var = c(var, names(latent_var[[i]]))
		index_var = c(index_var, rep(i, NCOL(latent_var[[i]])))
	}
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	# The names of the variables selected is going to be the names of start_latent_var_pop elements
	start_latent_var_pop  = vector("list", popsize)
	## Parameters of the Normal distribution with mean mu and covariance (sigma^2(I + vv^T))
	# mu (Mean)
	# s (e^2s = sigma^2)
	# c (e^c = length of u)
	# z (||z||=1, z is direction of u)
	# v (z*e^c, principal direction used to form covariance matrix)
	mu = rep(0, N_var)
	s = 0 # Lead to sigma^2 = 1
	c = 0 # Lead to length of 1
	z = rep(1/sqrt(N_var), N_var) # sqrt((1/sqrt(n))^2 + ... + (1/sqrt(n))^2) = sqrt(n * (1/n)) = sqrt(1) = 1
	v = z
	# Can't keep all X matrix because so big, so keeping the binary yes/no representing if a variable is kept in each subset of variable
	indexes_pop = vector("list", popsize)
	for (i in 1:popsize){
		indexes_pop[[i]] = vector("list", k)
		for (j in 1:k){
			indexes_pop[[i]][[j]] = rep(FALSE, NCOL(latent_var[[j]]))
		}
	}
	# cirterion of the individual
	crit = rep(0, popsize)
	crit_prev = rep(0, popsize)
	crit_ranked = rep(0, popsize)
	mu_diff = Inf
	# Generate population
	conv = FALSE
	iter = 0
	# Setting up parallel
	cl <- snow::makeCluster(n_cluster, outfile="")
	doSNOW::registerDoSNOW(cl)
	if (progress){
		if (runif(1) < .50) char="(^._.^) "
		else char="(=^o^=) "
		pb = utils::txtProgressBar(max = parallel_iter, style = 3, char=char)
		progress <- function(n) utils::setTxtProgressBar(pb, n)
		opts <- list(progress = progress)
	}
	else opts <- list()
	# Need to use this "with(c(),CODEHERE)" to prevent R check from returning a "no visible binding for global variable"
	with(c(),{
		results <- foreach::foreach(b = 1:parallel_iter, .options.snow = opts) %dopar% {
			if (!is.null(seed)) set.seed(seed+b)
			entropy = Inf
			while(entropy > entropy_threshold){
				# r if the frequency of appearance
				r_entropy = rep(0, N_var)
				iter = iter + 1
				if (print){
					toprint = paste0("i=", iter)
					for (i in 1:N_var) toprint = paste0(toprint, " mu[",i,"]: ", paste0(round(mu[i],2), collapse=" "))
					for (i in 1:N_var) toprint = paste0(toprint, " v[",i,"]: ", paste0(round(v[i],2), collapse=" "))
					toprint = paste0(toprint, " s: ", paste0(round(exp(s),2), collapse=" "))
					print(toprint)
					flush.console()
				}
				# List with the population generated samples (Need to keep them for gradient descent)
				X_pop = vector("list", popsize)
				W_pop = vector("list", popsize)
				for (p in 1:popsize){
					start_latent_var_pop[[p]] = vector("list", k)
					names(start_latent_var_pop[[p]]) = names(latent_var)
					# current latent_var (not to be retained, temporary)
					latent_var_pop = vector("list", k)
					names(latent_var_pop) = names(latent_var)
					curr = 0
					# Need to generate example and make sure its proper (i.e. that there is no empty latent variable)
					proper_example = FALSE
					while(!proper_example){
						proper_example = TRUE
						y = rnorm(N_var)
						a = rnorm(1)
						w = y+a*v
						W_pop[[p]] = w
						X_pop[[p]] = exp(s)*w + mu
						# Cutoff here to make into binary (variable selected, yes/no)
						indexes = X_pop[[p]] > 0
						for (i in 1:k) if(sum(indexes[index_var==i]) == 0) proper_example = FALSE
					}
					# Now that we have sample, we translate it into proper format and we generate a starting point
					for (i in 1:k){
						latent_var_pop[[i]] = latent_var[[i]][,indexes[index_var==i],drop=FALSE]
						# Keeping indexes 
						indexes_pop[[p]][[i]] = indexes[index_var==i]
						r_entropy = r_entropy + var %in% colnames(latent_var_pop[[i]])
						# Starting point (Dirichlet distribution)
						alpha_index = b %% length(b)
						if (alpha_index == 0) alpha_index = length(b)
						start_latent_var_pop[[p]][[i]] = rgamma(NCOL(latent_var_pop[[i]]),alpha[alpha_index])*(1-2*rbinom(NCOL(latent_var_pop[[i]]),1,.5))
						start_latent_var_pop[[p]][[i]] = start_latent_var_pop[[p]][[i]]/sum(abs(start_latent_var_pop[[p]][[i]]))
						names(start_latent_var_pop[[p]][[i]]) = colnames(latent_var_pop[[i]])
					}
					# Fit model
					fit = IMLEGIT(data=data, formula=formula, latent_var = latent_var_pop, start_latent_var=start_latent_var_pop[[p]], eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
					#fit = IMLEGIT(data=data, formula=formula, latent_var = latent_var_pop, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
					start_latent_var_pop[[p]] = lapply(fit$fit_latent_var,coef)
					if (search_criterion=="AIC") crit[p] = -fit$true_model_parameters$AIC
					else if (search_criterion=="AICc") crit[p] = -fit$true_model_parameters$AICc
					else if (search_criterion=="BIC") crit[p] = -fit$true_model_parameters$BIC
					else{
						fit_cv = IMLEGIT_cv(data=data, formula=formula, latent_var = latent_var_pop, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var_pop[[p]], eps=eps, maxiter=maxiter, family=family, ylim=ylim)
						if (search_criterion=="cv") crit[p] = mean(fit_cv$R2_cv)
						else if (search_criterion=="cv_Huber") crit[p] = -mean(fit_cv$Huber_cv)
						else if (search_criterion=="cv_L1") crit[p] = -mean(fit_cv$L1_cv)
						else if (search_criterion=="cv_AUC") crit[p] = mean(fit_cv$AUC)
					}
				}
				# Rank the criterion and apply fitness shaping function 2i-1 where i is rank in [0,1]
				# Make all low ranks to have the same value so that bad examples are not too penalized but good examples are more prioritized
				crit_new = crit
				crit = (rank(crit)-1)/(length(crit)-1)*2-1
				crit[crit < prop_ignored*2-1] = prop_ignored*2-1
				### Natural gradient descent
				# initialize gradients of J at 0
				mu_diff = rep(0, N_var)
				s_diff = 0
				c_diff = 0
				z_diff = rep(0, N_var)
				v_diff = rep(0, N_var)
				r_2 = sum(v^2)
				r = sqrt(r_2)
				for (p in 1:popsize){
					x = X_pop[[p]]
					w = W_pop[[p]]
					z = v/r
					w_w = c(w %*% w)
					w_z = c(w %*% z)
					# Gradients
					grad_mu = (x - mu)
					grad_s = (1/(2*(N_var-1)))*((w_w-N_var) - (((w_z)^2)-1))
					grad_v = ((1/(2*(N_var-1)*r))*((r_2 - N_var + 2)*((w_z)^2) - ((r_2+1)*w_w)))*z + ((w_z)/r)*w
					grad_c = (1/r)*c(grad_v %*% z)
					grad_z = (1/r)*(grad_v - c(grad_v %*% z)*z)
					# Gradient of J (diff)
					mu_diff = mu_diff + crit[p]*grad_mu
					s_diff = s_diff + crit[p]*grad_s
					v_diff = v_diff + crit[p]*grad_v
					c_diff = c_diff + crit[p]*grad_c
					z_diff = z_diff + crit[p]*grad_z
				}
				mu_diff = lr*mu_diff/popsize
				mu = mu + mu_diff
				s_diff = lr*s_diff/popsize
				s = s + s_diff
				c_diff = lr*c_diff/popsize
				v_diff = lr*v_diff/popsize
				z_diff = lr*z_diff/popsize
				if (c_diff < 0){
					c = c + c_diff
					z_new = z + z_diff
					z = z_new/sqrt(sum(z_new^2))
					v_old = v
					v = exp(c)*z
					v_diff = v - v_old
				}
				else{
					v = v + v_diff
					v_norm = sqrt(sum(v^2))
					c = log(v_norm)
					z = v/v_norm
				}
				r_ = r_entropy/popsize
				entropy = - (1/N_var)*sum(r_*log2(r_) + (1-r_)*log2(1-r_), na.rm=TRUE)
				if (print) print(entropy)
			}
			p = order(crit, decreasing=TRUE)[1]
			latent_var_best = vector("list", k)
			names(latent_var_best) = names(latent_var)
			for (i in 1:k) latent_var_best[[i]] = latent_var[[i]][,indexes_pop[[p]][[i]],drop=FALSE]
			start_latent_var_best = start_latent_var_pop[[p]]
			return(list(latent_var_best=latent_var_best,start_latent_var_best=start_latent_var_best, crit_best=max(crit)))
		}
		close(pb)
		snow::stopCluster(cl)
	})
	latent_var_best = results[[1]]$latent_var_best
	start_latent_var_best = results[[1]]$start_latent_var_best
	crit = results[[1]]$crit_best
	# Keeping the best ones
	if (parallel_iter > 1){
		for (i in 2:parallel_iter){
			if (results[[i]]$crit_best > crit){
				latent_var_best = results[[i]]$latent_var_best
				start_latent_var_best = results[[i]]$start_latent_var_best
				crit = results[[i]]$crit_best
			}
		}
	}
	best_fit = IMLEGIT(data=data, formula=formula, latent_var = latent_var_best, start_latent_var=start_latent_var_best, eps=eps, maxiter=maxiter, family=family, ylim=ylim, print=FALSE)
	best_fit_cv = IMLEGIT_cv(data=data, formula=formula, latent_var = latent_var_best, cv_iter=cv_iter, cv_folds=cv_folds, folds=folds, Huber_p=Huber_p, classification=classification, start_latent_var=start_latent_var_best, eps=eps, maxiter=maxiter, family=family, ylim=ylim)
	return(list(fit=best_fit, fit_cv=best_fit_cv, data=data, latent_var=latent_var_best, start_latent_var=start_latent_var_best))
}