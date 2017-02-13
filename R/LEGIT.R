#' @title Simulated example of a 2 way interaction GxE model.
#' @description Simulated example of a 2 way interaction GxE model. 
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
#' @return Returns a list containing, in the following order: data.frame with the observed outcome (with noise) and the true outcome (without noise), data.frame of the genetic variants, data.frame of the environments, vector of the true genetic coefficients, vector of the true environmental coefficients, vector of the true main model coefficients
#' @examples
#'	example_2way(5,1,logit=FALSE)
#'	example_2way(5,0,logit=TRUE)
#' @export
"example_2way"

#' @title Simulated example of a 3 way interaction GxE model
#' @description Simulated example of a 3 way interaction GxE model. 
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
#' @return Returns a list containing, in the following order: data.frame with the observed outcome (with noise), the true outcome (without noise) and \eqn{z}, data.frame of the genetic variants, data.frame of the environments, vector of the true genetic coefficients, vector of the true environmental coefficients, vector of the true main model coefficients
#' @examples
#'	example_3way(5,2.5,logit=FALSE)
#'	example_3way(5,0,logit=TRUE)
#' @export
"example_3way"

#' @title Latent Environmental & Genetic InTeraction (LEGIT) model
#' @description Constructs a generalized linear model (glm) with a weighted latent environmental score and weighted latent genetic score.
#' @param data data.frame of the dataset to be used. Do not include elements that are in the datasets \code{genes} and \code{env} and do not include manually coded interactions.
#' @param genes data.frame of the variables inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic).
#' @param env data.frame of the variables inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental).
#' @param formula Model formula. Use \emph{E} for the environmental score and \emph{G} for the genetic score. Do not manually code interactions, write them in the formula instead (ex: G*E*z or G:E:z).
#' @param start_genes Optional starting points for genetic score (must be same length as the number of columns of \code{genes}).
#' @param start_env Optional starting points for environmental score (must be same length as the number of columns of \code{env}).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param print If FALSE, nothing except warnings will be printed. (Default = TRUE).
#' @return Returns an object of the class "LEGIT" which is list containing, in the following order: a glm fit of the main model, a glm fit of the genetic score, a glm fit of the environmental score, a list of the true model parameters (AIC, BIC, rank, df.residual, null.deviance) for which the individual model parts (main, genetic, environmental) don't estimate properly.
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
#' @export
"LEGIT"

#' @title Predictions of LEGIT fits
#' @description Predictions of LEGIT fits.
#' @param object An object of class "LEGIT", usually, a result of a call to LEGIT.
#' @param data data.frame of the dataset to be used. Do not include elements that are in the datasets \code{genes} and \code{env} and do not include manually coded interactions.
#' @param genes data.frame of the variables inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic).
#' @param env data.frame of the variables inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental).
#' @param ... Further arguments passed to or from other methods.
#' @return Returns a vector with the predicted values.
#' @examples
#'	train = example_2way(250, 1, seed=777)
#'	test = example_2way(100, 1, seed=777)
#'	fit = LEGIT(train$data, train$G, train$E, y ~ G*E)
#'	ssres = sum((test$data$y - predict(fit, test$data, test$G, test$E))^2)
#'	sstotal = sum((test$data$y - mean(test$data$y))^2)
#'	R2 = 1 - ssres/sstotal
#' @export
"predict.LEGIT"

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

#' @title Cross-validation for the LEGIT model
#' @description Uses cross-validation on the LEGIT model. Note that this is not a very fast implementation since it was written in R.
#' @param data data.frame of the dataset to be used. Do not include elements that are in the datasets \code{genes} and \code{env} and do not include manually coded interactions.
#' @param genes data.frame of the variables inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic).
#' @param env data.frame of the variables inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental).
#' @param formula Model formula. Use \emph{E} for the environmental score and \emph{G} for the genetic score. Do not manually code interactions, write them in the formula instead (ex: G*E*z or G:E:z).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param classification Set to TRUE if you are doing classification (binary outcome).
#' @param start_genes Optional starting points for genetic score (must be same length as the number of columns of \code{genes}).
#' @param start_env Optional starting points for environmental score (must be same length as the number of columns of \code{env}).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param seed Seed for cross-validation folds.
#' @return If \code{classification} = FALSE, returns a list containing, in the following order: a vector of the cross-validated \eqn{R^2} at each iteration, a matrix of the possible outliers (standardized residuals > 2.5 or < -2.5) and their corresponding standardized residuals and standardized pearson residuals. If \code{classification} = TRUE, returns a list containing, in the following order: a vector of the cross-validated \eqn{R^2} at each iteration, a vector of the AUC at each iteration, a matrix of the best choice of threshold (based on Youden index) and the corresponding specificity and sensitivity at each iteration and a list of objects of class "roc" (to be able to make roc curve plots) at each iteration.
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
#' @export
"LEGIT_cv"

#' Internal function that does the forward step for the stepwise function.
#' @param empty_start_dataset If TRUE, the initial dataset is empty.
#' @param fit Current best fit.
#' @param ... Same parameters as in the stepwise function.
#' @return Returns fit, start_genes, start_env and genes_current, genes_toadd if search="genes" or env_current and env_toadd if search="env".
#' @keywords internal
"forward_step"

#' Internal function that does the backward step for the stepwise function.
#' @param empty_start_dataset If TRUE, the initial dataset is empty.
#' @param fit Current best fit.
#' @param ... Same parameters as in the stepwise function.
#' @return Returns fit, start_genes, start_env and genes_current, genes_toadd if search="genes" or env_current and env_toadd if search="env".
#' @keywords internal
"backward_step"

#' @title Stepwise search for the best subset of genetic variants or environments with the LEGIT model
#' @description Adds the best variable or drops the worst variable one at a time in the genetic (if \code{search="genes"}) or environmental score (if \code{search="env"}). For now, only \code{search_type="forward"} is implemented. You can select the desired search criterion (AIC, BIC, cross-validation error, cross-validation AUC) to determine which variable is the best/worst and should be added/dropped. If using cross-validation (\code{search_criterion="cv"} or \code{search_criterion="cv_AUC"}), to prevent cross-validating with each variable (extremely slow), we recommend setting a p-value threshold (\code{p_threshold}) and forcing the algorithm not to look at models with bigger AIC (\code{exclude_worse_AIC=TRUE}).
#' @param data data.frame of the dataset to be used. Do not include elements that are in the datasets \code{genes} and \code{env} and do not include manually coded interactions.
#' @param formula Model formula. Use \emph{E} for the environmental score and \emph{G} for the genetic score. Do not manually code interactions, write them in the formula instead (ex: G*E*z or G:E:z).
#' @param interactive_mode If TRUE, uses interactive mode. In interactive mode, at each iteration, the user is shown the AIC, BIC, p-value and also the cross-validation \eqn{R^2} if \code{search_criterion="cv"} and the cross-validation AUC if \code{search_criterion="cv_AUC"} for the best 5 variables. The user must then enter a number between 1 and 5 to select the variable to be added, entering anything else will stop the search.
#' @param genes_original data.frame of the variables inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic).
#' @param env_original data.frame of the variables inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental).
#' @param genes_extra data.frame of the additionnal variables to try including inside the genetic score \emph{G} (can be any sort of variable, doesn't even have to be genetic). If not NULL, \code{env_extra} should be NULL.
#' @param env_extra data.frame of the variables to try including inside the environmental score \emph{E} (can be any sort of variable, doesn't even have to be environmental). If not NULL, \code{genes_extra} should be NULL.
#' @param search_type If \code{search_type="forward"}, uses a forward search. If \code{search_type="bidirectional"}, uses bidirectional forward search. If \code{search_type="backward"}, uses backward search. For now only \code{search_type="forward"} is implemented (Default = "forward").
#' @param search If \code{search="genes"}, uses a stepwise search for the genetic score variables \code{genes_extra}, forcing \code{genes_original} to be included in the genetic score. If \code{search="env"}, uses a stepwise search for the environmental score variables \code{env_extra}, forcing \code{env_original} to be included in the genetic score (Default = "genes").
#' @param search_criterion Criterion used to determine which variable is the best to add or worst to drop. if \code{search_criterion="AIC"}, uses the AIC, if \code{search_criterion="BIC"}, uses the BIC, if \code{search_criterion="cv"}, uses the cross-validation error, if \code{search_criterion="cv_AUC"}, uses the cross-validated AUC (Default = "AIC").
#' @param forward_exclude_p_bigger If p-value > \code{forward_exclude_p_bigger}, we do not consider the variable for inclusion in the forward steps (Default = .20).
#' @param backward_exclude_p_smaller If p-value < \code{backward_exclude_p_smaller}, we do not consider the variable for removal in the backward steps (Default = .01).
#' @param exclude_worse_AIC If AIC with variable > AIC without variable, we ignore the variable (Default = TRUE).
#' @param max_steps Maximum number of steps taken (Default = 50).
#' @param cv_iter Number of cross-validation iterations (Default = 5).
#' @param cv_folds Number of cross-validation folds (Default = 10). Using \code{cv_folds=NROW(data)} will lead to leave-one-out cross-validation.
#' @param classification Set to TRUE if you are doing classification (binary outcome).
#' @param start_genes Optional starting points for genetic score (must be same length as the number of columns of \code{genes}).
#' @param start_env Optional starting points for environmental score (must be same length as the number of columns of \code{env}).
#' @param eps Threshold for convergence (.01 for quick batch simulations, .0001 for accurate results).
#' @param maxiter Maximum number of iterations.
#' @param family Outcome distribution and link function (Default = gaussian).
#' @param seed Seed for cross-validation folds.
#' @param print If TRUE, print all the steps and notes/warnings. Highly recommended unless you are batch running multiple stepwise searchs. (Default=TRUE).
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
	return(list(data=data.frame(y,y_true),G=data.frame(g1,g2,g3,g4,g1_g3,g2_g3),E=data.frame(e1,e2,e3),coef_G=c(.2,.15,-.3,.1,.05,.2),coef_E=c(-.45,.35,.2), coef_main=c(5,2,3,4)))
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
	return(list(data=data.frame(y,y_true,z),G=data.frame(g1,g2,g3,g4,g1_g3,g2_g3),E=data.frame(e1,e2,e3),coef_G=c(.2,.15,-.3,.1,.05,.2),coef_E=c(-.45,.35,.2), coef_main=c(5,2,3,1,5,1.5,2,2)))
}


LEGIT = function(data, genes, env, formula, start_genes=NULL, start_env=NULL, eps=.001, maxiter=50, family=gaussian, print=TRUE)
{
	if (maxiter <= 0) warning("maxiter must be > 0")
	if(!is.null(start_genes)){
		if (NCOL(genes)!=length(start_genes)) stop("start_genes must either be NULL or have the same length as the number of genes")
		}
	if(!is.null(start_env)){
		if (NCOL(env)!=length(start_env)) stop("start_env must either be NULL or have the same length as the number of environments")
	}
	if (class(data) != "data.frame" && class(data) != "matrix") stop("data must be a data.frame")

	# getting right formats
	data = data.frame(data)
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

	# remove missing data
	comp = stats::complete.cases(data,genes,env)
	data = data[comp,, drop=FALSE]
	genes = genes[comp,, drop=FALSE]
	env = env[comp,, drop=FALSE]
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	#Adding empty variables in main dataset for genes and env
	data[,colnames(genes)]=0
	data[,colnames(env)]=0
	data$R0_b=0
	data$R0_c=0

	# Setting up initial weighted scores
	if (is.null(start_genes)) weights_genes = rep(1/dim(genes)[2],dim(genes)[2])
	else if (sum(abs(start_genes))==0) weights_genes = rep(1/dim(genes)[2],dim(genes)[2])
	else weights_genes = start_genes/sum(abs(start_genes))

	if (is.null(start_env)) weights_env = rep(1/dim(env)[2],dim(env)[2])
	else if (sum(abs(start_env))==0) weights_env = rep(1/dim(env)[2],dim(env)[2])
	else weights_env = start_env/sum(abs(start_env))

	data$G = genes%*%weights_genes
	data$E = env%*%weights_env

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
	if (!(grepl("1",formula_elem_withG, fixed=TRUE) && TRUE)) formula_withG = paste0(formula_withG, " - 1")
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
	if (!(grepl("1",formula_elem_withE, fixed=TRUE) && TRUE)) formula_withE = paste0(formula_withE, " - 1")
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

			# Updating G estimates and checking convergence
			weights_genes_old = weights_genes
			weights_genes = weights_genes_/sum(abs(weights_genes_))
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

			## Step c : fit model for E
			fit_c = stats::glm(formula_c, data=data, family=family, y=FALSE, model=FALSE)
			weights_env_ = stats::coef(fit_c)

			# Updating E estimates and checking convergence
			weights_env_old = weights_env
			weights_env = weights_env_/sum(abs(weights_env_))
			data$E = env%*%weights_env
			if(sqrt(sum((weights_env_old-weights_env)^2)) < eps) conv_E = TRUE
			else conv_E = FALSE
		}
		else conv_E = TRUE

		if (conv_G & conv_E) break
	}

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

	fit_c = stats::glm(formula_c, data=data, family=family, y=FALSE, model=FALSE)
	data[,colnames(env)] = data[,colnames(env)]*sum(abs(stats::coef(fit_c)))
	fit_c = stats::glm(formula_c, data=data, family=family, y=FALSE, model=FALSE)

	if (!(((fit_a$deviance-fit_b$deviance)/fit_a$deviance)<.01 && ((fit_a$deviance-fit_b$deviance)/fit_a$deviance)<.01)) warning("Deviance differs by more than 1% between model parts. Make sure that everything was set up properly and try increasing the number of iterations (maxiter).")

	#Change some arguments so that we get the right AIC, BIC and dispersion for the model
	true_aic = fit_a$aic + 2*(fit_b$rank - 1) + 2*(fit_c$rank - 1)
	true_rank = fit_a$rank + (fit_b$rank - 1) + (fit_c$rank - 1)
	true_bic = true_aic - 2*true_rank + log(fit_a$df.null+1)*true_rank
	true_df.residual = (fit_a$df.null+1) - true_rank
	true_null.deviance = fit_a$null.deviance

	# print convergences stuff;
	if (conv_G & conv_E){
		if (print) cat(paste0("Converged in ",i, " iterations\n"))
	} 
	else{
		warning(paste0("Did not reach convergence in maxiter iterations. Try increasing maxiter or make eps smaller."))
	}

	result = list(fit_main = fit_a, fit_genes = fit_b, fit_env = fit_c, true_model_parameters=list(AIC = true_aic, BIC = true_bic, rank = true_rank, df.residual = true_df.residual, null.deviance=true_null.deviance))
	class(result) <- "LEGIT"
	return(result)
}

predict.LEGIT = function(object, data, genes, env, ...){
	data = data.frame(data)
	genes = as.matrix(genes, drop=FALSE)
	env = as.matrix(env, drop=FALSE)
	data$G = genes%*%stats::coef(object[[2]])
	data$E = env%*%stats::coef(object[[3]])
	return(stats::predict.glm(object[[1]], newdata=data, ...))
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

LEGIT_cv = function (data, genes, env, formula, cv_iter=5, cv_folds=10, classification=FALSE, start_genes=NULL, start_env=NULL, eps=.001, maxiter=50, family=gaussian, seed=NULL){

	# getting right formats
	data = data.frame(data)
	genes = as.matrix(genes, drop=FALSE)
	env = as.matrix(env, drop=FALSE)
	formula = stats::as.formula(formula)

	# remove missing data
	comp = stats::complete.cases(data,genes,env)
	data = data[comp,, drop=FALSE]
	genes = genes[comp,, drop=FALSE]
	env = env[comp,, drop=FALSE]
	if (dim(data)[1] <= 0) stop("no valid observation without missing values")

	formula_outcome = get.vars(formula)[1]
	R2_cv = c()
	AUC = c()
	best_threshold = c()
	roc_curve = list()
	residuals = rep(0,dim(data)[1])
	pearson_residuals = rep(0,dim(data)[1])

	for (j in 1:cv_iter){
		if (!is.null(seed)) set.seed(seed*j)
		# Folds
		s = sample(NROW(data))
		data_n = data[s,, drop=FALSE]
		genes_n = genes[s,, drop=FALSE]
		env_n = env[s,, drop=FALSE]
		id = cut(seq(1,NROW(data_n)),breaks=cv_folds,labels=FALSE)
		list = 1:cv_folds

		pred=c()
		y_test=c()

		for (i in 1:cv_folds){
			# Train and test datasets
			data_train = subset(data_n, id %in% list[-i], drop = FALSE)
			genes_train = subset(genes_n, id %in% list[-i], drop = FALSE)
	 		env_train = subset(env_n, id %in% list[-i], drop = FALSE)
	 		data_test = subset(data_n, id %in% c(i), drop = FALSE)
	 		genes_test = subset(genes_n, id %in% c(i), drop = FALSE)
	 		env_test = subset(env_n, id %in% c(i), drop = FALSE)
	 		y_test_new = data_test[,formula_outcome]

	 		# Fit model and add predictions
	 		fit_train = LEGIT(data=data_train, genes=genes_train, env=env_train, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, print=FALSE)
			pred_new = predict(fit_train, data=data_test,genes=genes_test,env=env_test,type="response")
			pred = c(pred,pred_new)
			y_test = c(y_test, y_test_new)
		}

		# Cross-validated R2
		ssres = sum((pred-y_test)^2)
		sstotal = sum((y_test-mean(y_test))^2)
		R2_cv = c(R2_cv, 1 - ssres/sstotal)

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
	possible_outliers_data = cbind(rownames(data)[possible_outliers],residuals[possible_outliers],pearson_residuals[possible_outliers])
	colnames(possible_outliers_data)=c("Observation","Standardized_residual","Standardized_pearson_residual")

	if (classification) return(list(R2_cv = R2_cv, AUC=AUC, best_threshold=best_threshold, roc_curve = roc_curve, possible_outliers = possible_outliers_data))
	return(list(R2_cv = R2_cv, possible_outliers = possible_outliers_data))
}


forward_step = function(empty_start_dataset, fit, data, formula, interactive_mode=FALSE, genes_current=NULL, env_current=NULL, genes_toadd=NULL, env_toadd=NULL, search="genes", search_criterion="AIC", p_threshold = .20, exclude_worse_AIC=TRUE, max_steps = 50, cv_iter=5, cv_folds=10, classification=FALSE, start_genes=NULL, start_env=NULL, eps=.01, maxiter=25, family=gaussian, seed=NULL, print=TRUE){
	# How much genes or env to add
	#print(head(data))
	#print(dim(data))
	#print(head(env_current))
	#print(dim(env_current))
	#print(head(env_toadd))
	#print(dim(env_toadd))
	if (search=="genes") elements_N = NCOL(genes_toadd)
	if (search=="env") elements_N = NCOL(env_toadd)
	if (elements_N == 0) return(NULL)
	# Vector which says which extra variables are "good" (worth exploring)
	good = rep(TRUE, elements_N)
	# Vector which says how much the criterion changed from including the variable
	criterion_before = rep(0, elements_N)
	criterion_after = rep(0, elements_N)
	criterion_diff = rep(0, elements_N)
	# In interactive model, we must keep track of every AIC, BIC, p-value, Cross-validated R2 and AUC to show the user at every iteration
	if (interactive_mode){
		if (search=="genes") interactive = data.frame(variable=colnames(genes_toadd), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N),cv_R2_old=rep(NA, elements_N),cv_R2_new=rep(NA, elements_N),cv_AUC_old=rep(NA, elements_N),cv_AUC_new=rep(NA, elements_N))
		if (search=="env") interactive = data.frame(variable=colnames(env_toadd), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N),cv_R2_old=rep(NA, elements_N),cv_R2_new=rep(NA, elements_N),cv_AUC_old=rep(NA, elements_N),cv_AUC_new=rep(NA, elements_N))
	}
	# Complete dataset
	if (NCOL(genes_current)==0) genes_current_nomiss = NULL
	else genes_current_nomiss = genes_current
	if (NCOL(env_current)==0) env_current_nomiss = NULL
	else env_current_nomiss = env_current
	comp_without = stats::complete.cases(data,genes_current_nomiss,env_current_nomiss)
	# Non-cross-validated models
	for (j in 1:elements_N){
		if (search=="genes") fit_with = LEGIT(data=data, genes=cbind(genes_current,genes_toadd[,j,drop=FALSE]), env=env_current, formula=formula, start_genes=c(start_genes,0), start_env=start_env, eps=eps, maxiter=maxiter, family=family, print=FALSE)
		if (search=="env") fit_with = LEGIT(data=data, genes=genes_current, env=cbind(env_current,env_toadd[,j,drop=FALSE]), formula=formula, start_genes=start_genes, start_env=c(start_env,0), eps=eps, maxiter=maxiter, family=family, print=FALSE)
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
			fit_without = LEGIT(data=data[comp,,drop=FALSE], genes=genes_current[comp,,drop=FALSE], env=env_current[comp,,drop=FALSE], formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, print=FALSE)
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
		if (search_criterion=="AIC" || search_criterion=="BIC"){
			if (min(criterion_diff) > 0){
				if (search=="genes" && print) cat("No gene added\n")
				if (search=="env" && print) cat("No environment added\n")
				return(NULL)
			}
			if (empty_start_dataset) best_var = which.min(criterion_after)
			else best_var = which.min(criterion_diff)
			if (search=="genes"){
				if (print) cat(paste0("Adding gene: ",colnames(genes_toadd)[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after= ",round(criterion_after[best_var],5),")\n"))
				genes_current = cbind(genes_current, genes_toadd[,best_var, drop=FALSE])
				genes_toadd = genes_toadd[,-best_var, drop=FALSE]
			}
			if (search=="env"){
				if (print) cat(paste0("Adding environment: ",colnames(env_toadd)[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after= ",round(criterion_after[best_var],5),")\n"))
				env_current = cbind(env_current, env_toadd[,best_var, drop=FALSE])
				env_toadd = env_toadd[,-best_var, drop=FALSE]						
			}
		}
	}
	# Cross-validated models
	if (search_criterion == "cv" || search_criterion == "cv_AUC"){
		# Dropping variables with p < threshold and worse AIC
		if (interactive_mode) interactive = interactive[good,]
		if (search=="genes") genes_toadd = genes_toadd[,good, drop=FALSE]
		if (search=="env") env_toadd = env_toadd[,good, drop=FALSE]
		if (search=="genes") elements_N = NCOL(genes_toadd)
		if (search=="env") elements_N = NCOL(env_toadd)
		# Vector which says how much the criterion changed from including the variable
		criterion_before = rep(0, elements_N)
		criterion_after = rep(0, elements_N)
		criterion_diff = rep(0, elements_N)
		# Set seed
		if (!is.null(seed)) current_seed = seed
		else current_seed = NULL
		if (!empty_start_dataset) fit_cv = LEGIT_cv(data=data, genes=genes_current, env=env_current, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=current_seed)
		for (j in 1:elements_N){
			if (search=="genes") fit_cv_with = LEGIT_cv(data=data, genes=cbind(genes_current,genes_toadd[,j,drop=FALSE]), env=env_current, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=c(start_genes,0), start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=current_seed)
			if (search=="env") fit_cv_with = LEGIT_cv(data=data, genes=genes_current, env=cbind(env_current,env_toadd[,j,drop=FALSE]), formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=start_genes, start_env=c(start_env,0), eps=eps, maxiter=maxiter, family=family, seed=current_seed)
			if (empty_start_dataset) fit_cv_without = NULL
			else{
				if (search=="genes") comp_with = stats::complete.cases(data,genes_current,genes_toadd[,j,drop=FALSE],env_current)
				if (search=="env") comp_with = stats::complete.cases(data,genes_current,env_toadd[,j,drop=FALSE],env_current)
				if (sum(comp_without) != sum(comp_with)) fit_cv_without = LEGIT_cv(data=data[comp_with,,drop=FALSE], genes=genes_current[comp_with,,drop=FALSE], env=env_current[comp_with,,drop=FALSE], formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=current_seed)
				else fit_cv_without = fit_cv
			}
			if (search_criterion=="cv"){
				if (empty_start_dataset) criterion_before[j] = 0
				else criterion_before[j] = mean(fit_cv_without$R2_cv)
				criterion_after[j] = mean(fit_cv_with$R2_cv)
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
				if(classification){
					if (empty_start_dataset) interactive$cv_AUC_old[j] = 0
					else interactive$cv_AUC_old[j] = mean(fit_cv_without$AUC)
					interactive$cv_AUC_new[j] = mean(fit_cv_with$AUC)
				}
			}
		}
		# Only do this if NOT in interactive mode, otherwise at the end of the algorithm, we show the the data.frame interactive and let the user choose
		if (!interactive_mode){
			if (max(criterion_diff) < 0){
				if (search=="genes" && print) cat("No gene added\n")
				if (search=="env" && print) cat("No environment added\n")
				return(NULL)
			}
			best_var = which.max(criterion_diff)
			if (search=="genes"){
				if (print) cat(paste0("Adding gene: ",colnames(genes_toadd)[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after= ",round(criterion_after[best_var],5),")\n"))
				genes_current = cbind(genes_current, genes_toadd[,best_var, drop=FALSE])
				genes_toadd = genes_toadd[,-best_var, drop=FALSE]
			}
			if (search=="env"){
				if (print) cat(paste0("Adding environment: ",colnames(env_toadd)[best_var], " (Criterion before = ",round(criterion_before[best_var],5), "; after= ",round(criterion_after[best_var],5),")\n"))
				env_current = cbind(env_current, env_toadd[,best_var, drop=FALSE])
				env_toadd = env_toadd[,-best_var, drop=FALSE]				
			}
		}
	}
	if (interactive_mode){
		interactive_n = min(5,elements_N)
		if (search_criterion=="AIC" && empty_start_dataset) neworder = order(interactive$AIC_new, decreasing=FALSE)
		if (search_criterion=="AIC" && !empty_start_dataset) neworder = order(interactive$AIC_new - interactive$AIC_old, decreasing=FALSE)
		if (search_criterion=="BIC" && empty_start_dataset) neworder = order(interactive$BIC_new, decreasing=FALSE)
		if (search_criterion=="BIC" && !empty_start_dataset) neworder = order(interactive$BIC_new - interactive$BIC_old, decreasing=FALSE)
		if (search_criterion=="cv" && empty_start_dataset) neworder = order(interactive$cv_R2_new, decreasing=TRUE)
		if (search_criterion=="cv" && !empty_start_dataset) neworder = order(interactive$cv_R2_new - interactive$cv_R2_old, decreasing=TRUE)
		if (search_criterion=="cv_AUC" && empty_start_dataset) neworder = order(interactive$cv_AUC_new, decreasing=TRUE)
		if (search_criterion=="cv_AUC" && !empty_start_dataset) neworder = order(interactive$cv_AUC_new - interactive$cv_AUC_old, decreasing=TRUE)
		interactive = interactive[neworder[1:interactive_n],]
		rownames(interactive)=1:interactive_n
		interactive = format(interactive,scientific=FALSE)
		print(interactive)
		input_user = readline(prompt="Enter the index of the variable to be added: ")
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
	fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, print=FALSE)
	start_genes = stats::coef(fit$fit_genes)
	start_env = stats::coef(fit$fit_env)
	if (search=="genes") return(list(fit=fit, start_genes=start_genes,start_env=start_env,genes_current=genes_current,genes_toadd=genes_toadd))
	if (search=="env") return(list(fit=fit, start_genes=start_genes,start_env=start_env,env_current=env_current,env_toadd=env_toadd))
}


backward_step = function(fit, data, formula, interactive_mode=FALSE, genes_current=NULL, env_current=NULL, genes_dropped=NULL, env_dropped=NULL, search="genes", search_criterion="AIC", p_threshold = .20, exclude_worse_AIC=TRUE, max_steps = 50, cv_iter=5, cv_folds=10, classification=FALSE, start_genes=NULL, start_env=NULL, eps=.01, maxiter=25, family=gaussian, seed=NULL, print=TRUE){
	# How much genes or env to add
	if (search=="genes") elements_N = NCOL(genes_current)
	if (search=="env") elements_N = NCOL(env_current)
	if (elements_N == 0) return(NULL)
	# Vector which says which variables are "good" (worth keeping)
	good = rep(FALSE, elements_N)
	# Vector which says how much the criterion changed from excluding the variable
	criterion_before = rep(0, elements_N)
	criterion_after = rep(0, elements_N)
	criterion_diff = rep(0, elements_N)
	# In interactive model, we must keep track of every AIC, BIC, p-value, Cross-validated R2 and AUC to show the user at every iteration
	if (interactive_mode){
		if (search=="genes") interactive = data.frame(variable=colnames(genes_current), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N),cv_R2_old=rep(NA, elements_N),cv_R2_new=rep(NA, elements_N),cv_AUC_old=rep(NA, elements_N),cv_AUC_new=rep(NA, elements_N))
		if (search=="env") interactive = data.frame(variable=colnames(env_current), N_old=rep(NA, elements_N), N_new=rep(NA, elements_N), p_value=rep(NA, elements_N),AIC_old=rep(NA, elements_N),AIC_new=rep(NA, elements_N),BIC_old=rep(NA, elements_N),BIC_new=rep(NA, elements_N),cv_R2_old=rep(NA, elements_N),cv_R2_new=rep(NA, elements_N),cv_AUC_old=rep(NA, elements_N),cv_AUC_new=rep(NA, elements_N))
	}
	# we need to always use this sample as it could be different when removing a variable
	comp_with = stats::complete.cases(data,genes_current,env_current)
	fit_with = fit
	p_value = stats::coef(summary(fit_with)$fit_genes)[,4]
	# Non-cross-validated models
	for (j in 1:elements_N){
		if (search=="genes"){
			comp_without = stats::complete.cases(data,genes_current[,-j,drop=FALSE],env_current)
			fit_without = LEGIT(data=data[comp_with,,drop=FALSE], genes=genes_current[comp_with,-j,drop=FALSE], env=env_current[comp_with,,drop=FALSE], formula=formula, start_genes=start_genes[-j], start_env=start_env, eps=eps, maxiter=maxiter, family=family, print=FALSE)
		}
		if (search=="env"){
			comp_without = stats::complete.cases(data,genes_current,env_current[,-j,drop=FALSE])
			fit_without = LEGIT(data=data[comp_with,,drop=FALSE], genes=genes_current[comp_with,,drop=FALSE], env=env_current[comp_with,-j,drop=FALSE], formula=formula, start_genes=start_genes, start_env=start_env[-j], eps=eps, maxiter=maxiter, family=family, print=FALSE)
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
		if (search_criterion=="AIC" || search_criterion=="BIC"){
			if (min(criterion_diff) > 0){
				if (search=="genes" && print) cat("No gene removed\n")
				if (search=="env" && print) cat("No environment removed\n")
				return(NULL)
			}
			worst_var = which.min(criterion_diff)
			if (search=="genes"){
				if (print) cat(paste0("Removing gene: ",colnames(genes_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after= ",round(criterion_after[worst_var],5),")\n"))
				genes_dropped = cbind(genes_dropped, genes_current[,worst_var, drop=FALSE])
				genes_current = genes_current[,-worst_var, drop=FALSE]
			}
			if (search=="env"){
				if (print) cat(paste0("Removing environment: ",colnames(env_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after= ",round(criterion_after[worst_var],5),")\n"))
				env_dropped = cbind(env_dropped, env_current[,worst_var, drop=FALSE])
				env_current = env_current[,-worst_var, drop=FALSE]
			}
		}
	}
	# Cross-validated models
	if (search_criterion == "cv" || search_criterion == "cv_AUC"){
		# Not looking at variables that labelled as good
		elements_N_cv = sum(!good)
		# Vector which says how much the criterion changed from removing the variable
		criterion_before = rep(0, elements_N)
		criterion_after = rep(0, elements_N)
		criterion_diff = rep(0, elements_N)
		# Set seed
		if (!is.null(seed)) current_seed = seed
		else current_seed = NULL
		fit_cv_with = LEGIT_cv(data=data, genes=genes_current, env=env_current, formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=current_seed)
		for (j in 1:elements_N){
			# Only do this if not labelled as good
			if (!good[j]){
				if (search=="genes") comp_without = stats::complete.cases(data,genes_current[,-j,drop=FALSE],env_current)
				if (search=="env") comp_without = stats::complete.cases(data,genes_current,env_current[,-j,drop=FALSE])
				if (search=="genes") fit_cv_without = LEGIT_cv(data=data[comp_with,,drop=FALSE], genes=genes_current[comp_with,-j,drop=FALSE], env=env_current[comp_with,,drop=FALSE], formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=start_genes[-j], start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=current_seed)
				if (search=="env") fit_cv_without = LEGIT_cv(data=data[comp_with,,drop=FALSE], genes=genes_current[comp_with,,drop=FALSE], env=env_current[comp_with,-j,drop=FALSE], formula=formula, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=start_genes, start_env=start_env[-j], eps=eps, maxiter=maxiter, family=family, seed=current_seed)
				if (search_criterion=="cv"){
					criterion_before[j] = mean(fit_cv_with$R2_cv)
					criterion_after[j] = mean(fit_cv_without$R2_cv)
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
					if(classification){
						interactive$cv_AUC_old[j] = mean(fit_cv_with$AUC)
						interactive$cv_AUC_new[j] = mean(fit_cv_without$AUC)
					}
				}
			}
		}
		# Only do this if NOT in interactive mode, otherwise at the end of the algorithm, we show the the data.frame interactive and let the user choose
		if (!interactive_mode){
			if (max(criterion_diff) < 0){
				if (search=="genes" && print) cat("No gene removed\n")
				if (search=="env" && print) cat("No environment removed\n")
				return(NULL)
			}
			worst_var = which.max(criterion_diff)
			if (search=="genes"){
				if (print) cat(paste0("Removing gene: ",colnames(genes_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after= ",round(criterion_after[worst_var],5),")\n"))
				genes_dropped = cbind(genes_dropped, genes_current[,worst_var, drop=FALSE])
				genes_current = genes_current[,-worst_var, drop=FALSE]
			}
			if (search=="env"){
				if (print) cat(paste0("Removing environment: ",colnames(env_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after= ",round(criterion_after[worst_var],5),")\n"))
				env_dropped = cbind(env_dropped, env_current[,worst_var, drop=FALSE])
				env_current = env_current[,-worst_var, drop=FALSE]
			}
		}
	}
	if (interactive_mode){
		if (search_criterion=="cv") interactive_n = min(5,elements_N_cv)
		else interactive_n = min(5,elements_N)
		if (search_criterion=="AIC") neworder = order(interactive$AIC_new - interactive$AIC_old, decreasing=FALSE)
		if (search_criterion=="BIC") neworder = order(interactive$BIC_new - interactive$BIC_old, decreasing=FALSE)
		if (search_criterion=="cv") neworder = order(interactive$cv_R2_new - interactive$cv_R2_old, decreasing=TRUE)
		if (search_criterion=="cv_AUC") neworder = order(interactive$cv_AUC_new - interactive$cv_AUC_old, decreasing=TRUE)
		interactive = interactive[neworder[1:interactive_n],]
		rownames(interactive)=1:interactive_n
		interactive = format(interactive,scientific=FALSE)
		print(interactive)
		input_user = readline(prompt="Enter the index of the variable to be removed: ")
		if (sum(input_user == rownames(interactive))==0){
			if (search=="genes" && print) cat("No gene removed\n")
			if (search=="env" && print) cat("No environment removed\n")
			return(NULL)
		}
		worst_var = neworder[1:interactive_n][input_user == rownames(interactive)]
		if (search=="genes"){
			if (print) cat(paste0("Removing gene: ",colnames(genes_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after= ",round(criterion_after[worst_var],5),")\n"))
			genes_dropped = cbind(genes_dropped, genes_current[,worst_var, drop=FALSE])
			genes_current = genes_current[,-worst_var, drop=FALSE]
		}
		if (search=="env"){
			if (print) cat(paste0("Removing environment: ",colnames(env_current)[worst_var], " (Criterion before = ",round(criterion_before[worst_var],5), "; after= ",round(criterion_after[worst_var],5),")\n"))
			env_dropped = cbind(env_dropped, env_current[,worst_var, drop=FALSE])
			env_current = env_current[,-worst_var, drop=FALSE]
		}
	}
	# Updated model and coefficients
	if (search=="genes") start_genes=start_genes[-worst_var]
	if (search=="env") start_env=start_env[-worst_var]
	fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, print=FALSE)
	start_genes = stats::coef(fit$fit_genes)
	start_env = stats::coef(fit$fit_env)
	if (search=="genes") return(list(fit=fit, start_genes=start_genes,start_env=start_env,genes_current=genes_current,genes_dropped=genes_dropped))
	if (search=="env") return(list(fit=fit, start_genes=start_genes,start_env=start_env,env_current=env_current,env_dropped=env_dropped))
}


stepwise_search = function(data, formula, interactive_mode=FALSE, genes_original=NULL, env_original=NULL, genes_extra=NULL, env_extra=NULL, search_type="bidirectional-forward", search="genes", search_criterion="AIC", forward_exclude_p_bigger = .20, backward_exclude_p_smaller = .01, exclude_worse_AIC=TRUE, max_steps = 50, cv_iter=5, cv_folds=10, classification=FALSE, start_genes=NULL, start_env=NULL, eps=.01, maxiter=25, family=gaussian, seed=NULL, print=TRUE){
	if (forward_exclude_p_bigger > 1 || forward_exclude_p_bigger <= 0) stop("forward_exclude_p_bigger must be between 0 and 1 (Set to 1 to ignore p-values in forward step)")
	if (backward_exclude_p_smaller >= 1 || backward_exclude_p_smaller < 0) stop("backward_exclude_p_smaller must be between 0 and 1 (Set to 0 to ignore p-values in backward step)")
	if (search_criterion=="AIC") string_choice="lowest AIC"
	else if (search_criterion=="BIC") string_choice="lowest BIC"
	else if (search_criterion=="cv") string_choice="lowest cross-validation error"
	else if (search_criterion=="cv_AUC") string_choice="biggest cross-validated area under the curve"
	else stop("Not a valid search_criterion, use: AIC, BIC, cv or cv_AUC")

	comp = stats::complete.cases(data, genes_original, env_original, genes_extra, env_extra)
	if (sum(comp)!=length(comp) && print) cat("Note: Missing data was found. This will increase computing time in forward steps \n and could lead to an incorrect subset of variable if the sample size change too much. \n")

	if (interactive_mode && print) cat("<<~ Interative mode enabled ~>>\n")

	# If true, then we start with no genes or envand must find the best one first
	if ((is.null(genes_original) && search=="genes") || (is.null(env_original) && search=="env")) empty_start_dataset = TRUE
	else empty_start_dataset = FALSE

	# Setting up initial weighted scores
	if (is.null(genes_original) && search=="genes") start_genes = c()
	else if (is.null(start_genes)) start_genes = rep(1/dim(genes_original)[2],dim(genes_original)[2])

	if (is.null(env_original) && search=="env") start_env = c()
	else if (is.null(start_env)) start_env = rep(1/dim(env_original)[2],dim(env_original)[2])

	fit = NULL

	if (search_type == "forward"){
		if (print){
			if (forward_exclude_p_bigger < 1 && (exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Keeping only variables with p-values smaller than ", forward_exclude_p_bigger, " and which inclusion decrease the AIC\n"))
			else if (forward_exclude_p_bigger < 1 && !(exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Keeping only variables with p-values smaller than ", forward_exclude_p_bigger,"\n"))
			else if (exclude_worse_AIC || search_criterion=="AIC") cat(paste0("Keeping only variables which inclusion decrease the AIC\n"))
			if ((search_criterion=="cv" || search_criterion=="cv_AUC") && (!exclude_worse_AIC || forward_exclude_p_bigger > .2)) cat("Note : We recommend using exclude_worse_AIC=TRUE and forward_exclude_p_bigger <= .20 to reduce the amount of cross-validations needed and to make the algorithm much faster.\n")
			if (search_criterion=="cv_AUC") classification=TRUE
			if (search=="genes") cat(paste0("Forward search of the genes to find the model with the ", string_choice,"\n"))
			if (search=="env") cat(paste0("Forward search of the environments to find the model with the ", string_choice,"\n"))
		}
		genes_current = genes_original
		env_current = env_original
		if (!empty_start_dataset){
			# Original model
			fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, print=FALSE)
			start_genes = stats::coef(fit$fit_genes)
			start_env = stats::coef(fit$fit_env)
		}
		else{
			# If dataset is NULL, then make the current dataset be an empty data frame of same size as the dataset with extra variables 
			if (is.null(genes_original)) genes_current = genes_extra[,-c(1:NCOL(genes_extra))]
			if (is.null(env_original)) env_current = env_extra[,-c(1:NCOL(env_extra))]
		}
		if (search=="genes") genes_toadd = genes_extra
		else genes_toadd = NULL
		if (search=="env") env_toadd = env_extra
		else env_toadd = NULL
		if (print) cat("\n")
		for (i in 1:max_steps){
			if (print) cat(paste0("[Iteration: ",i,"]\n"))
			results = forward_step(empty_start_dataset=empty_start_dataset, fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, genes_current=genes_current, env_current=env_current, genes_toadd=genes_toadd, env_toadd=env_toadd, search=search, search_criterion=search_criterion, p_threshold = forward_exclude_p_bigger, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
			if (is.null(results)) break
			# Resetting parameters based on iteration results
			empty_start_dataset = FALSE
			fit = results$fit
			start_genes=results$start_genes
			start_env=results$start_env
			if (search=="genes"){
				genes_current=results$genes_current
				genes_toadd=results$genes_toadd
			}
			if (search=="env"){
				env_current=results$env_current
				env_toadd=results$env_toadd
			}
		}
	}
	if (search_type == "backward"){
		if (print){
			if (backward_exclude_p_smaller < 1 && (exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Dropping only variables with p-values bigger than ", backward_exclude_p_smaller, " and which removal decrease the AIC\n"))
			else if (backward_exclude_p_smaller < 1 && !(exclude_worse_AIC || search_criterion=="AIC")) cat(paste0("Dropping only variables with p-values bigger than ", backward_exclude_p_smaller,"\n"))
			else if (exclude_worse_AIC || search_criterion=="AIC") cat(paste0("Dropping only variables which removal decrease the AIC\n"))
			if ((search_criterion=="cv" || search_criterion=="cv_AUC") && (!exclude_worse_AIC || backward_exclude_p_smaller < .01)) cat("Note : We recommend using exclude_worse_AIC=TRUE and backward_exclude_p_smaller >= .01 to reduce the amount of cross-validations needed and to make the algorithm much faster.\n")
			if (search_criterion=="cv_AUC") classification=TRUE
			if (search=="genes") cat(paste0("Backward search of the genes to find the model with the ", string_choice,"\n"))
			if (search=="env") cat(paste0("Backward search of the environments to find the model with the ", string_choice,"\n"))
		}
		genes_current = genes_original
		env_current = env_original
		# Original model
		fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, print=FALSE)
		start_genes = stats::coef(fit$fit_genes)
		start_env = stats::coef(fit$fit_env)
		# Creating empty data frames of same size as the other data frames
		genes_dropped = genes_original[,-c(1:NCOL(genes_original))]
		env_dropped = env_original[,-c(1:NCOL(env_original))]
		if (print) cat("\n")
		for (i in 1:max_steps){
			if (print) cat(paste0("[Iteration: ",i,"]\n"))
			results = backward_step(fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, genes_current=genes_current, genes_dropped=genes_dropped, env_current=env_current, env_dropped=env_dropped, search=search, search_criterion=search_criterion, p_threshold = backward_exclude_p_smaller, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
			if (is.null(results)) break
			# Resetting parameters based on iteration results
			fit = results$fit
			start_genes=results$start_genes
			start_env=results$start_env
			if (search=="genes"){
				genes_current=results$genes_current
				if (NCOL(genes_current)==1){
					if (print) cat("Only one gene remaining, stopping the algorithm\n")
					break
				}
				# Only for bidirectional
				genes_dropped=results$genes_dropped
			}
			if (search=="env"){
				env_current=results$env_current
				if (NCOL(env_current)==1){
					if (print) cat("Only one environment remaining, stopping the algorithm\n")
					break
				}
				# Only for bidirectional
				env_dropped=results$env_dropped
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
			if ((search_criterion=="cv" || search_criterion=="cv_AUC") && (!exclude_worse_AIC || forward_exclude_p_bigger > .2)) cat("Note : We recommend using exclude_worse_AIC=TRUE and forward_exclude_p_bigger <= .20 to reduce the amount of cross-validations needed and to make the algorithm much faster.\n")
			if (search_criterion=="cv_AUC") classification=TRUE
			if (search=="genes") cat(paste0("Bidirectional search of the genes to find the model with the ", string_choice,"\n"))
			if (search=="env") cat(paste0("Bidirectional search of the environments to find the model with the ", string_choice,"\n"))
		}
		if (search_type == "bidirectional-forward"){
			genes_current = genes_original
			env_current = env_original
			if (!empty_start_dataset){
				# Original model
				fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, print=FALSE)
				start_genes = stats::coef(fit$fit_genes)
				start_env = stats::coef(fit$fit_env)
			}
			else{
				# If dataset is NULL, then make the current dataset be an empty data frame of same size as the dataset with extra variables 
				if (is.null(genes_original)) genes_current = genes_extra[,-c(1:NCOL(genes_extra))]
				if (is.null(env_original)) env_current = env_extra[,-c(1:NCOL(env_extra))]
			}
			if (search=="genes") genes_toadd_ordrop = genes_extra
			else genes_toadd_ordrop = NULL
			if (search=="env") env_toadd_ordrop = env_extra
			else env_toadd_ordrop = NULL
			direction="forward"
			forward_failed=FALSE
			# Count as failed because we can't backward if forward doesn't work
			backward_failed=TRUE
		}
		if (search_type == "bidirectional-backward"){
			genes_current = genes_original
			env_current = env_original
			# Original model
			fit = LEGIT(data=data, genes=genes_current, env=env_current, formula=formula, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, print=FALSE)
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
		for (i in 1:max_steps){
			if (print) cat(paste0("[Iteration: ",i,"]\n"))
			if (direction=="forward"){
				results = forward_step(empty_start_dataset=empty_start_dataset, fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, genes_current=genes_current, env_current=env_current, genes_toadd=genes_toadd_ordrop, env_toadd=env_toadd_ordrop, search=search, search_criterion=search_criterion, p_threshold = forward_exclude_p_bigger, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
				if (is.null(results)) forward_failed=TRUE
				else{
					forward_failed=FALSE
					# Resetting parameters based on iteration results
					empty_start_dataset = FALSE
					fit = results$fit
					start_genes=results$start_genes
					start_env=results$start_env
					if (search=="genes"){
						genes_current=results$genes_current
						genes_toadd_ordrop=results$genes_toadd
						if (NCOL(genes_toadd_ordrop)==0){
							# Count as a fail
							forward_failed=TRUE
						}
					}
					if (search=="env"){
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
			if (forward_failed && backward_failed) break
			if (direction=="backward"){
				# Can't backward if only one remaining variable
				if (!((NCOL(genes_current)<=1 && search=="genes") || (NCOL(env_current)<=1 && search=="env"))){
					results = backward_step(fit=fit, data=data, formula=formula, interactive_mode=interactive_mode, genes_current=genes_current, genes_dropped=genes_toadd_ordrop, env_current=env_current, env_dropped=env_toadd_ordrop, search=search, search_criterion=search_criterion, p_threshold = backward_exclude_p_smaller, exclude_worse_AIC=exclude_worse_AIC, max_steps = max_steps, cv_iter=cv_iter, cv_folds=cv_folds, classification=classification, start_genes=start_genes, start_env=start_env, eps=eps, maxiter=maxiter, family=family, seed=seed, print=print)
					if (is.null(results)) backward_failed=TRUE
					else{
						backward_failed=FALSE
						# Resetting parameters based on iteration results
						fit = results$fit
						start_genes=results$start_genes
						start_env=results$start_env
						one_remain=FALSE
						if (search=="genes"){
							genes_current=results$genes_current
							if (NCOL(genes_current)==1){
								# Count as a fail
								backward_failed=TRUE
							}
							genes_toadd_ordrop=results$genes_dropped
						}
						if (search=="env"){
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
					if (search=="genes" && print) cat("No gene removed1\n")
					if (search=="env" && print) cat("No environment removed1\n")
				}
				direction="forward"
			}
			if (forward_failed && backward_failed) break
		}
	}
	return(fit)
}