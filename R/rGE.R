#' @title Gene-Environment correlation estimation and testing
#' @description Estimates the gene-environment correlation (rGE) and tests for a GxE using a residual environmental score. If there is an important correlation between G and E, the model is still valid prediction-wise but the interpretation is affected as the question becomes: is it really a GxE or a GxG since E is partially caused by G? To account for this, we remove the influence of G on E (If E = b0 + b1*G + e, we use E_resid = E - b1*G) and refit the model to see if the model parameters changed. The residual environmental score (E_resid) is uncorrelated with G. This does not account for passive rGE but only active rGE.
#' @param object An object of class "LEGIT" or "IMLEGIT".
#' @param ... Further arguments passed to or from other methods.
#' @export
"rGE"

#' @title Gene-Environment correlation estimation and testing of LEGIT models
#' @description Estimates the gene-environment correlation (rGE) and tests for a GxE using a residual environmental score. If there is an important correlation between G and E, the model is still valid prediction-wise but the interpretation is affected as the question becomes: is it really a GxE or a GxG since E is partially caused by G? To account for this, we remove the influence of G on E (If E = b0 + b1*G + e, we use E_resid = E - b1*G) and refit the model to see if the model parameters changed. The residual environmental score (E_resid) is uncorrelated with G. This does not account for passive rGE but only active rGE.
#' @param object An object of class "LEGIT", usually, a result of a call to LEGIT.
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param ... Further arguments passed to or from other methods.
#' @return Returns a list containing the Pearson correlation and Kendall tau correlation of G and E and a glm fit of the main model part when removing the influence of G on E so that E and G are now uncorrelated.
#' @examples
#'	# Note: These examples don't have G and E correlation so the model fit doesn't change
#'	# but this shows how to use the rGE function
#'	train = example_2way(500, 1, seed=777)
#'	fit = LEGIT(train$data, train$G, train$E, y ~ G*E)
#'	fit_rGE = rGE(fit, y ~ G*E)
#'	fit_rGE
#'	summary(fit_rGE$fit_main_resid)
#' @export
"rGE.LEGIT"

#' @title Gene-Environment correlation estimation and testing of IMLEGIT models
#' @description Estimates the gene-environment correlation (rGE) and tests for a GxE using a residual environmental score. If there is an important correlation between G and E, the model is still valid prediction-wise but the interpretation is affected as the question becomes: is it really a GxE or a GxG since E is partially caused by G? To account for this, we remove the influence of G on E (If E = b0 + b1*G + e, we use E_resid = E - b1*G) and refit the model to see if the model parameters changed. The residual environmental score (E_resid) is uncorrelated with G. This does not account for passive rGE but only active rGE.
#' @param object An object of class "IMLEGIT", usually, a result of a call to IMLEGIT.
#' @param formula Model formula. The names of \code{latent_var} can be used in the formula to represent the latent variables. If names(\code{latent_var}) is NULL, then L1, L2, ... can be used in formula to represent the latent variables. Do not manually code interactions, write them in the formula instead (ex: G*E1*E2 or G:E1:E2).
#' @param latent_var list of data.frame. The elements of the list are the datasets used to construct each latent variable. For interpretability and proper convergence, not using the same variable in more than one latent variable is highly recommended. It is recommended to set names to the list elements to prevent confusion because otherwise the latent variables will be named L1, L2, ... (See examples below for more details)
#' @param index_E vector or scalar representing the index of each latent variable that is part of the "environment"
#' @param index_G scalar representing the index of the latent variable for the "genetic" part
#' @param ... Further arguments passed to or from other methods.
#' @return Returns a list containing the Pearson correlation and Kendall tau correlation of G and E and a glm fit of the main model part when removing the influence of G on E so that E and G are now uncorrelated.
#' @examples
#'	# Note: These examples don't have G and E correlation so the model fit doesn't change 
#'	# but this shows how to use the rGE function
#'	train = example_3way_3latent(500, 1, seed=777)
#'	fit = IMLEGIT(train$data, train$latent_var, y ~ G*E*Z)
#'	# If we assume Z not to be an "environment"
#'	fit_rGE1 = rGE(fit, y ~ G*E, train$latent_var, 2, 1)
#'	fit_rGE1
#'	summary(fit_rGE1$fit_main_resid)
#'	# If we assume Z to be an "environment"
#'	fit_rGE2 = rGE(fit, y ~ G*E, train$latent_var, c(2,3), 1)
#'	fit_rGE2
#'	summary(fit_rGE2$fit_main_resid)
#' @export
"rGE.IMLEGIT"

rGE = function(object, ...) UseMethod("rGE")

rGE.LEGIT = function(object, formula, ...){
	fit_E = stats::lm(object$fit_main$data$E ~ object$fit_main$data$G)
	# replacing with residual E
	object$fit_main$data$E = (object$fit_main$data$E - coef(fit_E)[2]*object$fit_main$data$G)
	fit_resid = stats::glm(formula, data=object$fit_main$data, family=object$fit_main$family, y=FALSE, model=FALSE)
	#Changed AIC
	fit_resid$aic = fit_resid$aic + 2*(object$true_model_parameters$rank - fit_resid$rank)
	return(list(rGE_pearson = Hmisc::rcorr(as.matrix(object$fit_main$data[,c("G","E")]),type="pearson"), rGE_kendall = Hmisc::rcorr(as.matrix(object$fit_main$data[,c("G","E")]),type="spearman"), fit_main_resid=fit_resid))
}

rGE.IMLEGIT = function(object, formula, latent_var, index_E, index_G, ...){
	if (length(index_G) > 1) stop("index_G must be a single latent variable")
	if (sum(index_G > length(latent_var)) > 0 || sum(index_G > length(latent_var))> 0) stop("elements of index_G and index_E must be smaller than the number of latent_variable")

	for (i in 1:length(index_E)){
		fit_E = stats::lm(object$fit_main$data[,names(latent_var)[index_E[1]]] ~ object$fit_main$data[,names(latent_var)[index_G]])
		# replacing with residual E
		object$fit_main$data[,names(latent_var)[index_E[1]]] = (object$fit_main$data[,names(latent_var)[index_E[1]]] - coef(fit_E)[2]*object$fit_main$data[,names(latent_var)[index_G]])
	}
	fit_resid = stats::glm(formula, data=object$fit_main$data, family=object$fit_main$family, y=FALSE, model=FALSE)
	#Changed AIC
	fit_resid$aic = fit_resid$aic + 2*(object$true_model_parameters$rank - fit_resid$rank)
	return(list(rGE_pearson = Hmisc::rcorr(as.matrix(object$fit_main$data[,names(latent_var)]),type="pearson"), rGE_kendall = Hmisc::rcorr(as.matrix(object$fit_main$data[,names(latent_var)]),type="spearman"), fit_main_resid=fit_resid))
}