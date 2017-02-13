*--------------------------------------------------------------------------------*
*           Latent Environmental & Genetic InTeraction (LEGIT) Modelling   	     *
*--------------------------------------------------------------------------------*
* Author: Alexia Jolicoeur-Martineau                                             *
* Email : alexia.jolicoeur-martineau@mail.mcgill.ca                      		 *
*                                                                                *
* Created: 30 December 2016                                              		 *
* Revised: 13 February 2017                                              		 *
*                                                                                *
* Version: 1.0.0 																 *	
*                                                                                *	
* Macros:                                                                        *
*	1) LEGIT_glimmix: Constructs a generalized linear mixed model 				 *
*		(using PROC GLIMMIX) with a weighted latent environmental score and 	 *
*		weighted latent genetic score.. Output model, convergence statistics     *
*		and pseudo AIC.      													 *
*	1) LEGIT_mixed: Constructs a linear mixed model (using PROC MIXED) with 	 *
*		a weighted latent environmental score and weighted latent genetic 		 *
* 		score. Output model, convergence statistics and AIC and approximated 	 *
* 		cross-validated R^2.										 			 *
*	1) LEGIT_logistic: Constructs a logistic model (using PROC LOGISTIC) with 	 *
*		a weighted latent environmental score and weighted latent genetic 		 *
* 		score. Output model, convergence statistics and AIC and approximated 	 *
* 		cross-validated AUC.										 			 *
*	4) LEGIT_cv: Calculates the leave-one-out (removing individuals rather  	 *
*		than observations if repeated-measure model) cross-validated R^2. Also   *
*		shows the individuals with extreme residuals for outlier detection and   *
*		removal. 																 *
*   5) LEGIT_search: Add genes or environments, one at a time, to a GxE      	 *
*		model. Output models for which the added variable had a p-value smaller  *
* 		than the threshold and lower AIC. 										 *
*																				 *
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *
* R package vs SAS macro: The R package LEGIT provides similar functions 	 	 *
* but the implementation is incredibly faster. In R, it is also much easier and  *
* simpler to set the model properly with a single formula instead of in 4 parts. *
* The stepwise search function in R is also much better and does it for 		 *
* multiple steps (instead of just one step in SAS) and either does it 			 *
* automatically or let you enter an interactive mode where you choose which 	 *
* variable to add at every step based on the information shown to you.  		 *
* The only thing that the R package cannot do is mixed models. Therefore, for 	 *
* non-mixed models I recommend using R exclusively and for mixed models, I 		 *
* recommend doing stepwise search in R to find the best subset of variables 	 *
* assuming no random effects and then refitting the best model in SAS with the	 *
* random effects.																 * 
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *
*              																	 *
* Notes : 																	     *
*	- Although PROC GLIMMIX is more versatile, it can only provide               *
*	pseudo-likelihood (thus pseudo AIC) and provides no way to calculate quickly *
*	the cross-validated R^2 or AUC. It is thus recommended to use                *
*	LEGIT_mixed for linear mixed models with continuous outcome and use    	 	 *
*	LEGIT_logistic for logistic regression model.							 	 *
*	- The cross-validated R^2 and AUC approximated in macros LEGIT_mixed,    	 *
*	LEGIT_logistic and LEGIT_search are invalid, they assume that the  			 *
*	genetic score and environmental score are known. We only provide this        *
*	measure to guide the user toward the models with high potential              *
*	out-of-sample effect sizes. Make sure to use macro #4 LEGIT_cv to        	 *
*	find the true leave-one-out crossvalidated R^2 or AUC before reporting your  *
*	results. 											 						 *
*	- Binary/Categorical outcomes are modelled in descending order.              *
*   - If you stop a macro while option no_log = 1, your log will be gone,        *
*	  to get it back do :	 													 *
*			proc printto; 														 *
*			run;																 *
*																				 *
* Lastest version on github.com/AlexiaJM										 *
*                                                                                *
* Created: using SAS 9.4 1M2 on windows                                          * 		
*--------------------------------------------------------------------------------*
*---------------------------------*
* Notes about model specification *                                                   
*---------------------------------*--------------------------------------------------*
* model_noG_noE : model part without environmental score E and genetic score G       *
*	 Ex : model_noG_noE = Intercept + x -> Intercept + x                             *
* model_E_noG : model part without genetic score G but with environmental score E    *
*	 Ex : model_E_noG = Intercept + x -> E + E*x                                     *
* model_G_noE : model part without environmental score E but with genetic score G    *
*	 Ex : model_G_noE = Intercept + x -> G + G*x                                     *
* model_G_E : model part with both environmental score E and genetic score G         *
*   Ex : model_G_E = Intercept + x -> G*E + G*E*x                                    *
*																					 *
* model_noG_noE + model_E_noG + model_G_noE + model_G_E must lead to the final model *
*																					 *
* For a two-way model E + G + E*G,                                                   *
*	 we set model_noG_noE = model_E_noG = model_G_noE = model_G_E = Intercept        *
* For a three-way model E + G + x + E*G + E*x + G*x + E*G*x ,                        *
*	 we set model_noG_noE = model_E_noG = model_G_noE = model_G_E = Intercept + x    *
*------------------------------------------------------------------------------------*

********************************************************************************************************;


*----------------------------------------------------------------------------*
* Macro #1) LEGIT_glimmix: Constructs a generalized linear mixed model 		 *
*		(using PROC GLIMMIX) with a weighted latent environmental score and  *
*		weighted latent genetic score.. Output model, convergence statistics *
*		and pseudo AIC.      												 *
*----------------------------------------------------------------------------*
* Parameters necessary to run *
*-----------------------------*---------------------------------------------------------------------*
* Data : dataset to be used                                                                         *
* outcome : outcome variable                                                                        *
* genes : genes variables inside genetic score G 												    *
*          (can be any sort of variable, doesn't even have to be genetic)						    *
* env : environments variables inside environmental score E 										*
*          (can be any sort of variable, doesn't even have to be environmental)					    *
* model_noG_noE : model part without E and without G (See model specification above) 			    *
* model_E_noG : model part with E and without G (See model specification above)                     *
* model_G_noE : model part without E and with G (See model specification above)                     *
* model_G_E : model part with E and with G (See model specification above)                          *
* id : ID of individual (if no ID, just create it using id=_n_ in a data step)                      *
*-------------------------------------------------*												    *
* Optional parameters or parameters with defaults *												    *
*-------------------------------------------------*-------------------------------------------------*
* covs (optional) : extra covariates (Equivalent to adding them in model_noG_noE)				    *
* remove_miss : if 1 then remove missing data before running (Default=1)						    *
* time (optional) : time variable for repeated measure outcome   								    *
* start_genes (optional) : starting points for genetic score (must be same length as "genes")       *
* start_env (optional) : starting points for environmental score (must be same length as "env")     *
* eps : threshold for convergence (.01 for quick batch simulations, .0001 for accurate results)     *
* maxiter : Maximum number of iterations														    *
* print : If 1 then print all models in all iterations, use for debugging only (Default=0)		    *
* print_final : If 1 then print final model (Default=1)											    *
* ods_new : If 1 then close current ods output and start a new one (Default=0)                      *
* repeated : If 1 then the outcome is a repeated measure (Default=0)                                *
* repeated_type : covariance type for repeated measure (Default= un)                                *
* random_vars (optional) : variables of random effects                                              *
* random_sub (optional) : subject for randon effects                                                *
* random_type : covariance type for random effects (Default= vc)                                    *
* where (optional) : where                                                                          *
* ods_clear : If 1 then clear ods cluttering garbage (Default=0)                                    *
* dist : Outcome distribution (Default= normal)                                                     *
* link : GLM link (Default=identity)                                                                *
* method : Optimization method (Default=MSPL)                                                       *
*---------------------------------------------------------------------------------------------------*;
*----------*
* Examples *
*----------*
* two-way model : Intercept + g + e + ge + gender + ses 																								
* Four genes g1, g2, g3, g4. Three environments e1, e2, e3. Two gene by gene interaction g1_g2, g1_g4. 
* Continuous outcome y at four time-points time="3M", "6M", "18M" or "36M".	
*
* LEGIT_glimmix(data, y, g1 g2 g3 g4 g1_g2 g1_g4, e1 e2 e3, Intercept, Intercept, Intercept, Intercept, gender ses, id=id, time=time, repeated=1)		
*
* three-way model with binary outcome : Intercept + g + e + z + ge + gz + ez + gez + gender + ses 														
* Four genes g1, g2, g3, g4. Three environments e1, e2, e3. Two gene by gene interaction g1_g2, g1_g4. 
* Binary outcome y at one time-point.
*
* LEGIT_glimmix(data, y, g1 g2 g3 g4 g1_g2 g1_g4, e1 e2 e3, Intercept z, Intercept z, Intercept z, Intercept z, gender ses, id=id, dist=binomial, link=logit)


%macro LEGIT_glimmix(data, outcome, genes, env, model_noG_noE=, model_E_noG=, model_G_noE=, model_G_E=, covs=, id=id, remove_miss=1, time=, start_genes=, start_env=, eps = .0001, maxiter = 50, print = 0, print_final=1, ods_new=1, repeated=0, repeated_type=un, random_vars=, random_sub=, random_type=vc, where=, clear_ods=0, dist=normal, link=identity, method=MSPL);

* Won't give an error if trying to drop G or E when it doesn't exist;
options dkricond=warn;

* Counting how many variables;
%let genes_N = %sysfunc(countw(&genes));
%let env_N = %sysfunc(countw(&env));
%let model_noG_noE_N = %sysfunc(countw(&model_noG_noE,' '));
%let model_E_noG_N = %sysfunc(countw(&model_E_noG,' '));
%let model_G_noE_N = %sysfunc(countw(&model_G_noE,' '));
%let model_G_E_N = %sysfunc(countw(&model_G_E,' '));
%let covs_N = %sysfunc(countw(&covs,' '));

%if &ods_new eq 1 %then %do;
	ods html close;
	ods html;
%end;

* Sorting, setting up intercept and outlier removal; 
Proc sort data = &data;
	by &id 
%if &repeated eq 1 %then %do;
	&time 
%end;
;
run;
data Step_a;
	set &data;
	Intercept = 1;
run;

%if &remove_miss eq 1 %then %do;
	data Step_a (keep = &id 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		&outcome &genes &env &model_noG_noE &model_E_noG &model_G_noE &model_G_E &covs &random_vars &random_sub);
		set Step_a;
			where &where;
	run;

	data Step_a;
		set Step_a;
		if cmiss(of _all_) then delete;
	run;
%end;

* Setting up datasets with initial weighted scores;
%if %bquote(&start_genes) eq %then %do;
	data G_weights;
		G_weights_old = 1/&genes_N;
		%do i=1 %to &genes_N;
			output;
		%end;
	run;
%end;
%else %do;
	proc iml;
		G_weights_old = {&start_genes};
		create G_weights var {G_weights_old};
		append;
	quit;
%end;
%if %bquote(&start_env) eq %then %do;
	data E_weights;
		E_weights_old = 1/&env_N;
		%do i=1 %to &env_N;
			output;
		%end;
	run;
%end;
%else %do;
	proc iml;
		E_weights_old = {&start_env};
		create E_weights var {E_weights_old};
		append;
	quit;
%end;

* Setting up for first step;
proc iml;
	use Step_a;
		read all var {&genes} into Genes[colname=varnames];
		read all var {&env} into Env[colname=varnames];
		read all var {&id};
	close Step_a;
	use G_weights;
		read all;
		G = Genes*G_weights_old/(sum(abs(G_weights_old)));
	close G_weights;
	use E_weights;
		read all;
		E = Env*E_weights_old/(sum(abs(E_weights_old)));
	close E_weights;
	create Step_a_GE var {&id G E};
	append;
	%if &print eq 1 %then %do;
		print G E;
	%end;
quit;
data Step_a;
	* Merging and dropping M and E just in case that they already exist in original dataset;
	merge Step_a_GE Step_a(drop = G E);
run;

%do i=1 %to &maxiter;

	%if &print ne 1 %then %do;
		ods exclude all;
	%end;

	* Step a : Estimating betas and covariance matrix;
	proc glimmix data=Step_a method=&method IC=PQ;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model &outcome(DESC) = &model_noG_noE 
			%do j=1 %to &model_E_noG_N;
				E*%scan(&model_E_noG, &j,' ') 
			%end;
			%do j=1 %to &model_G_noE_N;
				G*%scan(&model_G_noE, &j,' ') 
			%end;
			%do j=1 %to &model_G_E_N;
				E*G*%scan(&model_G_E, &j,' ') 
			%end;
		&covs / solution noint dist=&dist link=&link nocenter;
		%if &random_vars ne %then %do;
			random &random_vars  
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		%if &repeated eq 1 %then %do;
			random &time / subject = &id type = &repeated_type residual;
		%end;
		%else %do;
			random Intercept / subject=&id residual;
		%end;
		ods output ParameterEstimates=SolutionF CovParms=Cov;
	run;

	*Variable modification and setting up for b step (estimating G);
	proc iml;
		use Step_a;
			read all var {&genes} into Genes[colname=varnames];
			read all var {E};
			read all var {&model_noG_noE} into X_model_noG_noE[colname=varnames];
			read all var {&model_E_noG} into X_model_E_noG[colname=varnames];
			read all var {&model_G_noE} into X_model_G_noE[colname=varnames];
			read all var {&model_G_E} into X_model_G_E[colname=varnames];
			read all var {&covs} into X_covs[colname=varnames];
			read all var {&outcome};
		close Step_a;
		use SolutionF;
			read all var{Estimate};
		close SolutionF;
		R0 = X_model_noG_noE*Estimate[1:(&model_noG_noE_N)] + (X_model_E_noG#E)*Estimate[(&model_noG_noE_N+1):(&model_noG_noE_N+&model_E_noG_N)] + X_covs*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+&covs_N)];
		R1 = X_model_G_noE*Estimate[(&model_noG_noE_N+&model_E_noG_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N)] + (X_model_G_E#E)*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N)];
		R1_Genes = Genes#R1;
		R = &outcome||R0||R1_Genes;
		create Step_b from R[colname={&outcome R0 &genes}];
		append from R;		
		%if &print eq 1 %then %do;
			print R0 R1 R1_genes Y_new;
		%end;
	quit;
	data Step_b;
		merge Step_a(keep = &id &random_vars &random_sub  
			%if &repeated eq 1 %then %do;
				&time
			%end;
		) Step_b;
	run;

	*Step b: Estimating G weights;
	proc glimmix data=Step_b method=&method IC=PQ;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model &outcome(DESC) = &genes / noint solution dist=&dist link=&link offset=R0 nocenter;
		%if &repeated eq 1 %then %do;
			random &time / subject = &id type = un residual;
		%end;
		%else %do;
			random Intercept / subject=&id residual;
		%end;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		parms / noiter PDATA=Cov;
		ods output ParameterEstimates=SolutionF_G;
	run;

	*Updating G;
	proc iml;
		use Step_a;
			read all var {G};
			read all var {&id};
			read all var {&genes} into Genes[colname=varnames];
		close Step_a;
		G_old = G;
		use SolutionF_G;
			read all;
			use G_weights;
				read all var {G_weights_old};
				G = Genes*Estimate/(sum(abs(Estimate)));
				G_weights_new = Estimate/(sum(abs(Estimate)));
				* term to multiply estimates;
				total_inv_G = 1/(sum(abs(Estimate)));
				call symputx("total_inv_G", total_inv_G);
				* Convergence stuff;
				diff_G = sqrt(ssq(G_weights_new - G_weights_old));
				iter = num(symget("i"));
				diff_threshold = num(symget("eps"));
				if (diff_G < diff_threshold) then call symputx("conv_half1", 1);
				else call symputx("conv_half1", 0);
				* Updating estimates;
				G_weights_old = G_weights_new;
			close G_weights;
			create G_weights var {G_weights_old};
				append;
			create Step_b_G var {&id G};
				append;
			* To be reported at end of macro as output for users;
			create conv_status1 var {diff_G diff_threshold iter total};
				append;
		close SolutionF_G;
		%if &print eq 1 %then %do;
			print G_weights_new diff_G diff_threshold iter;
		%end;
	quit;
	data Step_a;
		merge Step_b_G Step_a(drop = G);
	run;

	*Variable modification and setting up for c step (estimating E);
	proc iml;
		use Step_a;
			read all var {&env} into Env[colname=varnames];
			read all var {G};
			read all var {&model_noG_noE} into X_model_noG_noE[colname=varnames];
			read all var {&model_E_noG} into X_model_E_noG[colname=varnames];
			read all var {&model_G_noE} into X_model_G_noE[colname=varnames];
			read all var {&model_G_E} into X_model_G_E[colname=varnames];
			read all var {&covs} into X_covs[colname=varnames];
			read all var {&outcome};
		close Step_a;
		use SolutionF;
			read all var{Estimate};
		close SolutionF;
		* Setting special parametrization for next E step;
		* Estimates are currently ordered as &model_noG_noE -> &model_E_noG -> &model_G_noE -> &model_G_E -> covs;
		* We need them in this order : &model_noG_noE -> &model_G_noE -> &model_E_noG -> &model_G_E -> covs;
		R0 = X_model_noG_noE*Estimate[1:(&model_noG_noE_N)] + (X_model_G_noE#G)*Estimate[(&model_noG_noE_N+&model_E_noG_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N)] + X_covs*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+&covs_N)];
		R1 = X_model_E_noG*Estimate[(&model_noG_noE_N+1):(&model_noG_noE_N+&model_E_noG_N)] + (X_model_G_E#G)*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N)];
		R1_Env = Env#R1;
		R = &outcome||R0||R1_Env;
		create Step_c from R[colname={&outcome R0 &env}];
		append from R;
		%if &print eq 1 %then %do;
			print R0 R1 R1_Env Y_new;
		%end;
	quit;
	data Step_c;
		merge Step_a(keep = &id &random_vars &random_sub 
		%if &repeated eq 1 %then %do;
			&time
		%end;
		) Step_c;
	run;
	*Step c: Estimating E weights;
	proc glimmix data=Step_c method=&method IC=PQ;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model &outcome(DESC) = &env / noint solution dist=&dist link=&link offset=R0 nocenter;
		%if &repeated eq 1 %then %do;
			random &time / subject = &id type = un residual;
		%end;
		%else %do;
			random Intercept / subject=&id residual;
		%end;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		parms / noiter PDATA=Cov;
		ods output ParameterEstimates=SolutionF_E;
	run;

	* Normalizing weights and combing everything into dataset for step a;
	proc iml;
		use Step_a;
			read all var {&env} into Env[colname=varnames];
			read all var {&id};
			read all var {E};
		close Step_a;
		E_old = E;
		use SolutionF_E;
			read all;
			use E_weights;
				read all var {E_weights_old};
				E = Env*Estimate/(sum(abs(Estimate)));
				E_weights_new = Estimate/(sum(abs(Estimate)));
				* term to multiply estimates;
				total_inv_E = 1/(sum(abs(Estimate)));
				call symputx("total_inv_E", total_inv_E);
				* Looking for convergence;
				diff_E = sqrt(ssq(E_weights_new - E_weights_old));
				iter = num(symget("i"));
				diff_threshold = num(symget("eps"));
				if (diff_E < diff_threshold) then call symputx("conv_half2", 1);
				else call symputx("conv_half2", 0);
				* Updating estimates;
				E_weights_old = E_weights_new;
			close E_weights;
			* Creating datasets;
			create E_weights var {E_weights_old};
				append;
			create Step_c_E var {&id E};
				append;
			close Step_c_E;
			create conv_status2 var {diff_E diff_threshold iter total};
				append;
		close SolutionF_E;
		%if &print eq 1 %then %do;
			print E_weights_new diff_E diff_threshold iter;
		%end;
	quit;
	data Step_a;
		merge Step_c_E Step_a(drop = E);
	run;
	* If both E and G converged then stop;
	%if (&conv_half1 eq 1 and &conv_half2 eq 1) %then %let i=&maxiter;
%end;

	%if &print_final eq 1 %then %do;
		ods exclude none;
	%end;

	*Rerunning step a for AIC and predictions;
	proc glimmix data=Step_a method=&method IC=PQ;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model &outcome(DESC) = &model_noG_noE 
			%do j=1 %to &model_E_noG_N;
				E*%scan(&model_E_noG, &j,' ') 
			%end;
			%do j=1 %to &model_G_noE_N;
				G*%scan(&model_G_noE, &j,' ') 
			%end;
			%do j=1 %to &model_G_E_N;
				E*G*%scan(&model_G_E, &j,' ') 
			%end;
		&covs / solution noint dist=&dist link=&link nocenter;
		%if &repeated eq 1 %then %do;
			random &time / subject = &id type = un residual;
		%end;
		%else %do;
			random Intercept / subject=&id residual;
		%end;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		ods output FitStatistics=fit CovParms=Cov;
		output out=pred pred(noblup ilink)=pred variance(noblup ilink)=var;
	run;

	%if &print_final eq 1 %then %do;
		ods exclude all;
	%end;

	*Rerunning step b and c and storing models;
	proc glimmix data=Step_b method=&method IC=PQ 
		%if &repeated eq 1 %then %do;
			EMPIRICAL
		%end;
		;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model &outcome(DESC) = &genes / noint solution dist=&dist link=&link offset=R0 nocenter;
		%if &repeated eq 1 %then %do;
			random &time / subject = &id type = un residual;
		%end;
		%else %do;
			random Intercept / subject=&id residual;
		%end;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		store final_model_G;
		parms / noiter PDATA=Cov;
		ods output ParameterEstimates=SolutionF_G;
	run;

	proc glimmix data=Step_c method=&method IC=PQ 
		%if &repeated eq 1 %then %do;
			EMPIRICAL
		%end;
		;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model &outcome(DESC) = &env / noint solution dist=&dist link=&link offset=R0 nocenter;
		%if &repeated eq 1 %then %do;
			random &time / subject = &id type = un residual;
		%end;
		%else %do;
			random Intercept / subject=&id residual;
		%end;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		store final_model_E;
		parms / noiter PDATA=Cov;
		ods output ParameterEstimates=SolutionF_E;
	run;

ods exclude none;

	* AIC;
	proc iml;
		use fit;
		read all;
		AIC = value[2];
		AIC = AIC + 2*(&env_N-1);
		AIC = AIC + 2*(&genes_N-1);
		create AIC var {AIC};
			append;
		%if &print_final eq 1 %then %do;
			print AIC;
		%end;
	quit;

* print final results if print_final=1;
%if &print_final eq 1 %then %do;

	* G;
	proc plm source=final_model_G;
			estimate 
			%do i=1 %to &genes_N;
				"%scan(&genes,&i)" %do j=1 %to &genes_N;
					%if &i ne &j %then %do;
						%scan(&genes,&j) 0
					%end;
					%else %do;
						%scan(&genes,&j) &total_inv_G
					%end;
				%end;
				%if &i ne &genes_N %then ,;
			%end;
			;
	run;

	proc iml;
		use conv_status1;
		read all;
		print diff_G diff_threshold iter total;
	quit;

	* E;
	proc plm source=final_model_E;
			estimate 
			%do i=1 %to &env_N;
				"%scan(&env,&i)" %do j=1 %to &env_N;
					%if &i ne &j %then %do;
						%scan(&env,&j) 0
					%end;
					%else %do;
						%scan(&env,&j) &total_inv_E
					%end;
				%end;
				%if &i ne &env_N %then ,;
			%end;
			;
	run;

	* Convergence;
	proc iml;
		use conv_status2;
		read all;
		print diff_E diff_threshold iter total;
	quit;

	data pred_data;
		merge Step_a(keep = &id &outcome) pred;
	run;
	data pred_data;
		set pred_data;
		Resid = pred-&outcome;
	run;
	* R2 creation;
	proc iml;
		use pred_data;
			read all;
			SSTotal = ssq(&outcome-mean(&outcome));
			R2 = 1 - ssq(Resid)/SSTotal;
		close pred_data;
	create R2 var {R2};
	append;
	print R2;
	 quit;

	* AUC and ROC curve;
	options minoperator mlogic;
	%if &dist in bin,beta,binary,b,binomial,multinomial,mult,multi, %then %do;
		proc logistic data=pred_data DESC plots(only)=roc;
			model &outcome =  / outroc=rocstats;
			roc pred=pred;
        	roccontrast;
		run;
		data check;
			set rocstats;
			_SPECIF_ = (1 - _1MSPEC_);
			J = _SENSIT_ + _SPECIF_ - 1;
			D= Sqrt((1-_SENSIT_)**2 + (1-_SPECIF_)**2);
		run;
		proc sql noprint;
			create table cutoff as
			select _PROB_ , J, _SENSIT_, _SPECIF_, _FALPOS_, _FALNEG_, _POS_, _NEG_
			from check
			having J = max(J);
		run;
		proc sql noprint;
			create table cutoff1 as
			select _PROB_ , D, _SENSIT_, _SPECIF_, _FALPOS_, _FALNEG_, _POS_, _NEG_
			from check
			having D = min(D);
		run;
		proc print data = cutoff;
			title1 'Best choice of cutoff based on Youden Index';
			var _PROB_ J _SENSIT_ _SPECIF_ _FALPOS_ _FALNEG_ _POS_ _NEG_;
		run;
		proc print data = cutoff1;
			title1 'Best choice of cutoff based on Euclidean Distance';
			var _PROB_ D _SENSIT_ _SPECIF_ _FALPOS_ _FALNEG_ _POS_ _NEG_;
		run;
		title1;
	%end;

	%if &clear_ods eq 1 %then %do;
		dm odsresults "clear";
	%end;

%end;
%mend;




*--------------------------------------------------------------------------------*
* Macro #2) LEGIT_mixed: Constructs a linear mixed model (using PROC MIXED) with *
*		a weighted latent environmental score and weighted latent genetic 		 *
* 		score. Output model, convergence statistics and AIC and approximated 	 *
* 		cross-validated R^2.										 			 *
*--------------------------------------------------------------------------------*
* Parameters necessary to run *
*-----------------------------*---------------------------------------------------------------------*
* Data : dataset to be used                                                                         *
* outcome : outcome variable                                                                        *
* genes : genes variables inside genetic score G 												    *
          (can be any sort of variable, doesn't even have to be genetic)						    *
* env : environment variables inside environmental score E											*
          (can be any sort of variable, doesn't even have to be environmental)					    *
* model_noG_noE : model part without E and without G (See model specification above) 			    *
* model_E_noG : model part with E and without G (See model specification above)                     *
* model_G_noE : model part without E and with G (See model specification above)                     *
* model_G_E : model part with E and with G (See model specification above)                          *
* id : ID of individual (if no ID, just create it using id=_n_ in a data step)                      *
*-------------------------------------------------*												    *
* Optional parameters or parameters with defaults *												    *
*-------------------------------------------------*-------------------------------------------------*
* covs (optional) : extra covariates (Equivalent to adding them in model_noG_noE)				    *
* remove_miss : if 1 then remove missing data before running (Default=1)						    *
* time (optional) : time variable for repeated measure outcome   								    *
* start_genes (optional) : starting points for genetic score (must be same length as "genes")       *
* start_env (optional) : starting points for environmental score (must be same length as "env")     *
* eps : threshold for convergence (.01 for quick batch simulations, .0001 for accurate results)     *
* maxiter : Maximum number of iterations														    *
* print : If 1 then print all models in all iterations, use for debugging only (Default=0)		    *
* print_final : If 1 then print final model (Default=1)											    *
* ods_new : If 1 then close current ods output and start a new one (Default=0)                      *
* repeated : If 1 then the outcome is a repeated measure (Default=0)                                *
* repeated_type : covariance type for repeated measure (Default= un)                                *
* random_vars (optional) : variables of random effects                                              *
* random_sub (optional) : subject for randon effects                                                *
* random_type : covariance type for random effects (Default= vc)                                    *
* where (optional) : where                                                                          *
* ods_clear : If 1 then clear ods cluttering garbage (Default=0)                                    *
*---------------------------------------------------------------------------------------------------*;
*----------*
* Examples *
*----------*
* two-way model : Intercept + g + e + ge + gender + ses
* Four genes g1, g2, g3, g4. Three environments e1, e2, e3. Two gene by gene interaction g1_g2, g1_g4. 
* Continuous outcome y at four time-points time="3M", "6M", "18M" or "36M".
*
* LEGIT_mixed(data, y, g1 g2 g3 g4 g1_g2 g1_g4, e1 e2 e3, Intercept, Intercept, Intercept, Intercept, gender ses, id=id, time=time, repeated=1)			
*
* three-way model : Intercept + g + e + z + ge + gz + ez + gez + gender + ses
* Four genes g1, g2, g3, g4. Three environments e1, e2, e3. Two gene by gene interaction g1_g2, g1_g4. 
* Continuous outcome y at one time-point.
*
* LEGIT_mixed(data, y, g1 g2 g3 g4 g1_g2 g1_g4, e1 e2 e3, Intercept z, Intercept z, Intercept z, Intercept z, gender ses, id=id)				  		


%macro LEGIT_mixed(data, outcome, genes, env, model_noG_noE=, model_E_noG=, model_G_noE=, model_G_E=, covs=, id=PSCID, remove_miss=1, time=, start_genes=, start_env=, eps = .0001, maxiter = 50, print = 0, print_final=1, ods_new=1, repeated=1, repeated_type=un, random_vars=, random_sub=, random_type=vc, where=, clear_ods=0);

* Won't give an error if trying to drop G or E when it doesn't exist;
options dkricond=warn;

* Counting how many variables;
%let genes_N = %sysfunc(countw(&genes));
%let env_N = %sysfunc(countw(&env));
%let model_noG_noE_N = %sysfunc(countw(&model_noG_noE,' '));
%let model_E_noG_N = %sysfunc(countw(&model_E_noG,' '));
%let model_G_noE_N = %sysfunc(countw(&model_G_noE,' '));
%let model_G_E_N = %sysfunc(countw(&model_G_E,' '));
%let covs_N = %sysfunc(countw(&covs,' '));

%if &ods_new eq 1 %then %do;
	ods html close;
	ods html;
%end;

* Sorting, setting up intercept and outlier removal; 
Proc sort data = &data;
	by &id 
%if &repeated eq 1 %then %do;
	&time 
%end;
;
run;
data Step_a;
	set &data;
	Intercept = 1;
run;

%if &remove_miss eq 1 %then %do;
	data Step_a (keep = &id 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		&outcome &genes &env &model_noG_noE &model_E_noG &model_G_noE &model_G_E &covs &random_vars &random_sub);
		set Step_a;
			where &where;
	run;
	data Step_a;
		set Step_a;
		if cmiss(of _all_) then delete;
	run;
%end;

* Setting up datasets with initial weighted scores;
%if %bquote(&start_genes) eq %then %do;
	data G_weights;
		G_weights_old = 1/&genes_N;
		%do i=1 %to &genes_N;
			output;
		%end;
	run;
%end;
%else %do;
	proc iml;
		G_weights_old = {&start_genes};
		create G_weights var {G_weights_old};
		append;
	quit;
%end;
%if %bquote(&start_env) eq %then %do;
	data E_weights;
		E_weights_old = 1/&env_N;
		%do i=1 %to &env_N;
			output;
		%end;
	run;
%end;
%else %do;
	proc iml;
		E_weights_old = {&start_env};
		create E_weights var {E_weights_old};
		append;
	quit;
%end;

* Setting up for first step;
proc iml;
	use Step_a;
		read all var {&genes} into Genes[colname=varnames];
		read all var {&env} into Env[colname=varnames];
		read all var {&id};
	close Step_a;
	use G_weights;
		read all;
		G = Genes*G_weights_old/(sum(abs(G_weights_old)));
	close G_weights;
	use E_weights;
		read all;
		E = Env*E_weights_old/(sum(abs(E_weights_old)));
	close E_weights;
	create Step_a_GE var {&id G E};
	append;
	%if &print eq 1 %then %do;
		print G E;
	%end;
quit;
data Step_a;
	* Merging and dropping G and E just in case that they already exist in original dataset;
	merge Step_a_GE Step_a(drop = G E);
run;

%do i=1 %to &maxiter;

	%if &print ne 1 %then %do;
		ods exclude all;
	%end;

	* Step a : Estimating betas and covariance matrix;
	proc mixed data=Step_a method=ML;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model &outcome = &model_noG_noE 
			%do j=1 %to &model_E_noG_N;
				E*%scan(&model_E_noG, &j,' ') 
			%end;
			%do j=1 %to &model_G_noE_N;
				G*%scan(&model_G_noE, &j,' ') 
			%end;
			%do j=1 %to &model_G_E_N;
				E*G*%scan(&model_G_E, &j,' ') 
			%end;
		&covs / solution noint;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		%if &repeated eq 1 %then %do;
			repeated &time / subject = &id type = &repeated_type;
		%end;
		ods output SolutionF=SolutionF CovParms=Cov;
	run;

	*Variable modification and setting up for b step (estimating G);
	proc iml;
		use Step_a;
			read all var {&genes} into Genes[colname=varnames];
			read all var {E};
			read all var {&model_noG_noE} into X_model_noG_noE[colname=varnames];
			read all var {&model_E_noG} into X_model_E_noG[colname=varnames];
			read all var {&model_G_noE} into X_model_G_noE[colname=varnames];
			read all var {&model_G_E} into X_model_G_E[colname=varnames];
			read all var {&covs} into X_covs[colname=varnames];
			read all var {&outcome};
		close Step_a;
		use SolutionF;
			read all var{Estimate};
		close SolutionF;
		R0 = X_model_noG_noE*Estimate[1:(&model_noG_noE_N)] + (X_model_E_noG#E)*Estimate[(&model_noG_noE_N+1):(&model_noG_noE_N+&model_E_noG_N)] + X_covs*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+&covs_N)];
		R1 = X_model_G_noE*Estimate[(&model_noG_noE_N+&model_E_noG_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N)] + (X_model_G_E#E)*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N)];
		R1_Genes = Genes#R1;
		Y_new = &outcome - R0;
		R = Y_new||R1_Genes;
		create Step_b from R[colname={Y_new &genes}];
		append from R;		
		%if &print eq 1 %then %do;
			print R0 R1 R1_genes Y_new;
		%end;
	quit;
	data Step_b;
		merge Step_a(keep = &id &random_vars &random_sub  
			%if &repeated eq 1 %then %do;
				&time
			%end;
		) Step_b;
	run;

	*Step b: Estimating G weights;
	proc mixed data=Step_b method=ML;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model Y_new = &genes / noint solution;
		%if &repeated eq 1 %then %do;
			repeated &time / subject = &id type = un;
		%end;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		parms / noiter PDATA=Cov;
		ods output SolutionF=SolutionF_G;
	run;

	*Updating G;
	proc iml;
		use Step_a;
			read all var {G};
			read all var {&id};
			read all var {&genes} into Genes[colname=varnames];
		close Step_a;
		G_old = G;
		use SolutionF_G;
			read all;
			use G_weights;
				read all var {G_weights_old};
				G = Genes*Estimate/(sum(abs(Estimate)));
				G_weights_new = Estimate/(sum(abs(Estimate)));
				* term to multiply estimates;
				total_inv_G = 1/(sum(abs(Estimate)));
				call symputx("total_inv_G", total_inv_G);
				* Convergence stuff;
				diff_G = sqrt(ssq(G_weights_new - G_weights_old));
				iter = num(symget("i"));
				diff_threshold = num(symget("eps"));
				if (diff_G < diff_threshold) then call symputx("conv_half1", 1);
				else call symputx("conv_half1", 0);
				* Updating estimates;
				G_weights_old = G_weights_new;
			close G_weights;
			create G_weights var {G_weights_old};
				append;
			create Step_b_G var {&id G};
				append;
			* To be reported at end of macro as output for users;
			create conv_status1 var {diff_G diff_threshold iter total};
				append;
		close SolutionF_G;
		%if &print eq 1 %then %do;
			print G_weights_new diff_G diff_threshold iter G;
		%end;
	quit;
	data Step_a;
		merge Step_b_G Step_a(drop = G);
	run;

	*Variable modification and setting up for c step (estimating E);
	proc iml;
		use Step_a;
			read all var {&env} into Env[colname=varnames];
			read all var {G};
			read all var {&model_noG_noE} into X_model_noG_noE[colname=varnames];
			read all var {&model_E_noG} into X_model_E_noG[colname=varnames];
			read all var {&model_G_noE} into X_model_G_noE[colname=varnames];
			read all var {&model_G_E} into X_model_G_E[colname=varnames];
			read all var {&covs} into X_covs[colname=varnames];
			read all var {&outcome};
		close Step_a;
		use SolutionF;
			read all var{Estimate};
		close SolutionF;
		* Setting special parametrization for next E step;
		* Estimates are currently ordered as &model_noG_noE -> &model_E_noG -> &model_G_noE -> &model_G_E -> covs;
		* We need them in this order : &model_noG_noE -> &model_G_noE -> &model_E_noG -> &model_G_E -> covs;
		R0 = X_model_noG_noE*Estimate[1:(&model_noG_noE_N)] + (X_model_G_noE#G)*Estimate[(&model_noG_noE_N+&model_E_noG_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N)] + X_covs*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+&covs_N)];
		R1 = X_model_E_noG*Estimate[(&model_noG_noE_N+1):(&model_noG_noE_N+&model_E_noG_N)] + (X_model_G_E#G)*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N)];
		R1_Env = Env#R1;
		Y_new = &outcome - R0;
		R = Y_new||R1_Env;
		create Step_c from R[colname={Y_new &env}];
		append from R;
		%if &print eq 1 %then %do;
			print R0 R1 R1_Env Y_new;
		%end;
	quit;
	data Step_c;
		merge Step_a(keep = &id &random_vars &random_sub 
		%if &repeated eq 1 %then %do;
			&time
		%end;
		) Step_c;
	run;
	*Step c: Estimating E weights;
	proc mixed data=Step_c method=ML;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model Y_new = &env / noint solution;
		%if &repeated eq 1 %then %do;
			repeated &time / subject = &id type = un;
		%end;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		parms / noiter PDATA=Cov;
		ods output SolutionF=SolutionF_E;
	run;
	* Normalizing weights and combing everything into dataset for step a;
	proc iml;
		use Step_a;
			read all var {&env} into Env[colname=varnames];
			read all var {&id};
			read all var {E};
		close Step_a;
		E_old = E;
		use SolutionF_E;
			read all;
			use E_weights;
				read all var {E_weights_old};
				E = Env*Estimate/(sum(abs(Estimate)));
				E_weights_new = Estimate/(sum(abs(Estimate)));
				* term to multiply estimates;
				total_inv_E = 1/(sum(abs(Estimate)));
				call symputx("total_inv_E", total_inv_E);
				* Looking for convergence;
				diff_E = sqrt(ssq(E_weights_new - E_weights_old));
				iter = num(symget("i"));
				diff_threshold = num(symget("eps"));
				if (diff_E < diff_threshold) then call symputx("conv_half2", 1);
				else call symputx("conv_half2", 0);
				* Updating estimates;
				E_weights_old = E_weights_new;
			close E_weights;
			* Creating datasets;
			create E_weights var {E_weights_old};
				append;
			create Step_c_E var {&id E};
				append;
			close Step_c_E;
			create conv_status2 var {diff_E diff_threshold iter total};
				append;
		close SolutionF_E;
		%if &print eq 1 %then %do;
			print E_weights_new diff_E diff_threshold iter;
		%end;
	quit;
	data Step_a;
		merge Step_c_E Step_a(drop = E);
	run;
	* If both E and G converged then stop;
	%if (&conv_half1 eq 1 and &conv_half2 eq 1) %then %let i=&maxiter;
%end;

	%if &print_final eq 1 %then %do;
		ods exclude none;
	%end;

	*Rerunning step a for AIC;
	proc mixed data=Step_a method=ML COVTEST RATIO ;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model &outcome = &model_noG_noE 
			%do j=1 %to &model_E_noG_N;
				E*%scan(&model_E_noG, &j,' ') 
			%end;
			%do j=1 %to &model_G_noE_N;
				G*%scan(&model_G_noE, &j,' ') 
			%end;
			%do j=1 %to &model_G_E_N;
				E*G*%scan(&model_G_E, &j,' ') 
			%end;
		&covs / solution noint;
		%if &repeated eq 1 %then %do;
			repeated &time / subject = &id type = un;
		%end;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		ods output FitStatistics=fit CovParms=Cov;
	run;

	%if &print_final eq 1 %then %do;
		ods exclude all;
	%end;

	%if &print_final eq 1 %then %do;
		* Rerunning again step a again to get leave-one-out predictions for cross-validated R^2;
		proc mixed data=Step_a method=ML COVTEST RATIO ;
			id &id;
			class &id &random_vars 
			%if &repeated eq 1 %then %do;
				&time 
			%end;
			;
			model &outcome = &model_noG_noE 
				%do j=1 %to &model_E_noG_N;
					E*%scan(&model_E_noG, &j,' ') 
				%end;
				%do j=1 %to &model_G_noE_N;
					G*%scan(&model_G_noE, &j,' ') 
				%end;
				%do j=1 %to &model_G_E_N;
					E*G*%scan(&model_G_E, &j,' ') 
				%end;
			&covs / outp=pred RESIDUAL influence(ITER=0 EFFECT=&id) solution noint;
			%if &repeated eq 1 %then %do;
				repeated &time / subject = &id type = un;
			%end;
			%if &random_vars ne %then %do;
				random &random_vars / 
				%if &random_sub ne %then %do;
					subject = &random_sub 
				%end;
				type = &random_type;
			%end;
			ods output INFLUENCE=inf;
		run;
	%end;

	*Rerunning step b and c and storing models;
	proc mixed data=Step_b method=ML COVTEST RATIO 
		%if &repeated eq 1 %then %do;
			EMPIRICAL
		%end;
		;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model Y_new = &genes / noint solution;
		%if &repeated eq 1 %then %do;
			repeated &time / subject = &id type = un;
		%end;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		store final_model_G;
		parms / noiter PDATA=Cov;
		ods output SolutionF=SolutionF_G;
	run;

	proc mixed data=Step_c method=ML COVTEST RATIO 
		%if &repeated eq 1 %then %do;
			EMPIRICAL
		%end;
		;
		id &id;
		class &id &random_vars 
		%if &repeated eq 1 %then %do;
			&time 
		%end;
		;
		model Y_new = &env / noint solution;
		%if &repeated eq 1 %then %do;
			repeated &time / subject = &id type = un;
		%end;
		%if &random_vars ne %then %do;
			random &random_vars / 
			%if &random_sub ne %then %do;
				subject = &random_sub 
			%end;
			type = &random_type;
		%end;
		store final_model_E;
		parms / noiter PDATA=Cov;
		ods output SolutionF=SolutionF_E;
	run;

	ods exclude none;

	* AIC;
	proc iml;
		use fit;
		read all;
		AIC = value[2];
		AIC = AIC + 2*(&env_N-1);
		AIC = AIC + 2*(&genes_N-1);
		create AIC var {AIC};
			append;
		%if &print_final eq 1 %then %do;
			print AIC;
		%end;
	quit;

* print final results if print_final=1;
%if &print_final eq 1 %then %do;

	* G;
	proc plm source=final_model_G;
			estimate 
			%do i=1 %to &genes_N;
				"%scan(&genes,&i)" %do j=1 %to &genes_N;
					%if &i ne &j %then %do;
						%scan(&genes,&j) 0
					%end;
					%else %do;
						%scan(&genes,&j) &total_inv_G
					%end;
				%end;
				%if &i ne &genes_N %then ,;
			%end;
			;
	run;

	proc iml;
		use conv_status1;
		read all;
		print diff_G diff_threshold iter total;
	quit;

	* E;
	proc plm source=final_model_E;
			estimate 
			%do i=1 %to &env_N;
				"%scan(&env,&i)" %do j=1 %to &env_N;
					%if &i ne &j %then %do;
						%scan(&env,&j) 0
					%end;
					%else %do;
						%scan(&env,&j) &total_inv_E
					%end;
				%end;
				%if &i ne &env_N %then ,;
			%end;
			;
	run;

	* Convergence;
	proc iml;
		use conv_status2;
		read all;
		print diff_E diff_threshold iter total;
	quit;

	data pred_data;
		merge Step_a(keep = &id &outcome) pred(keep=&id pred);
	run;
	data pred_data;
		set pred_data;
		Resid = pred-&outcome;
	run;
	* R2 creation (not relevant for binary outcomes but will still be shown);
	proc iml;
		use inf;
			read all;
			PRESS_SSRes = sum(PRESS);
		close inf;
		use pred_data;
			read all;
			SSTotal = ssq(&outcome-mean(&outcome));
			R2 = 1 - ssq(Resid)/SSTotal;
			R2_leave_one_out = 1 - PRESS_SSRes/SSTotal;
		close pred_data;
	create R2 var {R2 R2_leave_one_out};
	append;
	title1 'R2 and approximated leave-one-out crossvalidated R2 (This assume that the genetic and environmentals scores are known, this is not the true cross-validated R2 but just a quick approximation, run macro LEGIT_cv to get the correct one (very slow)).';
	print R2 R2_leave_one_out;
	title1;
	quit;

	%if &clear_ods eq 1 %then %do;
		dm odsresults "clear";
	%end;

%end;
%mend;




*------------------------------------------------------------------------------------*
* Macro #3) LEGIT_logistic: Constructs a logistic model (using PROC LOGISTIC) with 	 *
*		a weighted latent environmental score and weighted latent genetic 		     *
* 		score. Output model, convergence statistics and AIC and approximated 	 	 *
* 		cross-validated AUC.										 			 	 *
*------------------------------------------------------------------------------------*
* Parameters necessary to run *
*-----------------------------*------------------------------------------------------------------*
* Data : dataset to be used                                                                      *
* outcome : outcome variable                                                                     *
* genes : genes variables inside genetic score G 												 *
          (can be any sort of variable, doesn't even have to be genetic)						 *
* env : environment variables inside environmental score E										 *
          (can be any sort of variable, doesn't even have to be environmental)					 *
* model_noG_noE : model part without E and without G (See model specification above) 			 *
* model_E_noG : model part with E and without G (See model specification above)                  *
* model_G_noE : model part without E and with G (See model specification above)                  *
* model_G_E : model part with E and with G (See model specification above)                       *
* id : ID of individual (if no ID, just create it using id=_n_ in a data step)                   *
*-------------------------------------------------*												 *
* Optional parameters or parameters with defaults *												 *
*-------------------------------------------------*----------------------------------------------*
* covs (optional) : extra covariates (Equivalent to adding them in model_noG_noE)				 *
* remove_miss : if 1 then remove missing data before running (Default=1)						 *
* start_genes (optional) : starting points for genetic score (must be same length as "genes")    *
* start_env (optional) : starting points for environmental score (must be same length as "env")  *
* eps : threshold for convergence (.01 for quick batch simulations, .0001 for accurate results)  *
* maxiter : Maximum number of iterations (Default=50)											 *
* print : If 1 then print all models in all iterations, use for debugging only (Default=0)		 *
* print_final : If 1 then print final model (Default=1)											 *
* ods_new : If 1 then close current ods output and start a new one (Default=0)                   *
* where (optional) : where                                                                       *
* ods_clear : If 1 then clear ods cluttering garbage (Default=0)                                 *
* link : GLM link function (Default=logit)														 *
*------------------------------------------------------------------------------------------------*;
*----------*
* Examples *
*----------*
* two-way model : Intercept + g + e + ge + gender + ses
* Four genes g1, g2, g3, g4. Three environments e1, e2, e3. Two gene by gene interaction g1_g2, g1_g4. 
* Binary outcome y at one time-point with logit link.
*
* LEGIT_logistic(data, y, g1 g2 g3 g4 g1_g2 g1_g4, e1 e2 e3, Intercept, Intercept, Intercept, Intercept, gender ses, id=id, link=logit)					
*
* three-way model : Intercept + g + e + z + ge + gz + ez + gez + gender + ses
* Four genes g1, g2, g3, g4. Three environments e1, e2, e3. Two gene by gene interaction g1_g2, g1_g4. 
* Binary outcome y at one time-point with loglog link.
*
* LEGIT_logistic(data, y, g1 g2 g3 g4 g1_g2 g1_g4, e1 e2 e3, Intercept z, Intercept z, Intercept z, Intercept z, gender ses, id=id, link=loglog)		


%macro LEGIT_logistic(data, outcome, genes, env, model_noG_noE=, model_E_noG=, model_G_noE=, model_G_E=, covs=, id=PSCID, remove_miss=1, start_genes=, start_env=, eps = .0001, maxiter = 50, print = 0, print_final=1, ods_new=1, where=, clear_ods=0, link=logit);

* Won't give an error if trying to drop G or E when it doesn't exist;
options dkricond=warn;

* Counting how many variables;
%let genes_N = %sysfunc(countw(&genes));
%let env_N = %sysfunc(countw(&env));
%let model_noG_noE_N = %sysfunc(countw(&model_noG_noE,' '));
%let model_E_noG_N = %sysfunc(countw(&model_E_noG,' '));
%let model_G_noE_N = %sysfunc(countw(&model_G_noE,' '));
%let model_G_E_N = %sysfunc(countw(&model_G_E,' '));
%let covs_N = %sysfunc(countw(&covs,' '));

%if &ods_new eq 1 %then %do;
	ods html close;
	ods html;
%end;

* Sorting, setting up intercept and outlier removal; 
Proc sort data = &data;
	by &id;
run;
data Step_a;
	set &data;
	Intercept = 1;
run;

%if &remove_miss eq 1 %then %do;
	data Step_a (keep = &id &outcome &genes &env &model_noG_noE &model_E_noG &model_G_noE &model_G_E &covs);
		set Step_a;
			where &where;
	run;
	data Step_a;
		set Step_a;
		if cmiss(of _all_) then delete;
	run;
%end;

* Setting up datasets with initial weighted scores;
%if %bquote(&start_genes) eq %then %do;
	data G_weights;
		G_weights_old = 1/&genes_N;
		%do i=1 %to &genes_N;
			output;
		%end;
	run;
%end;
%else %do;
	proc iml;
		G_weights_old = {&start_genes};
		create G_weights var {G_weights_old};
		append;
	quit;
%end;
%if %bquote(&start_env) eq %then %do;
	data E_weights;
		E_weights_old = 1/&env_N;
		%do i=1 %to &env_N;
			output;
		%end;
	run;
%end;
%else %do;
	proc iml;
		E_weights_old = {&start_env};
		create E_weights var {E_weights_old};
		append;
	quit;
%end;

* Setting up for first step;
proc iml;
	use Step_a;
		read all var {&genes} into Genes[colname=varnames];
		read all var {&env} into Env[colname=varnames];
		read all var {&id};
	close Step_a;
	use G_weights;
		read all;
		G = Genes*G_weights_old/(sum(abs(G_weights_old)));
	close G_weights;
	use E_weights;
		read all;
		E = Env*E_weights_old/(sum(abs(E_weights_old)));
	close E_weights;
	create Step_a_GE var {&id G E};
	append;
	%if &print eq 1 %then %do;
		print G E;
	%end;
quit;
data Step_a;
	* Merging and dropping G and E just in case that they already exist in original dataset;
	merge Step_a_GE Step_a(drop = G E);
run;

%do i=1 %to &maxiter;

	%if &print ne 1 %then %do;
		ods exclude all;
	%end;

	* Step a : Estimating betas and covariance matrix;
	proc logistic data=Step_a DESC;
		id &id;
		class &id;
		model &outcome = &model_noG_noE 
			%do j=1 %to &model_E_noG_N;
				E*%scan(&model_E_noG, &j,' ') 
			%end;
			%do j=1 %to &model_G_noE_N;
				G*%scan(&model_G_noE, &j,' ') 
			%end;
			%do j=1 %to &model_G_E_N;
				E*G*%scan(&model_G_E, &j,' ') 
			%end;
		&covs / noint link=&link;
		ods output ParameterEstimates=SolutionF;
	run;

	*Variable modification and setting up for b step (estimating G);
	proc iml;
		use Step_a;
			read all var {&genes} into Genes[colname=varnames];
			read all var {E};
			read all var {&model_noG_noE} into X_model_noG_noE[colname=varnames];
			read all var {&model_E_noG} into X_model_E_noG[colname=varnames];
			read all var {&model_G_noE} into X_model_G_noE[colname=varnames];
			read all var {&model_G_E} into X_model_G_E[colname=varnames];
			read all var {&covs} into X_covs[colname=varnames];
			read all var {&outcome};
		close Step_a;
		use SolutionF;
			read all var{Estimate};
		close SolutionF;
		R0 = X_model_noG_noE*Estimate[1:(&model_noG_noE_N)] + (X_model_E_noG#E)*Estimate[(&model_noG_noE_N+1):(&model_noG_noE_N+&model_E_noG_N)] + X_covs*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+&covs_N)];
		R1 = X_model_G_noE*Estimate[(&model_noG_noE_N+&model_E_noG_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N)] + (X_model_G_E#E)*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N)];
		R1_Genes = Genes#R1;
		R = &outcome||R0||R1_Genes;
		create Step_b from R[colname={&outcome R0 &genes}];
		append from R;		
		%if &print eq 1 %then %do;
			print R0 R1 R1_genes Y_new;
		%end;
	quit;
	data Step_b;
		merge Step_a(keep = &id) Step_b;
	run;

	*Step b: Estimating G weights;
	proc logistic data=Step_b DESC;
		id &id;
		class &id;
		model &outcome = &genes / noint link=&link offset=R0;
		ods output ParameterEstimates=SolutionF_G;
	run;
	*Removing offset from estimates;
	data SolutionF_G;
		set SolutionF_G end=last;
	    if last then delete;
	run;
	*Updating G;
	proc iml;
		use Step_a;
			read all var {G};
			read all var {&id};
			read all var {&genes} into Genes[colname=varnames];
		close Step_a;
		G_old = G;
		use SolutionF_G;
			read all;
			use G_weights;
				read all var {G_weights_old};
				G = Genes*Estimate/(sum(abs(Estimate)));
				G_weights_new = Estimate/(sum(abs(Estimate)));
				* term to multiply estimates;
				total_inv_G = 1/(sum(abs(Estimate)));
				call symputx("total_inv_G", total_inv_G);
				* Convergence stuff;
				diff_G = sqrt(ssq(G_weights_new - G_weights_old));
				iter = num(symget("i"));
				diff_threshold = num(symget("eps"));
				if (diff_G < diff_threshold) then call symputx("conv_half1", 1);
				else call symputx("conv_half1", 0);
				* Updating estimates;
				G_weights_old = G_weights_new;
			close G_weights;
			create G_weights var {G_weights_old};
				append;
			create Step_b_G var {&id G};
				append;
			* To be reported at end of macro as output for users;
			create conv_status1 var {diff_G diff_threshold iter total};
				append;
		close SolutionF_G;
		%if &print eq 1 %then %do;
			print G_weights_new diff_G diff_threshold iter;
		%end;
	quit;
	data Step_a;
		merge Step_b_G Step_a(drop = G);
	run;

	*Variable modification and setting up for c step (estimating E);
	proc iml;
		use Step_a;
			read all var {&env} into Env[colname=varnames];
			read all var {G};
			read all var {&model_noG_noE} into X_model_noG_noE[colname=varnames];
			read all var {&model_E_noG} into X_model_E_noG[colname=varnames];
			read all var {&model_G_noE} into X_model_G_noE[colname=varnames];
			read all var {&model_G_E} into X_model_G_E[colname=varnames];
			read all var {&covs} into X_covs[colname=varnames];
			read all var {&outcome};
		close Step_a;
		use SolutionF;
			read all var{Estimate};
		close SolutionF;
		* Setting special parametrization for next E step;
		* Estimates are currently ordered as &model_noG_noE -> &model_E_noG -> &model_G_noE -> &model_G_E -> covs;
		* We them in this order : &model_noG_noE -> &model_G_noE -> &model_E_noG -> &model_G_E -> covs;
		R0 = X_model_noG_noE*Estimate[1:(&model_noG_noE_N)] + (X_model_G_noE#G)*Estimate[(&model_noG_noE_N+&model_E_noG_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N)] + X_covs*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N+&covs_N)];
		R1 = X_model_E_noG*Estimate[(&model_noG_noE_N+1):(&model_noG_noE_N+&model_E_noG_N)] + (X_model_G_E#G)*Estimate[(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+1):(&model_noG_noE_N+&model_E_noG_N+&model_G_noE_N+&model_G_E_N)];
		R1_Env = Env#R1;
		R = &outcome||R0||R1_Env;
		create Step_c from R[colname={&outcome R0 &env}];
		append from R;
		%if &print eq 1 %then %do;
			print R0 R1 R1_Env Y_new;
		%end;
	quit;
	data Step_c;
		merge Step_a(keep = &id) Step_c;
	run;
	*Step c: Estimating E weights;
	proc logistic data=Step_c DESC;
		id &id;
		class &id;
		model &outcome = &env / noint link=&link offset=R0;
		ods output ParameterEstimates=SolutionF_E;
	run;
	*Removing offset from estimates;
	data SolutionF_E;
		set SolutionF_E end=last;
	    if last then delete;
	run;
	* Normalizing weights and combing everything into dataset for step a;
	proc iml;
		use Step_a;
			read all var {&env} into Env[colname=varnames];
			read all var {&id};
			read all var {E};
		close Step_a;
		E_old = E;
		use SolutionF_E;
			read all;
			use E_weights;
				read all var {E_weights_old};
				E = Env*Estimate/(sum(abs(Estimate)));
				E_weights_new = Estimate/(sum(abs(Estimate)));
				* term to multiply estimates;
				total_inv_E = 1/(sum(abs(Estimate)));
				call symputx("total_inv_E", total_inv_E);
				* Looking for convergence;
				diff_E = sqrt(ssq(E_weights_new - E_weights_old));
				iter = num(symget("i"));
				diff_threshold = num(symget("eps"));
				if (diff_E < diff_threshold) then call symputx("conv_half2", 1);
				else call symputx("conv_half2", 0);
				* Updating estimates;
				E_weights_old = E_weights_new;
			close E_weights;
			* Creating datasets;
			create E_weights var {E_weights_old};
				append;
			create Step_c_E var {&id E};
				append;
			close Step_c_E;
			create conv_status2 var {diff_E diff_threshold iter total};
				append;
		close SolutionF_E;
		%if &print eq 1 %then %do;
			print E_weights_new diff_E diff_threshold iter;
		%end;
	quit;
	data Step_a;
		merge Step_c_E Step_a(drop = E);
	run;
	* If both E and G converged then stop;
	%if (&conv_half1 eq 1 and &conv_half2 eq 1) %then %let i=&maxiter;
%end;

	%if &print_final eq 1 %then %do;
		ods exclude none;
	%end;
	
	*Rerunning step a for AIC;
	proc logistic data=Step_a DESC;
		id &id;
		class &id;
		model &outcome = &model_noG_noE 
			%do j=1 %to &model_E_noG_N;
				E*%scan(&model_E_noG, &j,' ') 
			%end;
			%do j=1 %to &model_G_noE_N;
				G*%scan(&model_G_noE, &j,' ') 
			%end;
			%do j=1 %to &model_G_E_N;
				E*G*%scan(&model_G_E, &j,' ') 
			%end;
		&covs / noint link=&link lackfit;
		ods output FitStatistics=fit;
		output out=pred PRED=pred 
		%if &print_final eq 1 %then %do;
			PREDPROBS=(CROSSVALIDATE) 
		%end;
		;
	run;

	%if &print_final eq 1 %then %do;
		ods exclude all;
	%end;

	*Rerunning step b and c and storing models;
	proc logistic data=Step_b DESC;
		id &id;
		class &id;
		model &outcome = &genes / noint link=&link offset=R0;
		store final_model_G;
		ods output ParameterEstimates=SolutionF_G;
	run;
	*Removing offset from estimates;
	data SolutionF_G;
		set SolutionF_G end=last;
	    if last then delete;
	run;

	proc logistic data=Step_c DESC;
		id &id;
		class &id;
		model &outcome = &env / noint link=&link offset=R0;
		store final_model_E;
		ods output ParameterEstimates=SolutionF_E;	
	run;
	*Removing offset from estimates;
	data SolutionF_E;
		set SolutionF_E end=last;
	    if last then delete;
	run;

ods exclude none;

	* AIC;
	proc iml;
		use fit;
		read all;
		AIC = WithCovariates[1];
		AIC = AIC + 2*(&env_N-1);
		AIC = AIC + 2*(&genes_N-1);
		create AIC var {AIC};
			append;
		%if &print_final eq 1 %then %do;
			print AIC;
		%end;
	quit;

* print final results if print_final=1;
%if &print_final eq 1 %then %do;

	* G;
	proc plm source=final_model_G;
			estimate 
			%do i=1 %to &genes_N;
				"%scan(&genes,&i)" %do j=1 %to &genes_N;
					%if &i ne &j %then %do;
						%scan(&genes,&j) 0
					%end;
					%else %do;
						%scan(&genes,&j) &total_inv_G
					%end;
				%end;
				%if &i ne &genes_N %then ,;
			%end;
			;
	run;

	proc iml;
		use conv_status1;
		read all;
		print diff_G diff_threshold iter total;
	quit;

	* E;
	proc plm source=final_model_E;
			estimate 
			%do i=1 %to &env_N;
				"%scan(&env,&i)" %do j=1 %to &env_N;
					%if &i ne &j %then %do;
						%scan(&env,&j) 0
					%end;
					%else %do;
						%scan(&env,&j) &total_inv_E
					%end;
				%end;
				%if &i ne &env_N %then ,;
			%end;
			;
	run;

	* Convergence;
	proc iml;
		use conv_status2;
		read all;
		print diff_E diff_threshold iter total;
	quit;

	* R2 creation;
	proc iml;
		use pred_data;
			read all;
			SSTotal = ssq(&outcome-mean(&outcome));
			R2 = 1 - ssq(Resid)/SSTotal;
		close pred_data;
	create R2 var {R2};
	append;
	print R2;
	 quit;

	* Combine predictions with step a dataset;
		data pred_data;
			merge Step_a(keep = &id &outcome) pred;
		run;
		data pred_data;
			set pred_data;
			Resid = pred-&outcome;
		run;

	* AUC and ROC curve;
		title1 'Note : This assume that the genetic and environmentals scores are known, this is not the true cross-validated ROC curve but just a quick approximation, run macro LEGIT_cv to get the correct one (very slow).';
		ods select ROCOverlay ROCAssociation ROCContrastTest;
		ods output Logistic.ROC1.ROCCurve=roccurve_cv;
		proc logistic data=pred_data DESC;
			id &id;
			class &id;
			model &outcome = &model_noG_noE 
			%do j=1 %to &model_E_noG_N;
				E*%scan(&model_E_noG, &j,' ') 
			%end;
			%do j=1 %to &model_G_noE_N;
				G*%scan(&model_G_noE, &j,' ') 
			%end;
			%do j=1 %to &model_G_E_N;
				E*G*%scan(&model_G_E, &j,' ') 
			%end;
			&covs / noint link=&link;
	        roc "Leave-one-out Crossvalidated" pred=xp_1;
	        roccontrast;
	    run;

	* Best cutoff non-crossvalidated;
		data check;
			set rocstats;
			_SPECIF_ = (1 - _1MSPEC_);
			J = _SENSIT_ + _SPECIF_ - 1;
			D= Sqrt((1-_SENSIT_)**2 + (1-_SPECIF_)**2);
		run;
		proc sql noprint;
			create table cutoff as
			select _PROB_ , J, _SENSIT_, _SPECIF_, _FALPOS_, _FALNEG_, _POS_, _NEG_
			from check
			having J = max(J);
		run;
		proc sql noprint;
			create table cutoff1 as
			select _PROB_ , D, _SENSIT_, _SPECIF_, _FALPOS_, _FALNEG_, _POS_, _NEG_
			from check
			having D = min(D);
		run;
		proc print data = cutoff;
			title1 'Best choice of cutoff based on Youden Index';
			var _PROB_ J _SENSIT_ _SPECIF_ _FALPOS_ _FALNEG_ _POS_ _NEG_;
		run;
		proc print data = cutoff1;
			title1 'Best choice of cutoff based on Euclidean Distance';
			var _PROB_ D _SENSIT_ _SPECIF_ _FALPOS_ _FALNEG_ _POS_ _NEG_;
		run;
		title1;

	* Best cutoff cross-validated;
	proc sql noprint;
		create table outcome_info as
		select sum(&outcome) as outcome_pos, sum(1-&outcome) as outcome_neg  
		from pred_data;
	quit;
	data outcome_info (drop = i);
		set outcome_info;
		do i=1 to (outcome_pos+outcome_neg+1);
			output;
		end;
	run;
	data check_cv;
		merge outcome_info roccurve_cv;
	run;
	data check_cv;
		set check_cv;
		_SPECIF_ = (1 - _1MSPEC_);
		J = _SENSIT_ + _SPECIF_ - 1;
		D= Sqrt((1-_SENSIT_)**2 + (1-_SPECIF_)**2);
		_FALPOS_ = (1-_SPECIF_)*outcome_neg;
		_FALNEG_ = (1-_SENSIT_)*outcome_pos;
		_POS_ = _SENSIT_*outcome_pos;
		_NEG_ = _SPECIF_*outcome_neg;
	run;
	proc sql noprint;
		create table cutoff_cv as
		select _PROB_ , J, _SENSIT_, _SPECIF_, _FALPOS_, _FALNEG_, _POS_,_NEG_
		from check_cv
		having J = max(J);
	run;
	proc sql noprint;
		create table cutoff1_cv as
		select _PROB_ , D, _SENSIT_, _SPECIF_, _FALPOS_, _FALNEG_, _POS_, _NEG_
		from check_cv
		having D = min(D);
	run;
	proc print data = cutoff_cv;
		title1 'Best choice of cutoff based on Youden Index (Crossvalidated)';
		var _PROB_ J _SENSIT_ _SPECIF_ _FALPOS_ _FALNEG_ _POS_ _NEG_;
	run;
	proc print data = cutoff1_cv;
		title1 'Best choice of cutoff based on Euclidean Distance (Crossvalidated)';
		var _PROB_ D _SENSIT_ _SPECIF_ _FALPOS_ _FALNEG_ _POS_ _NEG_;
	run;
	title1;

	%if &clear_ods eq 1 %then %do;
		dm odsresults "clear";
	%end;

%end;
 
%mend;




*--------------------------------------------------------------------------------*
* Macro #4) LEGIT_cv: Calculates the leave-one-out (removing individuals rather  *
*		than observations if repeated-measure model) cross-validated R^2. Also   *
*		shows the individuals with extreme residuals for outlier detection and   *
*		removal. 																 *
*--------------------------------------------------------------------------------*
* Parameters necessary to run *
*-----------------------------*------------------------------------------------------------------*
* Data : dataset to be used                                                                      *
* outcome : outcome variable                                                                     *
* genes : genes variables inside genetic score G 												 *
*         (can be any sort of variable, doesn't even have to be genetic)						 *
* env : environment variables inside environmental score E 										 *
*          (can be any sort of variable, doesn't even have to be environmental)					 *
* proc_used : glimmix = PROC GLIMMIX, logistic = PROC LOGISTIC, mixed = PROC MIXED               *
* model_noG_noE : model part without E and without G (See model specification above) 			 *
* model_E_noG : model part with E and without G (See model specification above)                  *
* model_G_noE : model part without E and with G (See model specification above)                  *
* model_G_E : model part with E and with G (See model specification above)                       *
* id : variable of the ID of each individual (ex : "Bob") 										 *
* true_id : variable of the ID of each observation (ex : "Bob_2010","Bob_2012","Bob_2016") 		 *
*			(same as "id" if the outcome is not a repeated measure)                      	     *
*-------------------------------------------------*												 *
* Optional parameters or parameters with defaults *												 *
*-------------------------------------------------*----------------------------------------------*
* covs (optional) : extra covariates (Equivalent to adding them in model_noG_noE)				 *
* time (optional) : time variable for repeated measure outcome   								 *
* start_genes (optional) : starting points for genetic score (must be same length as "genes")    *
* start_env (optional) : starting points for environmental score (must be same length as "env")  *
* eps : threshold for convergence (.01 for quick batch simulations, .0001 for accurate results)  *
* maxiter : Maximum number of iterations (Default=50)											 *
* where (optional) : where                                                                       *
* n_extreme : Number of extreme residuals to show (Deafult=10)                                   *
* repeated : If 1 then the outcome is a repeated measure (Default=0)                             *
* repeated_type : covariance type for repeated measure (Default= un)                             *
* random_vars (optional) : variables of random effects                                           *
* random_sub (optional) : subject for randon effects                                             *
* random_type : covariance type for random effects (Default= vc)                                 *
* dist : Outcome distribution (Default= normal)                                                  *
* link : GLM link (Default=identity)                                                             *
* method : Optimization method (Default=MSPL)                                                    *
* stop_short : If 1 then stop after 5 iterations, used for debugging (Default=0)                 *
* no_log : If 1 then remove the log for the duration of the function (Default=1)				 *
*			Warning, if you stop the macro abruptly, your log will be gone, to get it back do :	 *
*			proc printto; 																		 *
*			run;																				 *
*------------------------------------------------------------------------------------------------*;
*----------*
* Examples *
*----------*
* two-way model : Intercept + g + e + ge + gender + ses
* Four genes g1, g2, g3, g4. Three environments e1, e2, e3. Two gene by gene interaction g1_g2, g1_g4. 
* Continuous outcome y at 4 time-points time="3M", "6M", "18M" or "36M".
* Using PROC MIXED.
*
* LEGIT_cv(data, y, g1 g2 g3 g4 g1_g2 g1_g4, e1 e2 e3, mixed, Intercept, Intercept, Intercept, Intercept, gender ses, id=id, true_id=true_id)			
*
* three-way model : Intercept + g + e + z + ge + gz + ez + gez + gender + ses
* Four genes g1, g2, g3, g4. Three environments e1, e2, e3. Two gene by gene interaction g1_g2, g1_g4. 
* Binary outcome y at one time-point with logit link. 
* Using PROC LOGISTIC.
*
* LEGIT_cv(data, y, g1 g2 g3 g4 g1_g2 g1_g4, e1 e2 e3, logistic, Intercept z, Intercept z, Intercept z, Intercept z, gender ses, true_id=id, id=id, dist=binomial, link=logit)

%macro LEGIT_cv(data, outcome, genes, env, proc_used, model_noG_noE=, model_E_noG=, model_G_noE=, model_G_E=, covs=, id=PSCID, true_id=true_id, time=time, start_env=, start_genes=, eps = .01, maxiter = 50, where=, n_extreme=10, repeated=0, repeated_type=un, random_vars=, random_sub=, random_type=vc, link=identity, dist=normal, method=MSPL, stop_short=0, no_log=1);
	%if &no_log eq 1 %then %do;
		proc printto log="nul:"; 
		run;
	%end;
	%let _timer_start = %sysfunc(datetime());
	*Get unique id into macro variable to loop trough;
	data compress&data;
		set &data;
 		comp&id = compress(&id,'-');
	run;
	proc sort data=compress&data;
		  by &true_id;
	run;
	proc iml;
		use compress&data;
		read all;
		uni = unique(comp&id);
		uni_char = rowcat(uni);
		call symputx("mylist", uni_char);
	quit;
	* Dataset to storeweights of all folds;
	data G_weights_all;
	run;
	data E_weights_all;
	run;
	 ods exclude all;
	%let index = 1;
	%let val = %scan(&mylist,&index);
	%do %until(&val eq %nrstr( ));

		* Dataset without observation;
		data without_obs;
			set compress&data;
			if comp&id eq "&val" then &outcome = .;
		run;
		* LEGIT fit;
		%if &proc_used eq glimmix %then %LEGIT_glimmix(data=without_obs, outcome=&outcome, genes=&genes, env=&env, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=comp&id, remove_miss=0, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, print = 0, print_final=0, ods_new=0, where=&where, clear_ods=1, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, link=&link, dist=&dist, method=&method);
		%if &proc_used eq mixed %then %LEGIT_mixed(data=without_obs, outcome=&outcome, genes=&genes, env=&env, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=comp&id, remove_miss=0, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, print = 0, print_final=0, ods_new=0, where=&where, clear_ods=1, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type);
		%if &proc_used eq logistic %then %LEGIT_logistic(data=without_obs, outcome=&outcome, genes=&genes, env=&env, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=comp&id, remove_miss=0, start_genes=&start_genes, start_env=&start_env, eps=&eps, maxiter =&maxiter, print = 0, print_final=0, ods_new=0, where=&where, clear_ods=1, link=&link);
		* To see stats of betas;
		data G_weights;
			set G_weights;
			rename G_weights_OLD = G_weights&index;
		run;
		data G_weights_all;
			merge G_weights_all G_weights;
		run;
		data E_weights;
			set E_weights;
			rename E_weights_OLD = E_weights&index;
		run;
		data E_weights_all;
			merge E_weights_all E_weights;
		run;
		* Running final model without observations and getting predictions;
		%let model_E_noG_N = %sysfunc(countw(&model_E_noG,' '));
		%let model_G_noE_N = %sysfunc(countw(&model_G_noE,' '));
		%let model_G_E_N = %sysfunc(countw(&model_G_E,' '));
		%if &proc_used eq glimmix %then %do;
			proc glimmix data=Step_a method=&method IC=PQ;
				id comp&id &true_id;
				class comp&id &random_vars  
				%if &repeated eq 1 %then %do;
					&time 
				%end;
				;
				model &outcome(DESC) = &model_noG_noE 
				%do j=1 %to &model_E_noG_N;
					E*%scan(&model_E_noG, &j,' ') 
				%end;
				%do j=1 %to &model_G_noE_N;
					G*%scan(&model_G_noE, &j,' ') 
				%end;
				%do j=1 %to &model_G_E_N;
					E*G*%scan(&model_G_E, &j,' ') 
				%end;
				&covs / noint;
				%if &repeated eq 1 %then %do;
					random &time / subject=comp&id type=&repeated_type residual; 
				%end;
				%else %do;
					random Intercept / subject=comp&id residual;
				%end;
				%if &random_vars ne %then %do;
					random &random_vars / 
					%if &random_sub ne %then %do;
						subject = &random_sub 
					%end;
					type = &random_type;
				%end;
				output out=pred&index pred(noblup ilink)=pred variance(noblup ilink)=var;
			run;
		%end;
		%if &proc_used eq mixed %then %do;
			proc mixed data=Step_a method=ML;
				class comp&id &random_vars  
				%if &repeated eq 1 %then %do;
					&time 
				%end;
				;
				model &outcome = &model_noG_noE 
				%do j=1 %to &model_E_noG_N;
					E*%scan(&model_E_noG, &j,' ') 
				%end;
				%do j=1 %to &model_G_noE_N;
					G*%scan(&model_G_noE, &j,' ') 
				%end;
				%do j=1 %to &model_G_E_N;
					E*G*%scan(&model_G_E, &j,' ') 
				%end;
				&covs / outpm=pred&index noint;
				%if &repeated eq 1 %then %do;
					repeated &time / subject = comp&id type = &repeated_type; 
				%end;
				%if &random_vars ne %then %do;
					random &random_vars / 
					%if &random_sub ne %then %do;
						subject = &random_sub 
					%end;
					type = &random_type;
				%end;
			run;
		%end;
		%if &proc_used eq logistic %then %do;
			proc logistic data=Step_a DESC;
				class comp&id;
				model &outcome = &model_noG_noE 
				%do j=1 %to &model_E_noG_N;
					E*%scan(&model_E_noG, &j,' ') 
				%end;
				%do j=1 %to &model_G_noE_N;
					G*%scan(&model_G_noE, &j,' ') 
				%end;
				%do j=1 %to &model_G_E_N;
					E*G*%scan(&model_G_E, &j,' ') 
				%end;
				&covs / noint link=&link;
				output out=pred&index PRED=pred;
			run;
		%end;
		* Calculating cv predictions;
		data pred&index;
			set pred&index;
			if comp&id eq "&val"; /*only keep predictons from "left out" observations for that fold */
			keep &true_id pred 
			%if &proc_used eq glimmix %then %do;
				var 
			%end;
			;
		run;
		proc sort data=pred&index out=pred&index;
			by &true_id;
		run;

		%let index=%eval(&index+1);
		%let val=%scan(&mylist,&index);
		* For testing purposes, when stop_short=1, stop after 5 iterations;
		%if &stop_short eq 1 and &index > 5 %then %let val=;
	%end;

	* Merging cv predictions with main dataset;
	data final;
		merge compress&data 
		%do i=1 %to (&index-1);
			pred&i 
		%end;
		;
		by &true_id;
	run;

 	ods exclude none;

	* AUC and ROC curve;
	options minoperator mlogic;
	%if &proc_used eq glimmix and &dist in bin,beta,binary,b,binomial,multinomial,mult,multi, %then %do;
		title1 "Leave-one-out cross-validated ROC Curve";
		proc logistic data=final DESC plots(only)=roc;
			model &outcome =  / outroc=rocstats;
			roc pred=pred;
        	roccontrast;
		run;
		data check;
			set rocstats;
			_SPECIF_ = (1 - _1MSPEC_);
			J = _SENSIT_ + _SPECIF_ - 1;
			D= Sqrt((1-_SENSIT_)**2 + (1-_SPECIF_)**2);
		run;
		proc sql noprint;
			create table cutoff as
			select _PROB_ , J, _SENSIT_, _SPECIF_, _FALPOS_, _FALNEG_, _POS_, _NEG_
			from check
			having J = max(J);
		run;
		proc sql noprint;
			create table cutoff1 as
			select _PROB_ , D, _SENSIT_, _SPECIF_, _FALPOS_, _FALNEG_, _POS_, _NEG_
			from check
			having D = min(D);
		run;
		proc print data = cutoff;
			title1 'Best choice of cutoff based on Youden Index';
			var _PROB_ J _SENSIT_ _SPECIF_ _FALPOS_ _FALNEG_ _POS_ _NEG_;
		run;
		proc print data = cutoff1;
			title1 'Best choice of cutoff based on Euclidean Distance';
			var _PROB_ D _SENSIT_ _SPECIF_ _FALPOS_ _FALNEG_ _POS_ _NEG_;
		run;
		title1;
	%end;
	%if &proc_used eq logistic %then %do;
		title1 "Leave-one-out cross-validated ROC Curve";
		proc logistic data=final DESC;
			model &outcome = / outroc=rocstats;
			roc pred=pred;
        	roccontrast;
		run;
		data check;
			set rocstats;
			_SPECIF_ = (1 - _1MSPEC_);
			J = _SENSIT_ + _SPECIF_ - 1;
			D= Sqrt((1-_SENSIT_)**2 + (1-_SPECIF_)**2);
		run;
		proc sql noprint;
			create table cutoff as
			select _PROB_ , J, _SENSIT_, _SPECIF_, _FALPOS_, _FALNEG_, _POS_, _NEG_
			from check
			having J = max(J);
		run;
		proc sql noprint;
			create table cutoff1 as
			select _PROB_ , D, _SENSIT_, _SPECIF_, _FALPOS_, _FALNEG_, _POS_, _NEG_
			from check
			having D = min(D);
		run;
		proc print data = cutoff;
			title1 'Best choice of cutoff based on Youden Index';
			var _PROB_ J _SENSIT_ _SPECIF_ _FALPOS_ _FALNEG_ _POS_ _NEG_;
		run;
		proc print data = cutoff1;
			title1 'Best choice of cutoff based on Euclidean Distance';
			var _PROB_ D _SENSIT_ _SPECIF_ _FALPOS_ _FALNEG_ _POS_ _NEG_;
		run;
		title1;
	%end

  	* Print extreme weights; 
  	proc transpose data=G_weights_all out=G_weights_all_t prefix=m;
		run;
  	proc transpose data=E_weights_all out=E_weights_all_t prefix=e;
		run;
  	proc means data=G_weights_all_t min P5 P50 P95 max;
  	run;
  	proc means data=E_weights_all_t min P5 P50 P95 max;
  	run;

	* Residuals;
	data final;
		set final;
		Resid = pred-&outcome;
		%if &proc_used eq glimmix %then %do;
			Pearson_Resid = (pred-&outcome)/sqrt(var);
		%end;
		%if &proc_used eq logistic %then %do;
			Pearson_Resid = (pred-&outcome)/sqrt(pred*(1-pred));
		%end;
	run;

	* R2 creation;
	proc iml;
		use final;
			read all;
			SSTotal = ssq(&outcome-mean(&outcome));
			R2_leave_one_out = 1 - ssq(Resid)/SSTotal;
		close final;
	create R2_cv var {R2_leave_one_out};
	append;
	print R2_leave_one_out;
	quit;

	* Leave-one-out residuals for outlier inspection;
	title1 "Leave-one-out residuals for outlier inspection";
	%if &proc_used eq glimmix or &proc_used eq logistic %then %do;
		proc standard data = final mean = 0 std = 1 out = final_std;
			var Pearson_Resid;
		  run;
		 proc univariate data=final_std NEXTROBS=&n_extreme;
			id comp&id 
			%if &repeated eq 1 %then %do;
				&time 
			%end;
			;
			var Pearson_Resid;
		run;
		proc means data=final_std min P5 P50 P95 max mean var;
			var Pearson_Resid;
		run;
	%end;
	%if &proc_used eq mixed %then %do;
		proc standard data = final mean = 0 std = 1 out = final_std;
			var Resid;
		  run;
		 proc univariate data=final_std NEXTROBS=&n_extreme;
			id comp&id 
			%if &repeated eq 1 %then %do;
				&time 
			%end;
			;
			var Resid;
		run;
		proc means data=final_std min P5 P50 P95 max mean var;
			var Resid;
		run;
	%end;
	title1;

	%if &no_log eq 1 %then %do;
		proc printto; 
		run;
	%end;
	data _null_;
		dur = datetime() - &_timer_start;
		put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
	run;
  	dm odsresults "clear";

%mend;




*--------------------------------------------------------------------------------*
* Macro #5) LEGIT_search: Add genes or environments, one at a time, to a GxE     *
*		model. Output models for which the added variable had a p-value smaller  *
* 		than the threshold and lower AIC. 										 *
*--------------------------------------------------------------------------------*
* Parameters necessary to run *
*-----------------------------*---------------------------------------------------------------------*
* Data : dataset to be used                                                                         *
* outcome : outcome variable                                                                        *
* proc_used : glimmix = PROC GLIMMIX, logistic = PROC LOGISTIC, mixed = PROC MIXED                  *
* add_genes : If 1 then add genes, otherwise add environments (Default=1)                           *
* p_threshold : Only show models for which p-value < p_threshold (Default=.10)                      *
* 				Can be set to 1 if you don't want to look at p-values but only at AIC               *
* AIC : If 1 then only show models for which AIC decreased (Default=1)                              *
* genes_original : genes variables inside genetic score G 											*
*		(can be any sort of variable, doesn't even have to be genetic)						        *
* genes_extra : additionnal genes variables to be included, one at a time, inside genetic score G.	*
*		If genes_extra is not empty, env_extra should be empty.										*
* env_original : environment variables inside environmental score E								    *
*		(can be any sort of variable, doesn't even have to be environmental)					    *
* env_extra : additionnal environmental variables to be included, one at a time, inside             *
*		environmental score E. If env_extra is not empty, genes_extra should be empty 				*
* model_noG_noE : model part without E and without G (See model specification above) 			    *
* model_E_noG : model part with E and without G (See model specification above)                     *
* model_G_noE : model part without E and with G (See model specification above)                     *
* model_G_E : model part with E and with G (See model specification above)                          *
* id : ID of individual (if no ID, just create it using id=_n_ in a data step)                      *
*-------------------------------------------------*												    *
* Optional parameters or parameters with defaults *												    *
*-------------------------------------------------*-------------------------------------------------*
* covs (optional) : covariates (Equivalent to adding them in model_noG_noE)			        	    *
* time (optional) : time variable for repeated measure outcome   								    *
* start_genes (optional) : starting points for genetic score 										*
* (must be same length as "genes_original") 														*
* start_env (optional) : starting points for environmental score 									*
* (must be same length as "env_original") 															*
* eps : threshold for convergence (.01 for quick batch simulations, .0001 for accurate results)     *
* maxiter : Maximum number of iterations														    *
* repeated : If 1 then the outcome is a repeated measure (Default=0)                                *
* repeated_type : covariance type for repeated measure (Default= un)                                *
* random_vars (optional) : variables of random effects                                              *
* random_sub (optional) : subject for randon effects                                                *
* random_type : covariance type for random effects (Default= vc)                                    *
* where (optional) : where                                                                          *
* dist : Outcome distribution (Default= normal)                                                     *
* link : GLM link (Default=identity)                                                                *
* method : Optimization method (Default=MSPL)                                                       *
* ods_new : If 1 then close current ods output and start a new one (Default=1)                      *
* same_sample : If 1 then all extra variables in "genes_extra" or "env_extra" have the exact same	*
*		observations, meaning that all missing observations are the same for all extra variable. 	*
*		Otherwise, we have to refit the model without the extra variable at every iteration 		*
*		because AIC can only be compared on the same sample. (Default=0)							*
* no_log : If 1 then remove the log for the duration of the function (Default=1)					*
*			Warning, if you stop the macro abruptly, your log will be gone, to get it back do :	 	*
*			proc printto; 																		 	*
*			run;																				 	*
*---------------------------------------------------------------------------------------------------*;
*----------*
* Examples *
*----------*
* two-way model : Intercept + g + e + ge + gender + ses
* Four genes g1, g2, g3, g4. Three environments e1, e2, e3. Two gene by gene interaction g1_g2, g1_g4. 
* Continuous outcome y at 4 time-points time="3M", "6M", "18M" or "36M".
* Using PROC MIXED, searching for gene interactions.
*
* LEGIT_search(data, y, mixed, add_genes=1, genes_original=g1 g2 g3 g4, genes_extra=g1_g2 g1_g4, env_original=e1 e2 e3, model_noG_noE=Intercept, model_E_noG=Intercept, model_G_noE=Intercept, model_G_E=Intercept, covs=gender ses, id=id, repeated=1, time=time)
* 
* three-way model : Intercept + g + e + z + ge + gz + ez + gez + gender + ses
* Four genes g1, g2, g3, g4. Three environments e1, e2, e3. Two gene by gene interaction g1_g2, g1_g4. 
* Binary outcome y at one time-point with logit link.
* Using PROC LOGISTIC, searching for environments.
*
* LEGIT_search(data, y, logistic, add_genes=0, genes_original=g1 g2 g3 g4 g1_g2 g1_g4, env_original=e1, env_extra=e2 e3, model_noG_noE=Intercept z, model_E_noG=Intercept z, model_G_noE=Intercept z, model_G_E=Intercept z, covs=gender ses, id=id, repeated=0, dist=binomial, link=logit)


%macro LEGIT_search(data, outcome, proc_used, add_genes=1, p_threshold=.10, AIC=1, genes_original=, genes_extra=,env_original=,env_extra=, model_noG_noE=, model_E_noG=, model_G_noE=, model_G_E=, covs=, id=PSCID, time=, start_genes=, start_env=, eps = .001, maxiter = 50, where=, repeated=0, repeated_type=un, random_vars=, random_sub=, random_type=vc, link=identity, dist=normal, method=MSPL, ods_new=1, same_sample=0, no_log=1);
	%if &no_log eq 1 %then %do;
		proc printto log="nul:"; 
		run;
	%end;
	%let _timer_start = %sysfunc(datetime());
	%if &ods_new eq 1 %then %do;
		ods html close;
		ods html;
	%end;
	%if &genes_extra ne %then %let extra_genes_N = %sysfunc(countw(&genes_extra));
	%let original_genes_N = %sysfunc(countw(&genes_original));
	%if &env_extra ne %then %let extra_env_N = %sysfunc(countw(&env_extra));
	%let original_env_N = %sysfunc(countw(&env_original));
	%let model_noG_noE_N = %sysfunc(countw(&model_noG_noE,' '));
	%let model_E_noG_N = %sysfunc(countw(&model_E_noG,' '));
	%let model_G_noE_N = %sysfunc(countw(&model_G_noE,' '));
	%let model_G_E_N = %sysfunc(countw(&model_G_E,' '));

	%if &same_sample eq 1 %then %do;
		* Missing data is the same in all extra variables, we only need to run the original model once;
		* Run model one time without extra genes/env to get the baseline AIC;
		data data_nomiss;
			set &data;
			where %scan(&genes_extra,1) ne .;
		run;
		ods exclude all;
		%if &proc_used eq glimmix %then %do;
			%LEGIT_glimmix(data=data_nomiss, outcome=&outcome, genes=&genes_original, env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, link=&link, dist=&dist, method=&method, print_final=0);
		%end;
		%if &proc_used eq logistic %then %do;
			%LEGIT_logistic(data=data_nomiss, outcome=&outcome, genes=&genes_original, env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, link=&link, print_final=0)
		%end;
		%if &proc_used eq mixed %then %do;
			%LEGIT_mixed(data=data_nomiss, outcome=&outcome, genes=&genes_original, env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, print_final=0);
		%end;
		ods exclude none;
		data AIC_old;
			set AIC;
			rename AIC = AIC_old;
		run;
	%end;

	%if &add_genes eq 1 %then %do;
		%do iter=1 %to &extra_genes_N;
			%Let sign = 0;
			%Let AIC_lower = 0;
			ods exclude all;
			%if &proc_used eq glimmix %then %do;
				%LEGIT_glimmix(data=&data, outcome=&outcome, genes=&genes_original %scan(&genes_extra,&iter), env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, link=&link, dist=&dist, method=&method, print_final=0);
			%end;
			%if &proc_used eq logistic %then %do;
				%LEGIT_logistic(data=&data, outcome=&outcome, genes=&genes_original %scan(&genes_extra,&iter), env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, link=&link, print_final=0)
			%end;
			%if &proc_used eq mixed %then %do;
				%LEGIT_mixed(data=&data, outcome=&outcome, genes=&genes_original %scan(&genes_extra,&iter), env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, print_final=0);
			%end;
			* keep best weights;
			Data G_weights_best;
				set G_weights;
			run;
			Data E_weights_best;
				set E_weights;
			run;
			ods exclude none;
			* Verify is p-value < p_threshold;
			data Solutionf_G;
				set Solutionf_G;
				%if &proc_used eq logistic %then %do;
					if (_N_ eq %eval(&original_genes_N+1) and ProbChiSq < &p_threshold) then call symput('sign',1);
				%end;
				%else %do;
					if (_N_ eq %eval(&original_genes_N+1) and Probt < &p_threshold) then call symput('sign',1);
				%end;
			run;
			%if &sign eq 1 %then %do;
				%if &AIC eq 1 %then %do;
					data AIC_new;
						rename AIC = AIC_new;
						set AIC;
					run;
					%if &same_sample eq 0 %then %do;
						* Run model one time without extra genes/env to get the baseline AIC;
						data data_nomiss;
							set &data;
							where %scan(&genes_extra,&iter) ne .;
						run;
						ods exclude all;
						%if &proc_used eq glimmix %then %do;
							%LEGIT_glimmix(data=data_nomiss, outcome=&outcome, genes=&genes_original, env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, link=&link, dist=&dist, method=&method, print_final=0);
						%end;
						%if &proc_used eq logistic %then %do;
							%LEGIT_logistic(data=data_nomiss, outcome=&outcome, genes=&genes_original, env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, link=&link, print_final=0)
						%end;
						%if &proc_used eq mixed %then %do;
							%LEGIT_mixed(data=data_nomiss, outcome=&outcome, genes=&genes_original, env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, print_final=0);
						%end;
						ods exclude none;
						* Check if AIC is better;
						data AIC_old;
							set AIC;
							rename AIC = AIC_old;
						run;
					%end;
					data AIC_new;
						merge AIC_new AIC_old;
					run;
					data AIC_new;
						set AIC_new;
						if AIC_new < AIC_old then call symput('AIC_lower',1);
					run;
					%if &AIC_lower eq 1 %then %do;
						proc print data=AIC_new;
							var AIC_new AIC_old;
						run;
					%end;
				%end;
				%if (&AIC eq 1 and &AIC_lower eq 1) or (&AIC eq 0) %then %do;
					* Transform best weights into macro variables;
					proc iml;
						use G_weights_best;
						read all;
						use E_weights_best;
						read all;
						G_weights = rowcat(" " + shape(char(G_WEIGHTS_OLD),1));
						E_weights = rowcat(" " + shape(char(E_WEIGHTS_OLD),1));
						call symputx("genes_start_final",G_weights);
						call symputx("env_start_final",E_weights); 
					quit;
					* Refits LEGIT model for one iteration at final estimates and output it;
					%if &proc_used eq glimmix %then %do;
						%LEGIT_glimmix(data=&data, outcome=&outcome, genes=&genes_original %scan(&genes_extra,&iter), env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&genes_start_final, start_env=&env_start_final, time=&time, eps=&eps, maxiter = 1, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, link=&link, dist=&dist, method=&method, print_final=1);
					%end;
					%if &proc_used eq logistic %then %do;
						%LEGIT_logistic(data=&data, outcome=&outcome, genes=&genes_original %scan(&genes_extra,&iter), env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&genes_start_final, start_env=&env_start_final, eps=&eps, maxiter = 1, ods_new=0, where=&where, clear_ods=0, link=&link, print_final=1)
					%end;
					%if &proc_used eq mixed %then %do;
						%LEGIT_mixed(data=&data, outcome=&outcome, genes=&genes_original %scan(&genes_extra,&iter), env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&genes_start_final, start_env=&env_start_final, time=&time, eps=&eps, maxiter = 1, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, print_final=1);
					%end;
				%end;
			%end;
		%end;
	%end;
	%else %do;
		%do iter=1 %to &extra_env_N;
			%Let sign = 0;
			%Let AIC_lower = 0;
			ods exclude all;
			%if &proc_used eq glimmix %then %do;
				%LEGIT_glimmix(data=&data, outcome=&outcome, genes=&genes_original, env=&env_original %scan(&env_extra,&iter), model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, link=&link, dist=&dist, method=&method, print_final=0);
			%end;
			%if &proc_used eq logistic %then %do;
				%LEGIT_logistic(data=&data, outcome=&outcome, genes=&genes_original, env=&env_original %scan(&env_extra,&iter), model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, link=&link, print_final=0)
			%end;
			%if &proc_used eq mixed %then %do;
				%LEGIT_mixed(data=&data, outcome=&outcome, genes=&genes_original, env=&env_original %scan(&env_extra,&iter), model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, print_final=0);
			%end;
			* keep best weights;
			Data G_weights_best;
				set G_weights;
			run;
			Data E_weights_best;
				set E_weights;
			run;
			ods exclude none;
			data Solutionf_E;
				set Solutionf_E;
				%if &proc_used eq logistic %then %do;
					if (_N_ eq %eval(&original_env_N+1) and ProbChiSq < &p_threshold) then call symput('sign',1);
				%end;
				%else %do;
					if (_N_ eq %eval(&original_env_N+1) and Probt < &p_threshold) then call symput('sign',1);
				%end;
			run;
			* If it is significant rerun model and print it;
			%if &sign eq 1 %then %do;
				data AIC_new;
					set AIC;
					rename AIC = AIC_new;
				run;
				%if &same_sample eq 0 %then %do;
					* Run model one time without extra genes/env to get the baseline AIC;
					data data_nomiss;
						set &data;
						where %scan(&env_extra,&iter) ne .;
					run;
					ods exclude all;
					%if &proc_used eq glimmix %then %do;
						%LEGIT_glimmix(data=data_nomiss, outcome=&outcome, genes=&genes_original, env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, link=&link, dist=&dist, method=&method, print_final=0);
					%end;
					%if &proc_used eq logistic %then %do;
						%LEGIT_logistic(data=data_nomiss, outcome=&outcome, genes=&genes_original, env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, link=&link, print_final=0)
					%end;
					%if &proc_used eq mixed %then %do;
						%LEGIT_mixed(data=data_nomiss, outcome=&outcome, genes=&genes_original, env=&env_original, model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&start_genes, start_env=&start_env, time=&time, eps=&eps, maxiter =&maxiter, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, print_final=0);
					%end;
					ods exclude none;
					* Check if AIC is better;
					data AIC_old;
						set AIC;
						rename AIC = AIC_old;
					run;
				%end;
				data AIC_new;
					merge AIC_new AIC_old;
				run;
				data AIC_new;
					set AIC_new;
					if AIC_new < AIC_old then call symput('AIC_lower',1);
				run;
				proc print data=AIC_new;
					var AIC_new AIC_old;
				run;
				%if (&AIC eq 1 and &AIC_lower eq 1) or (&AIC eq 0) %then %do;
					* Transform best weights into macro variables;
					proc iml;
						use G_weights_best;
						read all;
						use E_weights_best;
						read all;
						G_weights = rowcat(" " + shape(char(G_WEIGHTS_OLD),1));
						E_weights = rowcat(" " + shape(char(E_WEIGHTS_OLD),1));
						call symputx("genes_start_final",G_weights);
						call symputx("env_start_final",E_weights); 
					quit;
					* Refits LEGIT model and output it;
					%if &proc_used eq glimmix %then %do;
						%LEGIT_glimmix(data=&data, outcome=&outcome, genes=&genes_original, env=&env_original %scan(&env_extra,&iter), model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&genes_start_final, start_env=&env_start_final, time=&time, eps=&eps, maxiter =1, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, link=&link, dist=&dist, method=&method, print_final=1);
					%end;
					%if &proc_used eq logistic %then %do;
						%LEGIT_logistic(data=&data, outcome=&outcome, genes=&genes_original, env=&env_original %scan(&env_extra,&iter), model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&genes_start_final, start_env=&env_start_final, eps=&eps, maxiter =1, ods_new=0, where=&where, clear_ods=0, link=&link, print_final=1)
					%end;
					%if &proc_used eq mixed %then %do;
						%LEGIT_mixed(data=&data, outcome=&outcome, genes=&genes_original, env=&env_original %scan(&env_extra,&iter), model_noG_noE=&model_noG_noE, model_E_noG=&model_E_noG, model_G_noE=&model_G_noE, model_G_E=&model_G_E, covs=&covs, id=&id, start_genes=&genes_start_final, start_env=&env_start_final, time=&time, eps=&eps, maxiter =1, ods_new=0, where=&where, clear_ods=0, repeated=&repeated, repeated_type=&repeated_type, random_vars=&random_vars, random_sub=&random_sub, random_type=&random_type, print_final=1);
					%end;
				%end;
			%end;
		%end;
	%end;
	%if &no_log eq 1 %then %do;
		proc printto; 
		run;
	%end;

	data _null_;
		dur = datetime() - &_timer_start;
		put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
	run;
	dm odsresults "clear";
%mend;