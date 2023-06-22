* This do file demonstrates the use of the command ivcloglog.ado to estimate the
* control-function proportional hazards model found in the paper
* "An IV Hazard Model of Loan Default with an Application to Subprime Mortgage Cohorts"
* by Christopher Palmer, cjpalmer@mit.edu

* open a sample dataset
webuse laborsup, clear

* pretend that the failure variable is fem_work
* and that the key endogenous explanatory variable is fem_educ
* the excluded instrument is kids
* male_educ is an included exogenous control

* for demonstration purposes, first estimate the control function estimator by hand

* estimate first stage residuals
reg fem_educ male_educ kids
predict vhat, residuals

* estimate the second stage controlling for the control function, aka the first-stage residuals
cloglog fem_work male_educ fem_educ vhat

* the problem with the above coefficients is that the standard errors aren't adjusted for
* the residuals being a generated regressor with estimation error
* as discussed in the paper, one solution is to use the GMM variance formula
* with the moment conditions for both stages stacked, similar to Newey (1987)
* this is what the ivcloglog.ado command does, using Stata's built-in GMM command
* and feeding in the start values from the "by hand" estimation since optimization on the GMM
* objective function is extremely inefficient.

* note that the RHS endogenous variable is specified only in the endogenous option
* and not in the varlist in the initial command call.
* included controls need to be manualy included in both the first and second stages
* the nogenerate option specifies that the first stage residuals aren't to be saved after
ivcloglog fem_work male_educ, endog(fem_educ = kids male_educ) vhatname(vhat) nogenerate

* you can also do robust clustered standard errors with vce(robust) or vce(cluster clustvar) as options

* you can see that the coefficients are the same between the second stage cloglog command above
* and the ivcloglog command, as well as between the first stage regress command and the ivcloglog command
* however, the standard errors are much larger once correcting for the uncertainty in the residuals

* end of do file
