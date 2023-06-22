/*
AUTHOR:
(1)	William Liu (刘威廉)

ACKNOWLEDGEMENTS:
(1)	My -ivcloglog- command is inspired by Dr. Enrique Pinzon's (2020; Stata Conference) presentation on implementing control functions in Stata.
	His presentation has helped teach me how to write a -gmm- moment evaluator program (i.e., Stata subroutine).
	
	In addition, Dr. Pinzon shared the code for his -cfunction- command with me.
	Whilst my -ivcloglog- code diverges greatly from the -cfunction- code, I have found it extremely useful for learning how to write a Stata command.
	
	Many thanks to Dr. Enrique Pinzon, whose assistance has been invaluable for the development of -ivcloglog-.
	
	Link:
		https://www.stata.com/meeting/uk20/slides/UK20_Pinzon.pdf

GUIDE:
(1)	Syntax:
		ivcloglog outcomevar controlvars ..., ENDOgenous(endovars = 1ststageregressorvars, noCONstant) VHATname(vhatname) [noCONstant order(...) vce() NOGENerate]

(2)	In this implementation, the first-stage regressors are assumed to be the same for all endogenous variables.
(3)	-vce(robust)- is the default vce option.
(4)	-vhatname()- specifies the name for the vhat(s), the residuals of the first stage(s). ""
(5)	-order()- specifies the order of the control function polynomial and is 1 by default.
(6)	-nogenerate- is used to stop the powers of the vhat(s), which comprise the control functions, from being added to the dataset.

COMMENTS:
(1)	The first-stage standard errors from the -gmm- call will be larger than those from the -regress- call.
	This is because they are adjusted for the fact that we are also doing second-stage estimates.
	However, this is unnecessary because the first-stage estimates are unaffected by the second-stage estimates.
	The SEs from -regress- should be preferentially used.
*/

program ivcloglog, eclass byable(onecall)
	version 16
	
	if _by() {
		local BY `"by `_byvars'`_byrc0':"'
	}

	`BY' _vce_parserun ivcloglog, noeqlist jkopts(eclass): `0'
		
	`BY' Estimate `0'
end

program define Estimate, eclass byable(recall)
	version 16
	syntax varlist(fv numeric ts) [if] [in]								///
								  [fweight iweight pweight aweight],	///
								  VHATname(string)						///
								  [										///
								  noCONstant                            ///
								  _get_diopts							///
								  order(integer 1)						///
								  vce(string)							///
								  NOGENerate							///
								  *										///
								  ]
	_get_diopts diopts other, `options'
	
	if ("`weight'" != "") {
		local wgt [`weight'`exp']
	}
    
	gettoken depvar indepvars : varlist
	_fv_check_depvar `depvar', k(1)				// Stops depvar being a factor variable
	
	* Parse -endogenous(...)- option
	local nocon2 `constant'						// For better readability
	parse_remaining_opt, `options'
	local options `s(options)'
	local num_eqs =`s(eq)'
	local allinst `s(allinst)'
	local endovars `s(endovars)'
	local nocon1 `s(nocons1)'
	
	marksample touse
	
	* Check for perfect predictors of the outcome and collinearity amongst the regressors
	iden_check, lhs(`depvar') exog(`indepvars') endo(`endovars') allinst(`allinst') touse(`touse') wgt(`wgt') asis(`asis') nocon1(`nocon1') nocon2(`nocon2')
	local indepvars `r(exog)'
	local allinst `r(allinst)'
	
	* Formatting
	tempname A1 A2 B F
	forval i = 1/`num_eqs' {
		forval p = 1/`order' {
			local resid`i'_`p' `vhatname'`i'_`p'		// For better readability
			qui gen double `resid`i'_`p'' = .			// Previously, for some reason, initializing the higher powers with missing values didn't work when not using -from(...)-
			local resid`i'_powers = "`resid`i'_powers' `resid`i'_`p''"
		}
		
		local endovar`i' : word `i' of `endovars'
		init_values_1 `endovar`i'' `allinst' if `touse' `wgt', `nocon1' predictvars(`resid`i'_powers') order(`order')		// Store powers of vhat for endo. variable `i' in local macro `resid`i'_powers'
		matrix `B' = r(A)
		matrix `A1' = nullmat(`A1'), `B'		// -nullmat()- relaxes the restriction that `A' must exist
		local yparam1 "`r(yparam)'"
		local yparam1_combined "`yparam1_combined' `yparam1'"
		
		forval p = 1/`order' {
			local resid_powers_combo "`resid_powers_combo' `resid`i'_`p''"		// Powers of vhat from all endo. variables
		}
		
		local inst_opt1 "`inst_opt1' instruments(`endovar`i'':`allinst', `nocon1')"
	}
	
	local indepvars "`indepvars' `resid_powers_combo'"
	init_values_2 `depvar' `endovars' `indepvars' if `touse' `wgt', `nocon2'
	local yparam2 "`r(yparam)'"

	matrix `A2' = r(A)
	matrix `F' = `A1', `A2'
	
	local inst_opt2 "instruments(`depvar': `endovars' `indepvars', `nocon2')"
	
	* Run GMM
	// The undocumented option -iter(0)- prevents numerical iterations -- we only want the sandwich formula.
	// -haslfderivatives- differs from -hasderivatives- in that you don't need to specify each parameter, only the parameter groupings.
	// These are what I call the parameters that form linear combinations.
	// They are specified using the curly braces in the -parameter()- option or through residual equation names.
	gmm gmm_iv_cloglog if `touse', equations(`endovars' `depvar')								///
								   parameters("`yparam1_combined' `yparam2'")					///
								   y1(`depvar') y2(`endovars') 									///		// Wooldridge-style notation for this line; the 1 and 2 on other lines represent stages
								   vhatname(`vhatname') order(`order')							///
								   `inst_opt1' `inst_opt2'										///
								   winitial(identity)											///
								   haslfderivatives onestep iter(0) from(`F')					///
								   vce(`vce')													///
								   valueid("EE criterion") `diopts'
								   
	* Drop powers of vhat if requested with -nogenerate-
	// Yes, -nogenerate- is a misnomer because we do actually generate the powers of vhat.
	// We do this because it's easier than renaming the columns of ("column-restriping") the matrix holding the coefficient estimates,
	// which is necessary to make the powers of vhat be called "`vhatname'`i'_`p'".
	if "`nogenerate'" != "" {
		drop `resid_powers_combo'
	}
end

program parse_remaining_opt, sclass		// Parses -endogenous()- and also checks for remaining options (should be none)
	version 16
	syntax [, ENDOgenous(string) *]
	
	* Check for remaining options (should be none)
	if ("`options'"!="") {
		di as err "option(s) {bf:`options'} not allowed"
		exit 198
	}
	
	* Parse contents of -endogenous()-
	gettoken endogenous temp: endogenous, parse(",")		// Delimiter is stored in `temp' (which is not str.) and is then deleted
	local options = substr("`temp'", 2, .)
	parse_subopt_noconstant, `options'
	local nocons1 `r(constant)'								// Don't directly store `options' in `nocons1' since the former can be abbreviated
	
	gettoken temp_endovars temp: endogenous, parse("=")
	local temp_allinst = substr("`temp'", 2, .)
	
	* Count and parse endogenous variables
	local eq = 0
	local k: list sizeof temp_endovars
	forval i = 1/`k' {
		local ++eq
	    local temp_word : word `i' of `temp_endovars'		// Format is endogenous(endovar_1 ... endovar_k = ...)
		_fv_check_depvar `temp_word', k(1)					// Stops the dependent variables in the auxiliary equations being factor variables
	}
	
	sreturn local eq `eq'
	sreturn local options `options'
	sreturn local allinst "`temp_allinst'"
	sreturn local endovars "`temp_endovars'"
	sreturn local nocons1 "`nocons1'"
end

program define parse_subopt_noconstant, rclass
	version 16
	syntax [anything], [noCONstant * ]
	if ("`options'"!="") {
		di as err "suboption(s) {bf:`options'} not allowed in option {bf:endogenous()}"
		exit 198
	}
	
	return local constant "`constant'"
end

program define init_values_1, rclass		// First stage
	version 16
	syntax [anything][if][in]							///
					 [fweight iweight pweight aweight], ///
					 [									///
					 noCONstant 						///
					 predictvars(string)				///
					 order(integer 1)					///
					 ]
		
		marksample touse 	
		if ("`weight'" != "") {
			local wgt [`weight'`exp']
        }
		
		tempname A
		qui regress `anything' if `touse' `wgt', `constant'
		
		matrix `A' = e(b)		// e(b) from -regress- DOES NOT have the format "depvar: ...", so we must change the format
		reformat_params			// (Renaming the columns is unnecessary because -from()- extracts positionally and ignores names.)
		local yparam "`r(yparam)'"
		
		tempvar prediction
		qui predict double `prediction' if `touse', resid
		forval p = 1/`order' {
			local predictvar_`p' : word `p' of `predictvars'
			qui replace `predictvar_`p'' = `prediction'^`p'
		}
		
		return matrix A = `A'
		return local yparam "`yparam'"
end

program define init_values_2, rclass		// Second stage
	version 16 
	syntax [anything][if][in]							///
					 [fweight iweight pweight aweight], ///
					 [									///
					 noCONstant 						///
					 ]
					 
		marksample touse
		if ("`weight'" != "") {
			local wgt [`weight'`exp']
		}
		
		tempname A
		qui cloglog `anything' if `touse' `wgt', `constant'
		
		matrix `A' = e(b)		// e(b) from -cloglog- DOES have the format "depvar: ..."
		local yparam: colfullnames e(b)

		return matrix A = `A'
		return local yparam "`yparam'"
end 

program define reformat_params, rclass
	version 16
	syntax [anything], [subtract(string)]
	
	local paramy: colfullnames e(b)
	local yvar "`e(depvar)'"
	
	foreach var in `paramy' {
		local yparam "`yparam' `yvar':`var'"
	}

	return local yparam "`yparam'"
end

program define iden_check, rclass		// Based off ivprobit.ado
	syntax, lhs(varname fv ts) [exog(varlist fv ts)] endo(varlist fv ts) allinst(varlist fv ts) touse(name) [wgt(name) asis(name) nocon1(string) nocon2(string)]
	
	* Primary model
	di as text "Checking primary model for perfect predictors of the outcome and regressor collinearity..."
	cap noi _rmcoll `lhs' `endo' `exog' if `touse' `wgt', touse(`touse') logit expand `asis' `nocon2' noskipline		// Changed: putting the exogenous covariates at the end instead prioritizes them for omission
	if c(rc) {							// Equivalent to using _rc ()
		error c(rc)
	}

	* Extract the exogenous controls post-flagging
	local vlist "`r(varlist)'"
	gettoken lhs vlist : vlist
	
	fvexpand `endo'						// In the -ivprobit- code, the -expand- option on -_rmcoll- is actually redundant because the variables were already expanded beforehand.
	local endo_expanded `r(varlist)'	// Not the case here. So, we need to use the correct number of (endogenous) covariates that we get AFTER expanding factor variables!
	
	local n : list sizeof endo_expanded	// In -ivprobit-, this is -local n : list sizeof exog- instead, which uses the number of covariates before -_rmcoll, expand- is called.
	local endo							// Clear this local macro since it's already been used for something else
	forval i = 1/`n' {
		gettoken v vlist : vlist
		local endo `endo' `v'
	}
	local exog : copy local vlist

	* Disallow dropping endogenous variables
	CheckEndogDropped, endog(`endo')
		
	* Auxiliary model
	di as text "Checking auxiliary model for regressor collinearity..."
	cap noi _rmcoll `allinst' if `touse' `wgt', touse(`touse') expand `asis' `nocon1' noskipline		// Note: all instruments are shared between the auxiliary models
	if c(rc) {
		error c(rc)
	}
	
	local vlist "`r(varlist)'"
	
	return local exog "`exog'"
	return local allinst "`vlist'"
end

program define CheckEndogDropped		// Adapted from ivprobit.ado
	syntax, endog(string)

	while "`endog'" != "" {
		gettoken var endog : endog
		_ms_parse_parts `var'			// Parses the matrix stripe `var' and returns the token parts in r()

		if r(omit) {					// From _ms_parse_parts
			di as err "tried to drop an endogenous variable, which is not allowed"
			exit 498
		}
	}
end