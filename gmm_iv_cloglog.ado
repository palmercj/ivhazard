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

program define gmm_iv_cloglog
	version 16
	syntax varlist if [fweight iweight pweight],	///
					at(name)						///			// This contains the parameter estimates
					y1(varname)						///
					y2(varlist)						///
					vhatname(string)				///
					order(integer)					///
					[								///
					derivatives(varlist)			///			// Must be optional because -gmm- presumably calls the moment evaluator program multiple times but does not request derivatives every time
					]
	
	local num_endo: list sizeof y2
	tokenize `varlist'		// Get the names of the variables storing the residual functions
	
	* Reduced form equations ("first stage")
	forval i = 1/`num_endo' {
		local y2_`i' : word `i' of `y2'
		
		tempvar zp_`i'
		local mfeval_reduced`i' ``i''
		matrix score double `zp_`i'' = `at' `if', eq(#`i')				// Get projection
		qui replace `mfeval_reduced`i'' = `y2_`i'' - `zp_`i'' `if'		// Gets returned as evaluated residual function
	}
	
	* Main equation ("second stage") moment functions except for CFs
	tempvar xb
	local mfeval_main ``=`num_endo' + 1''
	matrix score double `xb' = `at' `if', eq(#`=`num_endo' + 1') 
	qui replace `mfeval_main' = cond(`y1', exp(`xb'-exp(`xb'))/(-expm1(-exp(`xb'))), -exp(`xb')) `if'
	// These functions are "residual functions" (more specifically, the contribution of each observation).
	// Stata multiplies them with the instruments (specified for a given "residual equation") to get the moment functions.
	// Each residual function is used by Stata to create w moment functions, where w is the number of instruments for that residual function.
	// If no instruments are specified for a given residual function (like this section), then it directly becomes a moment function.
	
	* Calculate derivatives only if requested
	if "`derivatives'" == "" {
		exit
	}
	
	* Get the names of the variables storing the derivatives
	forval w = 1/`=(`num_endo' + 1)^2' {
		local d`w' : word `w' of `derivatives'
	}
	
	local c = 0
	
	* First-stage residual function(s)
	forval eq = 1/`num_endo' {
		forval pg = 1/`=`eq' - 1' {		// Note that loops do not run if the second number is less than the first
			local ++c
			replace `d`c'' = 0 `if'
		}
		
		local ++c
		replace `d`c'' = -1 `if'
		
		forval pg = 1/`=`num_endo' - `eq' + 1' {
			local ++c
			replace `d`c'' = 0 `if'
		}
	}
		
	* Second-stage residual function
	tempvar C_i
	gen double `C_i' = cond(`y1', exp(`xb')*(exp(exp(`xb'))-exp(`xb'+exp(`xb'))-1)/(expm1(exp(`xb')))^2, -exp(`xb'))
	
	forval eq = 1/`num_endo' {		// Get the needed rhos (coefficients on vhat powers)
		forval p = 1/`order' {
			scalar rho`eq'_`p' = `at'[1, colnumb(`at', "`y1': `vhatname'`eq'_`p'")]
		}
	}
	
	forval eq = 1/`num_endo' {
		forval p = 2/`order' {
			local temp`eq' "`temp`eq'' + rho`eq'_`p' * `vhatname'`eq'_`=`p' - 1'"
		}
		local resid`eq'_polyderiv "rho`eq'_1 `temp`eq''"
		
		local ++c
		replace `d`c'' = -`C_i' * (`resid`eq'_polyderiv') `if'
	}
	
	local ++c
	replace `d`c'' = `C_i' `if'
	
	// We must supply Stata with the derivatives of the residual functions, so there are
	//     <no. of residual functions>*<no. of parameter groups>
	// derivatives to supply, NOT
	//     <no. of moment functions>*<no. of parameter groups>
end