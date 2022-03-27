// This is a simple utility program that corrects the equation and variable labels in the 
// coefficient and VCE matrices in the marginal effect estimation after multinomial estimations, 
// such as mprobit and mlogit. 
// Typically, the format of resulting estimation matrices should have alternatives as the equation 
// names and variables as the variable names, which the current version of Stata has backwards. 

// It also allows the user to specify the two alternatives, when the differences are of interest. 
// If the option is left blank, it computes the results from all the outcomes.

// Usage: 
// fixit2[, diff(Trump Clinton)]

// Written by Won-ho Park Jan. 2022.


capture prog drop fixit2
program fixit2, eclass
      syntax [, DIFF (namelist min=2 max=2)]

* Check whether the model is either mlogit of mprobit
local command =  e(est_cmd)
if ("`command'" ~= "mlogit")&("`command'" ~= "mprobit"){
      di "This only works with <margins, dydx(*) post> after mprobit and mlogit."
      error 9022
}
* Check whether the 
local numchoices=e(k_predict)
if (`numchoices'<2){
      di "If you are only predicting one outcome, no need for this program."
      error 9022
}

mat B = e(b)
*mat list B
mat V = e(V)

local varlist: coleq B
local eqnames: colnames B

*mat coleq B = `eqnames'
mat colnames B = `varlist'

*mat list B
local alternatives 
      forvalues i=1/`numchoices'{
            local phrase = e(predict`i'_label)
            tokenize `phrase', parse("==" ")")
            local alternatives "`alternatives'  `3'"
 *           di "`alternatives'"
      }

local xvars= e(xvars)
*di "`xvars'"
local newcollect
foreach x in `alternatives'{
      foreach y in `xvars'{
            local newcollect "`newcollect' `x':`y'"
*            di "`newcollect'"
      }
}
*di "`newcollect'"

local numvar = colsof(B)/`numchoices'
*di "`numvar'"
local collect
      forvalues j = 1/`numvar'{
            local collect "`collect' `alternatives'"
*            di "`j' `collect'"
      }


mat coleq B = `collect'

*mat list B
matselrc B b, c(`newcollect')
*mat list b

*local newfullname: colfullnames B
*di "`newfullname'"
mat roweq V = `collect'
mat rownames V = `varlist'
mat coleq V = `collect'
mat colnames V = `varlist'

foreach alt in `alternatives'{
      matselrc b `alt', c(`alt':)
      ereturn mat `alt'=`alt'
}


matselrc V v, c(`newcollect') r(`newcollect')
mat V=v
*mat list b
*mat list V

ereturn repost b=b V=V, resize
ereturn local alts `alternatives'
ereturn local xvars `xvars'

if "`diff'" ~= ""{
      tokenize `diff'
      local alt1 `1'
      local alt2 `2'
      di "You chose to compute the difference between `alt1' and `alt2'. "

      matrix bMNL = J(1,`numvar',0)
      matrix sMNL = J(1,`numvar',0)
      scalar j = 0
      foreach var in `xvars' {
                  scalar j = j + 1
      *            di "`var'"
                  qui lincom ([`alt1']`var' - [`alt2']`var')/2
                  matrix bMNL[1,j] = r(estimate)
                  matrix sMNL[1,j] = r(se)
      *            mat list sMNL
            }
      mat coleq bMNL = "Diff"
      mat colnames bMNL = `xvars'
      mat coleq sMNL = "Diff"
      mat colnames sMNL = `xvars'   
      mat V = diag(sMNL)*diag(sMNL)
      mat err = e(error)[1, 1..`numvar']
      mat obs = e(_N)[1, 1..`numvar']
      ereturn mat error = err
      ereturn mat _N = obs
      ereturn repost b=bMNL V=V, resize
}

/*    mat b = (e(`alt1')-e(`alt2'))/2
*     mat list b
      mat coleq "Diff"
      matselrc V V_1, c(`alt1':) r(`alt1':)
      matselrc V V_2, c(`alt2':) r(`alt2':)
      matselrc V cov, c(`alt1':) r(`alt2':)
      mat v = (I(`numvar')*V_1 + I(`numvar')*V_2 - 2*cov)/4
      mat V=v
      mat list V 
      * Note this does not work very well. Maybe come back to this later. For now, just use the se vector instead. 
      */
* }


ereturn display
end
