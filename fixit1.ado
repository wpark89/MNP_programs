// This is a simple utility program that can be used after estimation of marginal effects in
// multinomial probit models with correlated errors, such as cmmprobit. It requires the active
// estimates results are from "margins, dydx(*)" after cmmprobit. It does the following: 
// 1) 'trims' the dimensions and corrects the names of coefficients and the VCE matrix 
// after a cmmprobit estimation with a "fake" unity alternative-specific variable is used. The name of 
// this "fake" unity variables is required. 
// 2) recomputes the differences of coefficients on two alternatives when the differences are of interest. 
// If the option is left blank, it computes the results from all the outcomes.
// Usage: 
// fixit1 varname [, diff(alternative1 alternative2)]
// where "varname" should be a constant unity variable. 

// Written by Won-ho Park Jan. 2022.


capture prog drop fixit1
program fixit1, eclass
      syntax varname(numeric)[, DIFF (namelist min=2 max=2)]

confirm variable `varlist'

qui sum `varlist'
local rmean = r(mean)
local rvar = r(Var)
if `rmean'~=1|`rvar'~=0{
      di "Specify the Dummy unity variable in the form <fixit1 unity>"
      error 9022
}

capture mat drop B V b V2 alts
      local xvars = e(xvars)
      mat B = e(b)
      mat V= e(V)
      tokenize `xvars'
      macro shift
      local vars `*'
*      di "`vars'"
local eqs
      foreach var in `vars'{
            local eqs "`eqs' `var': "
      }
*      di "`eqs'"
*      mat list B
matselrc B B2, c(`eqs')
*mat list B2
* Now switch coleqname and colname

local varlist: coleq B2
local eqnames: colnames B2
mat colnames B2 = `varlist'
*mat list B2

mat alts = e(altvals)
local numchoices=colsof(alts)
local alternatives: colnames alts
*di "`alternatives'"
*local xvars= r(xvars)
local newcollect
foreach x in `alternatives'{
      foreach y in `vars'{
            local newcollect "`newcollect' `x':`y'"
*            di "`newcollect'"
      }
}
*di "`newcollect'"

local numvar = colsof(B2)/`numchoices'
*di "`numvar'"
local collect
      forvalues j = 1/`numvar'{
            local collect "`collect' `alternatives'"
*            di "`j' `collect'"
      }

mat coleq B2 = `collect'
*mat list e(b)
*mat list B
*mat list B2
*di "`newcollect'"
matselrc B2 b, c(`newcollect')
* mat list b

* Now take care of the VCE matrix. 
*mat list V
* matselrc V v, c(`newcollect') r(`newcollect')
matselrc V V2, c(`eqs') r(`eqs')
* mat list V2
mat roweq V2 = `collect'
mat rownames V2 = `varlist'
mat coleq V2 = `collect'
mat colnames V2 = `varlist'
* mat list V2
matselrc V2 V3, c(`newcollect') r(`newcollect')
* mat list V3
mat V=V3
*mat list b
*mat list V
ereturn repost b=b V=V, resize
mat err = e(error)[1, 1.. colsof(B2)]
mat obs = e(_N)[1, 1.. colsof(B2)]
ereturn mat error = err
ereturn mat _N = obs


if "`diff'" ~= ""{
      tokenize `diff'
      local alt1 `1'
      local alt2 `2'
      di "You chose to compute the difference between `alt1' and `alt2'. "

      matrix bMNP = J(1,`numvar',0)
      matrix sMNP = J(1,`numvar',0)
      scalar j = 0
      foreach var in `vars' {
                  scalar j = j + 1
                  *di "`var'"
                  qui lincom ([`alt1']`var' - [`alt2']`var')/2
                  matrix bMNP[1,j] = r(estimate)
                  matrix sMNP[1,j] = r(se)
      *            mat list sMNP
            }
      mat coleq bMNP = "Diff"
      mat colnames bMNP = `vars'
      mat coleq sMNP = "Diff"
      mat colnames sMNP = `vars'
*      mat list bMNP
*      mat list sMNP      
      mat V = diag(sMNP)*diag(sMNP)
      mat err = e(error)[1, 1..`numvar']
      mat obs = e(_N)[1, 1..`numvar']
      ereturn mat error = err
      ereturn mat _N = obs
      ereturn repost b=bMNP V=V, resize
}

ereturn display

end
