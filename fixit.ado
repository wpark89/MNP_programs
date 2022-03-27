// This is a simple utility program that just 'trims' the dimensions and corrects the names of coefficients 
// and VCE matrix, after a cmmprobit estimation with a "fake" unity alternative-specific variable. 
// This program should only be used after such estimation. Also note that it erases all the 
// stored "ereturned" estimations and replace them with just a new set of b and V. 
// Written by Won-ho Park Jan. 2022. 


capture prog drop fixit
program fixit, eclass
capture mat drop B V b V2 alts
      mat B = e(b)
      mat V2= e(V)
      mat alts = e(altvals)
      local cut = colsof(alts)+1
      mat define B = B[1, `cut' ...]
      mat V2 = V2[`cut' ..., `cut' ...]

      local varlist: coleq B
      matrix coleq B = ""
      matrix colnames B = `varlist'

      matrix coleq V2 = ""
      matrix roweq V2 = ""
      matrix colnames V2=`varlist'
      mat rownames V2 = `varlist'
      *ereturn repost b = B V = V2, rename
      mat b = B
      mat V = V2
	  
	  mat list b
	  mat list V
      ereturn post b V
 end
