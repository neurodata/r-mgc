# mgc
R package for MGC code


Installation:
open MGC.Rproj in R studio; 
library(MGC)


Demo:
Once installed, type 
> run_demo() 
which takes < 10 seconds to run, which shows the actual MGC and p-value for testing on linear and quadratic data.


For any two corresponded data X and Y of size n*d1 and size n*d2,
type 
> result=MGCPermutationTest(X,Y) 
to do dependence testing. Refer to the document for more details on input and output.