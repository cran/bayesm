cgetC = function(e,k) 
{
# purpose: get a list of cutoffs for use with scale usage problems
#
# arguments:
#   e: the "e" parameter from the paper
#   k: the point scale, eg. items are rated from 1,2,...k
# output:
#   vector of grid points
temp = (1:(k-1))+.5
m1 = sum(temp)
m2 = sum(temp^2)
return(.C('getC',as.double(e),as.integer(k),as.double(m1),as.double(m2),cc=double(k+1))$cc)
}
