rdirichlet = 
function(alpha)
{
#
# Purpose:
# draw from Dirichlet(alpha)
#
dim = length(alpha)
y=rep(0,dim)
for(i in 1:dim) y[i] = rgamma(1,alpha[i])
return(y/sum(y))
}
