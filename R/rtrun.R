rtrun=
function(mu,sigma,a,b){
#
# function to draw from univariate truncated norm
# a is vector of lower bounds for truncation
# b is vector of upper bounds for truncation
#
FA=pnorm(((a-mu)/sigma))
FB=pnorm(((b-mu)/sigma))
return(mu+sigma*qnorm(runif(length(mu))*(FB-FA)+FA))
}
