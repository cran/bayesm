ghkvec = 
function(L,trunpt,above,r){
#
# R interface to GHK code -- allows for a vector of truncation points
# revision history-
#    P. Rossi 4/05
#
   dim=length(above)
   n=length(trunpt)/dim
return(.C('ghk_vec',as.integer(n),as.double(L),as.double(trunpt),
   as.integer(above),as.integer(dim), as.integer(r),res=double(n))$res)
}
