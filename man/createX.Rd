\name{createX}
\alias{createX}

\title{ Create X Matrix for Use in Multinomial Logit and Probit Routines }
\description{
  \code{createX} makes up an X matrix in the form expected by Multinomial
  Logit (\code{\link{rmnlIndepMetrop}}) and Multinomial Probit (\code{\link{rmnpGibbs}})
  routines.  Requires an array of alternative specfic variables and/or an
  array of "demographics" or variables constant across alternatives which
  may vary across choice occasions. 
}
\usage{
createX(p, na, nd, Xa, Xd, INT = TRUE, DIFF = FALSE, base = p)
}

\arguments{
  \item{p}{ integer - number of choice alternatives }
  \item{na}{ integer - number of alternative-specific vars in Xa  }
  \item{nd}{ integer - number of "demographics" vars }
  \item{Xa}{ n x p*na matrix of alternative-specific vars }
  \item{Xd}{ n x nd matrix of "demographic" vars }
  \item{INT}{ logical flag for inclusion of intercepts }
  \item{DIFF}{ logical flag for differencing wrt to base alternative }
  \item{base}{ integer - index of base choice alternative }
  note: na,nd,Xa,Xd can be NULL to indicate lack of Xa or Xd variables.
}

\value{
  X matrix -- n*(p-DIFF) x [(INT+nd)*(p-1) + na] matrix.
}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\seealso{ \code{\link{rmnlIndepMetrop}}, \code{\link{rmnpGibbs}} }
\examples{
na=2; nd=1; p=3
vec=c(1,1.5,.5,2,3,1,3,4.5,1.5)
Xa=matrix(vec,byrow=TRUE,ncol=3)
Xa=cbind(Xa,-Xa)
Xd=matrix(c(-1,-2,-3),ncol=1)
createX(p=p,na=na,nd=nd,Xa=Xa,Xd=Xd)
createX(p=p,na=na,nd=nd,Xa=Xa,Xd=Xd,base=1)
createX(p=p,na=na,nd=nd,Xa=Xa,Xd=Xd,DIFF=TRUE)
createX(p=p,na=na,nd=nd,Xa=Xa,Xd=Xd,DIFF=TRUE,base=2)
createX(p=p,na=na,nd=NULL,Xa=Xa,Xd=NULL)
createX(p=p,na=NULL,nd=nd,Xa=NULL,Xd=Xd)
}
\keyword{ array }
\keyword{ utilities }
