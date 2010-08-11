\name{mnpProb}
\alias{mnpProb}
\concept{MNP}
\concept{Multinomial Probit Model}
\concept{GHK}
\concept{market share simulator}
\title{ Compute MNP Probabilities }
\description{
  \code{mnpProb} computes MNP probabilities for a given X matrix corresponding to one 
   observation.  This function can be used with output from \code{rmnpGibbs} to simulate
   the posterior distribution of market shares or fitted probabilties.
}
\usage{
mnpProb(beta, Sigma, X, r)
}
\arguments{
  \item{beta}{ MNP coefficients }
  \item{Sigma}{ Covariance matrix of latents }
  \item{X}{ X array for one observation -- use \code{createX} to make }
  \item{r}{ number of draws used in GHK (def: 100)}
}
\details{
  see \code{\link{rmnpGibbs}} for definition of the model and the interpretation of
  the beta, Sigma parameters. Uses the GHK method to compute choice probabilities.
  To simulate a distribution of probabilities, loop over the beta, Sigma draws from
  \code{rmnpGibbs} output.
}
\value{
 p x 1 vector of choice probabilites 
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi,Allenby and McCulloch,  Chapters 2 and 4. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\seealso{ \code{\link{rmnpGibbs}}, \code{\link{createX}} }
\examples{
##
## example of computing MNP probabilites
##  here I'm thinking of Xa as having the prices of each of the 3 alternatives
Xa=matrix(c(1,.5,1.5),nrow=1)
X=createX(p=3,na=1,nd=NULL,Xa=Xa,Xd=NULL,DIFF=TRUE)
beta=c(1,-1,-2)  ## beta contains two intercepts and the price coefficient
Sigma=matrix(c(1,.5,.5,1),ncol=2)
mnpProb(beta,Sigma,X)
}
\keyword{ models }
