\name{mnpProb}
\alias{mnpProb}
\concept{MNP}
\concept{Multinomial Probit Model}
\concept{GHK}
\concept{market share simulator}

\title{Compute MNP Probabilities}

\description{
\code{mnpProb} computes MNP probabilities for a given \eqn{X} matrix corresponding to one observation. This function can be used with output from \code{rmnpGibbs} to simulate the posterior distribution of market shares or fitted probabilties.
}

\usage{mnpProb(beta, Sigma, X, r)}

\arguments{
  \item{beta}{ MNP coefficients }
  \item{Sigma}{ Covariance matrix of latents }
  \item{X}{ \eqn{X} array for one observation -- use \code{createX} to make }
  \item{r}{ number of draws used in GHK (def: 100)}
}

\details{
See \code{\link{rmnpGibbs}} for definition of the model and the interpretation of the beta and Sigma parameters. Uses the GHK method to compute choice probabilities. To simulate a distribution of probabilities, loop over the beta and Sigma draws from \code{rmnpGibbs} output.
}

\value{\eqn{p x 1} vector of choice probabilites}

\author{ Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{ For further discussion, see Chapters 2 and 4, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\seealso{ \code{\link{rmnpGibbs}}, \code{\link{createX}} }

\examples{
## example of computing MNP probabilites
## here Xa has the prices of each of the 3 alternatives

Xa    = matrix(c(1,.5,1.5), nrow=1)
X     = createX(p=3, na=1, nd=NULL, Xa=Xa, Xd=NULL, DIFF=TRUE)
beta  = c(1,-1,-2)  ## beta contains two intercepts and the price coefficient
Sigma = matrix(c(1, 0.5, 0.5, 1), ncol=2)

mnpProb(beta, Sigma, X)
}

\keyword{models}
