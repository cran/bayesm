\name{customerSat}
\alias{customerSat}
\docType{data}
\title{ Customer Satisfaction Data}
\description{
  Responses to a satisfaction survey for a Yellow Pages advertising product.
  All responses are on a 10 point scale from 1 to 10 (10 is "Excellent"
  and 1 is "Poor")
}
\usage{data(customerSat)}
\format{
  A data frame with 1811 observations on the following 10 variables.
  \describe{
    \item{\code{q1}}{Overall Satisfaction}
    \item{\code{q2}}{Setting Competitive Prices}
    \item{\code{q3}}{Holding Price Increase to a Minimum}
    \item{\code{q4}}{Appropriate Pricing given Volume}
    \item{\code{q5}}{Demonstrating Effectiveness of Purchase}
    \item{\code{q6}}{Reach a Large # of Customers}
    \item{\code{q7}}{Reach of Advertising}
    \item{\code{q8}}{Long-term Exposure}
    \item{\code{q9}}{Distribution}
    \item{\code{q10}}{Distribution to Right Geographic Areas}
  }
}
\source{
  Rossi et al (2001), "Overcoming Scale Usage Heterogeneity,"
  \emph{JASA} 96, 20-31.
}
\references{
  Case Study 3, \emph{Bayesian Statistics and Marketing} by Rossi et al.\cr 
  \url{http://www.perossi.org/home/bsm-1}
}
\examples{
data(customerSat)
apply(as.matrix(customerSat),2,table)
}
\keyword{datasets}
