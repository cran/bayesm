\name{customerSat}
\alias{customerSat}
\docType{data}

\title{Customer Satisfaction Data}

\description{
Responses to a satisfaction survey for a Yellow Pages advertising product. All responses are on a 10 point scale from 1 to 10 (1 is "Poor" and 10 is "Excellent").
}

\usage{data(customerSat)}

\format{
A data frame with 1811 observations on the following 10 variables:
  \tabular{ll}{
    \ldots\code{$q1 } \tab Overall Satisfaction \cr
    \ldots\code{$q2 } \tab Setting Competitive Prices \cr
    \ldots\code{$q3 } \tab Holding Price Increase to a Minimum \cr
    \ldots\code{$q4 } \tab Appropriate Pricing given Volume \cr
    \ldots\code{$q5 } \tab Demonstrating Effectiveness of Purchase \cr
    \ldots\code{$q6 } \tab Reach a Large Number of Customers \cr
    \ldots\code{$q7 } \tab Reach of Advertising \cr
    \ldots\code{$q8 } \tab Long-term Exposure \cr
    \ldots\code{$q9 } \tab Distribution \cr
    \ldots\code{$q10} \tab Distribution to Right Geographic Areas 
  }
}

\source{Rossi, Peter, Zvi Gilula, and Greg Allenby (2001), "Overcoming Scale Usage Heterogeneity," \emph{Journal of the American Statistical Association} 96, 20--31.}

\references{Case Study 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
data(customerSat)
apply(as.matrix(customerSat),2,table)
## see also examples for 'rscaleUsage'
}

\keyword{datasets}
