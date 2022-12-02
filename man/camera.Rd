\name{camera}
\alias{camera}
\docType{data}

\title{Conjoint Survey Data for Digital Cameras}

\description{
Panel dataset from a conjoint survey for digital cameras with 332 respondents. Data exclude respondents that always answered none, always picked the same brand, always selected the highest priced offering, or who appeared to be answering randomly.
}

\usage{data(camera)}

\format{
A list of lists. Each inner list corresponds to one survey respondent and contains a numeric vector (\code{y}) of choice indicators and a numeric matrix (\code{X}) of covariates. Each respondent participated in 16 choice scenarios each including 4 camera options (and an outside option) for a total of 80 rows per respondent.
}

\details{
The covariates included in each \code{X} matrix are:
  \tabular{ll}{
    \ldots\code{$canon     } \tab an indicator for brand Canon \cr
    \ldots\code{$sony      } \tab an indicator for brand Sony \cr
    \ldots\code{$nikon     } \tab an indicator for brand Nikon \cr
    \ldots\code{$panasonic } \tab an indicator for brand Panasonic \cr
    \ldots\code{$pixels    } \tab an indicator for a higher pixel count \cr
    \ldots\code{$zoom      } \tab an indicator for a higher level of zoom \cr
    \ldots\code{$video     } \tab an indicator for the ability to capture video \cr
    \ldots\code{$swivel    } \tab an indicator for a swivel video display \cr
    \ldots\code{$wifi      } \tab an indicator for wifi capability \cr
    \ldots\code{$price     } \tab in hundreds of U.S. dollars 
  }
}

\source{
Allenby, Greg, Jeff Brazell, John Howell, and Peter Rossi (2014), "Economic Valuation of Product Features," \emph{Quantitative Marketing and Economics} 12, 421--456.

Allenby, Greg, Jeff Brazell, John Howell, and Peter Rossi (2014), "Valuation of Patented Product Features," \emph{Journal of Law and Economics} 57, 629--663.
}

\references{
For analysis of a similar dataset, see Case Study 4, \emph{Bayesian Statistics and Marketing} Rossi, Allenby, and McCulloch.
}

\keyword{datasets}
