\name{bank}
\alias{bank}
\docType{data}

\title{Bank Card Conjoint Data}

\description{
A panel dataset from a conjoint experiment in which two partial profiles of credit cards were presented to 946 respondents from a regional bank wanting to offer credit cards to customers outside of its normal operating region. Each respondent was presented with between 13 and 17 paired comparisons. The bank and attribute levels are disguised to protect the proprietary interests of the cooperating firm.
}

\usage{data(bank)}

\format{
  The \code{bank} object is a list containing two data frames. The first, \code{choiceAtt}, provides choice attributes for the partial credit card profiles. The second, \code{demo}, provides demographic information on the respondents.
}

\details{
In the \code{choiceAtt} data frame:
  \tabular{ll}{
    \ldots\code{$id           } \tab respondent id \cr
    \ldots\code{$choice       } \tab profile chosen \cr
    \ldots\code{$Med_FInt     } \tab medium fixed interest rate \cr
    \ldots\code{$Low_FInt     } \tab low fixed interest rate\cr
    \ldots\code{$Med_VInt     } \tab variable interest rate\cr
    \ldots\code{$Rewrd_2      } \tab reward level 2 \cr
    \ldots\code{$Rewrd_3      } \tab reward level 3 \cr
    \ldots\code{$Rewrd_4      } \tab reward level 4 \cr
    \ldots\code{$Med_Fee      } \tab medium annual fee level \cr
    \ldots\code{$Low_Fee      } \tab low annual fee level \cr
    \ldots\code{$Bank_B       } \tab bank offering the credit card \cr
    \ldots\code{$Out_State    } \tab location of the bank offering the credit card \cr
    \ldots\code{$Med_Rebate   } \tab medium rebate level \cr
    \ldots\code{$High_Rebate  } \tab high rebate level \cr
    \ldots\code{$High_CredLine} \tab high credit line level \cr
    \ldots\code{$Long_Grace   } \tab grace period 
  }

The profiles are coded as the difference in attribute levels. Thus, that a "-1" means the profile coded as a choice of "0" has the attribute. A value of 0 means that the attribute was not present in the comparison.

In the \code{demo} data frame:
  \tabular{ll}{
    \ldots\code{$id    } \tab respondent id \cr
    \ldots\code{$age   } \tab respondent age in years \cr
    \ldots\code{$income} \tab respondent income category \cr
    \ldots\code{$gender} \tab female=1 
  }
}

\source{Allenby, Gregg and James Ginter (1995), "Using Extremes to Design Products and Segment Markets," \emph{Journal of Marketing Research}, 392--403.}

\references{Appendix A, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
data(bank)
cat(" table of Binary Dep Var", fill=TRUE)
print(table(bank$choiceAtt[,2]))
cat(" table of Attribute Variables", fill=TRUE)
mat = apply(as.matrix(bank$choiceAtt[,3:16]), 2, table)
print(mat)
cat(" means of Demographic Variables", fill=TRUE)
mat=apply(as.matrix(bank$demo[,2:3]), 2, mean)
print(mat)

## example of processing for use with rhierBinLogit
if(0) {

  choiceAtt = bank$choiceAtt
  Z = bank$demo
  
  ## center demo data so that mean of random-effects
  ## distribution can be interpreted as the average respondent
  Z[,1] = rep(1,nrow(Z))
  Z[,2] = Z[,2] - mean(Z[,2])
  Z[,3] = Z[,3] - mean(Z[,3])
  Z[,4] = Z[,4] - mean(Z[,4])
  Z = as.matrix(Z)
  
  hh = levels(factor(choiceAtt$id))
  nhh = length(hh)
  lgtdata = NULL
  for (i in 1:nhh) {
    y = choiceAtt[choiceAtt[,1]==hh[i], 2]
    nobs = length(y)
    X = as.matrix(choiceAtt[choiceAtt[,1]==hh[i], c(3:16)])
    lgtdata[[i]] = list(y=y, X=X)
  }
  cat("Finished Reading data", fill=TRUE)
  
  Data = list(lgtdata=lgtdata, Z=Z)
  Mcmc = list(R=10000, sbeta=0.2, keep=20)
  
  set.seed(66)
  out = rhierBinLogit(Data=Data, Mcmc=Mcmc)
  
  begin = 5000/20
  summary(out$Deltadraw, burnin=begin)
  summary(out$Vbetadraw, burnin=begin)
  
  ## plotting examples
  if(0){
    
    ## plot grand means of random effects distribution (first row of Delta)
    index = 4*c(0:13)+1
    matplot(out$Deltadraw[,index], type="l", xlab="Iterations/20", ylab="",
            main="Average Respondent Part-Worths")
    
    ## plot hierarchical coefs
    plot(out$betadraw)
    
    ## plot log-likelihood
    plot(out$llike, type="l", xlab="Iterations/20", ylab="", 
         main="Log Likelihood")
  }
}
}

\keyword{datasets}
