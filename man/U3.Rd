\name{U3, ZU3}
\alias{U3}
\alias{ZU3}
\title{U3, ZU3 person-fit statistics}
\description{van der Flier's U3 and ZU3 person-fit statistics.}
\usage{
U3(matrix)

ZU3(matrix)
}
\arguments{
  \item{matrix}{A data matrix of dichotomous item scores: Persons as rows, items as columns, item scores are either 0 or 1, no missing data.}
}
\details{
Suppose the items are ordered in decreasing proportion-correct score, \eqn{p_1>p_2>\ldots>p_I}{p1>p2>...>pI} (\eqn{I} = number of items). Given response vector \eqn{(X_1,X_2,\ldots,X_I)}{(X1,X2,...,XI)} with total score \eqn{S=\sum_iX_i}{S=sum(Xi; i=1,...,I)}, van der Flier (1980, 1982) defined the U3 statistic as
\deqn{\frac{\sum_{i=1}^S p_i-\sum_{i=1}^I p_i}{\sum_{i=1}^S p_i-\sum_{i=I-S+1}^I p_i}.}{
[sum(pi; i=1,...,S) - sum(pi; i=1,...,I)] / [sum(pi; i=1,...,S) - sum(pi; i=I-S+1,...,I)].}

U3 varies from 0 for perfect Guttman response vectors (i.e., with only the easiest items correct) through 1 for reversed Guttman response vectors (i.e., with only the hardest items correct). Hence, increasingly large U3 values provide stronger indications of answering misfit. 

U3 scores are dependent on the number-correct score, hence van der Flier proposed ZU3 as a standardization (formulas to compute E(U3) and Var(U3) can be found in van der Flier, 1982). ZU3 is supposed to be asymptotically approximated by the standard normal distribution, but this approximation is not without problems (see Emons, Meijer, and Sijtsma, 2002).

Response vectors (i.e., rows of \code{matrix}) that only consist of 0s or 1s need to be removed before the function can be used.
}
\value{
An object of class "PerFit", which is a list with elements:
  \item{PFscores }{A vector of length \eqn{N} (number of respondents) with the values of the person-fit statistic.}
  \item{PFstatistic}{The name of the person-fit statistic used.}
}
\references{
Emons, W. M., Meijer, R. R., and Sijtsma, K. (2002). Comparing simulated and theoretical sampling distributions of the U3 person-fit statistic. \emph{Applied Psychological Measurement}, \bold{26(1)}, 88--108.

Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.

Meijer, R. R., and Sijtsma, K. (2001) Methodology review: Evaluating person fit. \emph{Applied Psychological Measurement}, \bold{25(2)}, 107--135.

van der Flier, H. (1980) \emph{Vergelijkbaarheid van individuele testprestaties} [Comparability of individual test performance]. Lisse: The Netherlands.

van der Flier, H. (1982) Deviant response patterns and comparability of test scores. \emph{Journal of Cross-Cultural Psychology}, \bold{13(3)}, 267--298.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{\code{\link{U3poly}}}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# Compute the U3 scores:
U3(InadequacyData);

# Compute the ZU3 scores:
ZU3(InadequacyData);
}
\keyword{univar}
