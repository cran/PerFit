\name{U3poly}
\alias{U3poly}
\title{U3poly person-fit statistic}
\description{Generalization of van der Flier's U3 person-fit statistic to polytomously scored items.}
\usage{U3poly(matrix,Ncat)}
\arguments{
  \item{matrix}{A data matrix of polytomous item scores: Persons as rows, items as columns, item scores are integers between 0 and (Ncat-1), no missing data.}
  \item{Ncat}{Number of answer options for each item.}
}
\details{
Emons (2008) generalized the U3 statistic (van der Flier, 1980, 1982) to polytomous items. The idea is based on the so-called item-step difficulty, which is the probability of moving from answer category \eqn{(c)} to answer category \eqn{(c+1)} (\eqn{c=0,\ldots,Ncat-2}{c=0,...,Ncat-2}).

U3poly varies from 0 (no misfit) through 1 (extreme misfit). Hence, increasingly large U3poly values provide stronger indications of answering misfit. 

The number of answer options, \code{Ncat}, is the same for all items. 

\code{U3poly} reduces to \code{U3} when \code{Ncat=2}.
}
\value{
An object of class "PerFit", which is a list with elements:
  \item{PFscores }{A vector of length \eqn{N} (number of respondents) with the values of the person-fit statistic.}
  \item{PFstatistic}{The name of the person-fit statistic used.}
}
\references{
Emons, W. M. (2008) Nonparametric person-fit analysis of polytomous item scores. \emph{Applied Psychological Measurement}, \bold{32(3)}, 224--247. 

Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.

Meijer, R. R., and Sijtsma, K. (2001) Methodology review: Evaluating person fit. \emph{Applied Psychological Measurement}, \bold{25(2)}, 107--135.

van der Flier, H. (1980) \emph{Vergelijkbaarheid van individuele testprestaties} [Comparability of individual test performance]. Lisse: The Netherlands.

van der Flier, H. (1982) Deviant response patterns and comparability of test scores. \emph{Journal of Cross-Cultural Psychology}, \bold{13(3)}, 267--298.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{\code{\link{U3}}, \code{\link{ZU3}}}
\examples{
# Load the physical functioning data (polytomous item scores):
data(PhysFuncData);

# Compute the U3poly scores:
U3poly(PhysFuncData,Ncat=3);
}
\keyword{univar}
