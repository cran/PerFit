\name{NCI statistic}
\alias{NCI}
\title{NCI person-fit statistic}
\description{Tatsuoka and Tatsuoka's NCI statistic.}
\usage{NCI(matrix)}
\arguments{
  \item{matrix}{A data matrix of dichotomous item scores: Persons as rows, items as columns, item scores are either 0 or 1, no missing data.}
}
\details{
The NCI person-fit statistic was introduced by Tatsuoka and Tatsuoka (1982, 1983). It is perfectly linearly related to van der Flier's (1977) Gnormed statistic (\eqn{NCI = 1-2Gnormed}), see \code{\link{G}} for mathematical details. 

NCI equals 1 for perfect Guttman vectors (i.e., when only the easiest items are answered correctly, given the total score) and equals -1 for reversed Guttman vectors (i.e., when only the hardest items are answered correctly, given the total score). Hence, (potentially) aberrant response behavior is indicated by small values of NCI (i.e., in the left tail of the sampling distribution).
}
\value{
An object of class "PerFit", which is a list with elements:
  \item{PFscores }{A vector of length \eqn{N} (number of respondents) with the values of the person-fit statistic.}
  \item{PFstatistic}{The name of the person-fit statistic used.}
}
\references{
Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.

Meijer, R. R., and Sijtsma, K. (2001) Methodology review: Evaluating person fit. \emph{Applied Psychological Measurement}, \bold{25(2)}, 107--135.

Tatsuoka, K. K., and Tatsuoka, M. M. (1982) Detection of aberrant response patterns and their effect on dimensionality. \emph{Journal of Educational Statistics}, \bold{7(3)}, 215--231.

Tatsuoka, K. K., and Tatsuoka, M. M. (1983) Spotting erroneous rules of operation by the individual consistency index. \emph{Journal of Educational Measurement}, \bold{20(3)}, 221--230.

van der Flier, H. (1977) Environmental factors and deviant response patterns. In Y. H. Poortinga (Ed.), \emph{Basic problems in cross-cultural psychology.} Amsterdam: The Netherlands.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{\code{\link{G}}}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# Compute the NCI scores:
NCI(InadequacyData);
}
\keyword{univar}
