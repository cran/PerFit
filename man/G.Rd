\name{G, Gnormed}
\alias{G}
\alias{Gnormed}
\title{Number of Guttman errors}
\description{van der Flier's statistics based on the number of Guttman errors.}
\usage{
G(matrix)

Gnormed(matrix)
}
\arguments{
  \item{matrix}{A data matrix of dichotomous item scores: Persons as rows, items as columns, item scores are either 0 or 1, no missing data.}
}
\details{
Consider the items' proportion-correct scores, \eqn{p_1,p_2,\ldots,p_I}{p1,p2,...,pI} (\eqn{I} = number of items). A Guttman error consists of an item score pair \eqn{(X_i,X_j)=(0,1)}{(Xi,Xj)=(0,1)} with \eqn{p_i>p_j}{pi>pj}. Hence, there is a Guttman error when an easier item is answered incorrectly and a more difficult item is answered correctly.

\code{G} counts the number of (0,1) pairs given that the items are ordered in decreasing proportion-correct scores order. However, \code{G} depends on the total number of items for a given number-correct score. In particular, for a number-correct \eqn{S}, \code{G} has maximum equal to \eqn{S(I-S)}. \code{Gnormed} was created to bound \code{G} between 0 and 1 by dividing it by its maximum (conditional on the number-correct score). Hence, (potentially) aberrant response behavior is indicated by large values of G/Gnormed (i.e., in the right tail of the sampling distribution).

\code{Gnormed} is perfectly linearly related to Tatsuoka and Tatsuoka's (1982, 1983) NCI statistic (\eqn{NCI = 1-2Gnormed}).
}
\value{
An object of class "PerFit", which is a list with elements:
  \item{PFscores }{A vector of length \eqn{N} (number of respondents) with the values of the person-fit statistic.}
  \item{PFstatistic}{The name of the person-fit statistic used.}
}
\references{
Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.

Meijer, R. R. (1994) The number of Guttman errors as a simple and powerful person-fit statistic. \emph{Applied Psychological Measurement}, \bold{18(4)}, 311--314.

Meijer, R. R., and Sijtsma, K. (2001) Methodology review: Evaluating person fit. \emph{Applied Psychological Measurement}, \bold{25(2)}, 107--135.

Tatsuoka, K. K., and Tatsuoka, M. M. (1982) Detection of aberrant response patterns and their effect on dimensionality. \emph{Journal of Educational Statistics}, \bold{7(3)}, 215--231.

Tatsuoka, K. K., and Tatsuoka, M. M. (1983) Spotting erroneous rules of operation by the individual consistency index. \emph{Journal of Educational Measurement}, \bold{20(3)}, 221--230.

van der Flier, H. (1977) Environmental factors and deviant response patterns. In Y. H. Poortinga (Ed.), \emph{Basic problems in cross-cultural psychology.} Amsterdam: The Netherlands.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{\code{\link{NCI}}}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# Compute the G scores:
G(InadequacyData);

# Compute the Gnormed scores:
Gnormed(InadequacyData);
}
\keyword{univar}
