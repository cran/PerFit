\name{A, D, E statistics}
\alias{ADE.KB}
\alias{A.KB}
\alias{D.KB}
\alias{E.KB}
\title{Agreement, disagreement, and dependability statistics}
\description{Kane and Brennan's person-fit statistics.}
\usage{
A.KB(matrix)

D.KB(matrix)

E.KB(matrix)
}
\arguments{
  \item{matrix}{A data matrix of dichotomous item scores: Persons as rows, items as columns, item scores are either 0 or 1, no missing data.}
}
\details{
Kane and Brennan (1980) discussed the agreement, disagreement, and dependability statistics. Assume that the items are ordered in increasing difficulty order (i.e., according to decreasing proportion-correct score). The agreement statistic for respondent \eqn{n} is
\deqn{A.KB=\sum_i{X_i}p_i,}{A.KB = sum(Xi*pi),}
where \eqn{X_i}{Xi} is the 0-1 score of respondent \eqn{n} on item \eqn{i} and \eqn{p_i}{pi} is the proportion-correct score of item \eqn{i}.

The disagreement statistic is
\deqn{D.KB = A.KB(max) - A.KB,}
where \eqn{A.KB(max)} is the maximum value of A.KB given respondent \eqn{n}'s total score.

The dependability statistic is
\deqn{E.KB = A.KB / A.KB(max).}

Small values of A.Kb and E.KB (i.e., in the left tail of the sampling distribution) are (potentially) indicative of aberrant response behavior. Large values of D.Kb (i.e., in the right tail of the sampling distribution) are (potentially) indicative of aberrant response behavior.

Response vectors (i.e., rows of \code{matrix}) that only consist of 0s or 1s need to be removed before the function can be used.
}
\value{
An object of class "PerFit", which is a list with elements:
  \item{PFscores }{A vector of length \eqn{N} (number of respondents) with the values of the person-fit statistic.}
  \item{PFstatistic}{The name of the person-fit statistic used.}
}
\references{
Kane, M. T., and Brennan, R. L. (1980) Agreement coefficients as indices of dependability for domain-referenced tests. \emph{Applied Psychological Measurement}, \bold{4(1)}, 105--126.

Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.

Meijer, R. R., and Sijtsma, K. (2001) Methodology review: Evaluating person fit. \emph{Applied Psychological Measurement}, \bold{25(2)}, 107--135.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# Compute the A.KB, D.KB, and E.KB scores:
A.KB(InadequacyData);
D.KB(InadequacyData);
E.KB(InadequacyData);
}
\keyword{univar}
