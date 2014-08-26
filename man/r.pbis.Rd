\name{Personal biserial statistic}
\alias{r.pbis}
\title{Personal biserial statistic}
\description{Donlon and Fischer's personal biserial statistic.}
\usage{r.pbis(matrix)}
\arguments{
  \item{matrix}{A data matrix of dichotomous item scores: Persons as rows, items as columns, item scores are either 0 or 1, no missing data.}
}
\details{
Donlon and Fischer (1968) suggested to use the correlation between a respondent's score vector and the item proportion-correct scores in the sample as a measure of person fit. Low values should be indicative of misfit of the response vector with respect to the group of respondents.

Response vectors (i.e., rows of \code{matrix}) that only consist of 0s or 1s need to be removed before the function can be used.
}
\value{
An object of class "PerFit", which is a list with elements:
  \item{PFscores }{A vector of length \eqn{N} (number of respondents) with the values of the person-fit statistic.}
  \item{PFstatistic}{The name of the person-fit statistic used.}
}
\references{
Donlon, T. F., and Fischer, F. E. (1968) An index of an individual's agreement with group-defined item difficulties. \emph{Educational and Psychological Measurement}, \bold{28(1)}, 105--113.

Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.

Meijer, R. R., and Sijtsma, K. (2001) Methodology review: Evaluating person fit. \emph{Applied Psychological Measurement}, \bold{25(2)}, 107--135.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# Compute the r.pbis scores:
r.pbis(InadequacyData);
}
\keyword{univar}
