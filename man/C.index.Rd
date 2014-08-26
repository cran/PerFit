\name{The caution statistic}
\alias{Cstar}
\alias{C.Sato}
\title{C.Sato, Cstar person-fit statistics}
\description{Computes the caution statistic C.Sato and the modified caution statistic Cstar.}
\usage{
C.Sato(matrix)

Cstar(matrix)
}
\arguments{
  \item{matrix}{A data matrix of dichotomous item scores: Persons as rows, items as columns, item scores are either 0 or 1, no missing data.}
}
\details{
The C.Sato statistic (also refered to as C in the literature) was proposed by Sato (1975):
\deqn{C.Sato = 1-\frac{cov(x_n,p)}{cov(x_n^*,p)},}{C.Sato = 1-cov(Xn,p)/cov(Xn*,p),}
where \eqn{x_n}{Xn} is the 0-1 response vector of respondent \eqn{n}, \eqn{p} is the vector of item proportions-correct, and \eqn{x_n^*}{Xn*} is the so-called Guttman vector containing correct answers for the easiest items (i.e., with the largest proportion-correct values) only. C.Sato is zero for Guttman vectors and its value tends to increase for response vectors that depart from the group's answering pattern, hence warning the researcher to be \emph{cautious} about interpreting such item scores. Therefore, (potentially) aberrant response behavior is indicated by large values of C.Sato (i.e., in the right tail of the sampling distribution).

Harnisch and Linn (1981) proposed a modified version of the caution statistic which bounds the caution statistic between 0 and 1 (also referred to as C* or MCI  in the literature):
\deqn{Cstar = \frac{cov(x_n^*,p)-cov(x_n,p)}{cov(x_n^*,p)-cov(x_n',p)},}{Cstar = [cov(Xn*,p)-cov(Xn,p)] / [cov(Xn*,p)-cov(Xn',p)],}
where \eqn{x_n'}{Xn'} is the reversed Guttman vector containing correct answers for the hardest items (i.e., with the smallest proportion-correct values) only. Cstar is sensitive to the so-called Guttman errors. A Guttman error is a pair of scores (0,1), where the 0-score pertains to the easiest item and the 1-score pertains to the hardest item. Cstar ranges between 0 (perfect Guttman vector) and 1 (reversed Guttman error), thus larger values indicate potential aberrant response behavior.
}
\value{
An object of class "PerFit", which is a list with elements:
  \item{PFscores }{A vector of length \eqn{N} (number of respondents) with the values of the person-fit statistic.}
  \item{PFstatistic}{The name of the person-fit statistic used.}
}
\references{
Harnisch, D. L., and Linn, R. L. (1981) Analysis of item response patterns: Questionable test data and dissimilar curriculum practices. \emph{Journal of Educational Measurement}, \bold{18(3)}, 133--146.

Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.

Meijer, R. R., and Sijtsma, K. (2001) Methodology review: Evaluating person fit. \emph{Applied Psychological Measurement}, \bold{25(2)}, 107--135.

Sato, T. (1975) \emph{The construction and interpretation of S-P tables.} Tokyo: Meiji Tosho.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{\code{\link{Ht}}}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# Compute the C.Sato scores:
C.Sato(InadequacyData);

# Compute the Cstar scores:
Cstar(InadequacyData);
}
\keyword{univar}
