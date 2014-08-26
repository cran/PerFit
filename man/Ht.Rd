\name{Ht statistic}
\alias{Ht}
\title{Ht person-fit statistic}
\description{Sijtsma's Ht person-fit statistic.}
\usage{Ht(matrix)}
\arguments{
  \item{matrix}{A data matrix of dichotomous item scores: Persons as rows, items as columns, item scores are either 0 or 1, no missing data.}
}
\details{
Sijtsma (1986) adapted a statistic introduced by Mokken (1971) that originally allowed assessing the scalability of an item to the Guttman (1944, 1950) model. The same statistic was applied by Sijtsma to the \emph{transposed} data in order to detect respondents that would not comply with the Guttman model. Assume, without loss of generality, that the rows of the data matrix are ordered by increasing order of total score \eqn{S_n}{Sn} (\eqn{n=1,\ldots,N}{n = 1,..., N}). The statistic formula is given by the ratio
\deqn{Ht=\frac{cov(x_n,r_{(n)})}{cov_{max}(x_n,r_{(n)})},}{Ht=Cov(Xn,R(n)) / Cov.max(Xn,R(n)),}
where \eqn{r_{(n)}}{R(n)} is the vector of total item scores computed \emph{excluding} respondent \eqn{n} and the denominator is the maximum covariance given the marginal. Hence, Ht is actually similar to Sato's \link{C.Sato}. 

Ht is maximum 1 for respondent \eqn{n} when no respondent with a total score smaller/larger than \eqn{S_n}{Sn} can answer an item correctly/incorrectly that respondent \eqn{n} has answered incorrectly/correctly, respectively. Ht equals zero when the average covariance of the response pattern of respondent \eqn{n} with the other response patterns equals zero. Hence, (potentially) aberrant response behavior is indicated by small values of Ht (i.e., in the left tail of the sampling distribution).

The Ht statistic was shown to perform relatively well in several simulation studies (Karabatsos, 2003; Sijtsma, 1986; Sijtsma and Meijer, 1992, Tendeiro and Meijer, 2013).

Response vectors (i.e., rows of \code{matrix}) that only consist of 0s or 1s need to be removed before the function can be used.
}
\value{
An object of class "PerFit", which is a list with elements:
  \item{PFscores }{A vector of length \eqn{N} (number of respondents) with the values of the person-fit statistic.}
  \item{PFstatistic}{The name of the person-fit statistic used.}
}
\references{
Guttman, L. (1944) A basis for scaling qualitative data. \emph{American Sociological Review}, \bold{9}, 139-150.

Guttman, L. (1950) The basis for scalogram analysis. In S. A. Stouffer, L. Guttman, E. A. Suchman, P. F. Lazarsfeld, S. A. Star & J. A. Claussen (Eds.), \emph{Measurement and precision} (pp. 60-90). Princeton NJ: Princeton University Press.

Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.

Mokken, R. J. (1971) \emph{A theory and procedure of scale analysis.} Berlin, Germany: De Gruyter.

Sijtsma, K. (1986) A coefficient of deviance of response patterns. \emph{Kwantitatieve Methoden}, \bold{7}, 131--145.

Sijtsma, K., and Meijer, R. R. (1992) A method for investigating the intersection of item response functions in Mokken's nonparametric IRT model. \emph{Applied Psychological Measurement}, \bold{16(2)}, 149-157.

Tendeiro, J. N., and Meijer, R. R. (2013) Detection of Invalid Test Scores: The Usefulness of Simple Nonparametric Statistics. Submitted.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{\code{\link{C.Sato}}}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# Compute the Ht scores:
Ht(InadequacyData);
}
\keyword{univar}
