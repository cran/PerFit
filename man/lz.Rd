\name{lz, lzstar}
\alias{lz}
\alias{lzstar}
\title{lz and lzstar person-fit statistics}
\description{Compute the lz (Drasgow, Levine, and Williams, 1985) and the lzstar (Snijders, 2001) person-fit statistics.}
\usage{lz(matrix, ip = NA, model = "2PL", ability = NA, method = "ML", mu = 0, sigma = 1)

lzstar(matrix, ip = NA, model = "2PL", method = "ML", mu = 0, sigma = 1)
}
\arguments{
  \item{matrix}{A data matrix of dichotomous item scores: Persons as rows, items as columns, item scores are either 0 or 1, no missing data.
  }
  \item{ip}{Matrix with previously estimated item parameters: One row per item, and three columns ([,1] item discrimination; [,2] item difficulty; [,3] lower-asymptote, also referred to as pseudo-guessing parameter).
  
  In case no item parameters are available then \code{ip=NA}.}
  \item{model}{Specify the IRT model to use in order to estimate the item parameters (only if \code{ip=NA}). The options available are \code{"1PL"}, \code{"2PL"} (default), and \code{"3PL"}.}
  \item{ability}{Vector with previoulsy estimated latent ability parameters, one per respondent, following the order of the row index of \code{matrix}.
  
  In case no ability parameters are available then \code{ability=NA}.}
  \item{method}{Specify the method to use in order to estimate the latent ability parameters (only if \code{ability=NA} for \code{lz}). The options available are \code{"ML"} (default), \code{"BM"}, and \code{"WL"}.}
  \item{mu}{Mean of the apriori distribution. Only used when \code{method="BM"}. Default is 0.}
  \item{sigma}{Standard deviation of the apriori distribution. Only used when \code{method="BM"}. Default is 1.}
}
\details{
Drasgow et al. (1985) introduced one of the most used person-fit statistics, \code{lz}. This statistic is the standardized log-likelihood of the respondent's response vector. \code{lz} is (supposed to be) asymptotically standard normally distributed. 

The computation of \code{lz} requires that both item and ability parameters are available. Function \code{lz} allows to user to enter his/her own item and ability parameter estimates (variables \code{ip} and \code{ability}, respectively). Alternatively, \code{lz} relies on functions available through the \code{\link{irtoys}} package for estimating the parameters. Specifically, the user can choose one from three possible IRT models to fit the data: \code{model="1PL"}, \code{model="2PL"}, or \code{model="3PL"}. As for estimating the ability parameters there are three possible methods: \code{method="ML"} (maximum likelihood), \code{method="BM"} (Bayes modal), or \code{method="WL"} (weighted likelihood).

It was later observed by several researchers (e.g., Molenaar and Hoijtink, 1990) that the asymptotic approximation only holds when true ability values are used. This limitation was overcome by Snijders (2001), who further developed \code{lz} into the \code{lzstar} statistic. An accessible paper that thoroughly explains the basic principles behind \code{lzstar} is Magis, Raiche, and Beland (2012). It is important to realize that not all item and/or ability estimation procedures can be used when computing \code{lzstar}. In particular, the estimation of the ability parameters is constrained (see Snijders, 2001, Equation 5). The \code{lzstar} algorithm internally estimates the ability parameters accordingly for one of three possible methods: \code{method="ML"} (maximum likelihood), \code{method="BM"} (Bayes modal), or \code{method="WL"} (weighted likelihood), see Magis et al. (2012).

Aberrant response behavior is (potentially) indicated by small values of lz/lzstar (i.e., in the left tail of the sampling distribution).
}
\value{
An object of class "PerFit", which is a list with elements:
  \item{PFscores }{A vector of length \eqn{N} (number of respondents) with the values of the person-fit statistic.}
  \item{PFstatistic}{The name of the person-fit statistic used.}
}
\references{
Drasgow, F., Levine, M. V., and Williams, E. A. (1985) Appropriateness measurement with polychotomous item response models and standardized indices. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{38(1)}, 67--86.

Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.

Magis, D., Raiche, G., and Beland, S. (2012) A didactic presentation of Snijders's l[sub]z[/sub] index of person fit with emphasis on response model selection and ability estimation. \emph{Journal of Educational and Behavioral Statistics}, \bold{37(1)}, 57--81.

Meijer, R. R., and Sijtsma, K. (2001) Methodology review: Evaluating person fit. \emph{Applied Psychological Measurement}, \bold{25(2)}, 107--135.

Molenaar, I. W., and Hoijtink, H. (1990) The many null distributions of person fit indices. \emph{Psychometrika}, \bold{55(1)}, 75--106. 

Snijders, T. B. (2001) Asymptotic null distribution of person fit statistics with estimated person parameter. \emph{Psychometrika}, \bold{66(3)}, 331--342. 
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{\code{\link{lzpoly}}}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# Compute the lz scores using a subsample of the first 200 response vectors:
lz(InadequacyData[1:200,]);

# Compute the lzstar scores using a subsample of the first 200 response vectors:
lzstar(InadequacyData[1:200,]);
}
\keyword{univar}
