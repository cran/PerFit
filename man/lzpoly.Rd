\name{lzpoly}
\alias{lzpoly}
\title{lzpoly person-fit statistic}
\description{Compute the lzpoly (Drasgow, Levine, and Williams, 1985) person-fit statistic.}
\usage{lzpoly(matrix,Ncat,ip = NA, model = "GRM", ability = NA, method = "EAP")}
\arguments{
  \item{matrix}{A data matrix of polytomous item scores: Persons as rows, items as columns, item scores are integers between 0 and (Ncat-1), no missing data.
  }
  \item{Ncat}{Number of answer options for each item.}
  \item{ip}{Matrix with previously estimated item parameters: One row per item. The first (Ncat-1) columns contain the between-categories threshold parameters (for the GRM) or the item step difficulties (for the PCM and the GPCM). The last, Ncat-th, column has the slopes.
  
  If \code{ip} is provided then \code{ability} must also be provided.
  
  In case no item parameters are available then \code{ip=NA}.}
  \item{model}{Specify the IRT model to use in order to estimate the item parameters. The options available are \code{"PCM"}, \code{"GPCM"}, and \code{"GRM"} (default).}
  \item{ability}{Vector with previoulsy estimated latent ability parameters, one per respondent, following the order of the row index of \code{matrix}.
  
  In case no ability parameters are available then \code{ability=NA}.}
  \item{method}{Specify the method to use in order to estimate the latent ability parameters. The options available are \code{"EB"}, \code{"EAP"} (default), and \code{"MI"}.}
}
\details{
Statistic \code{lzpoly} is the natural extension of \code{\link{lz}} to polytomously scores items. In this case the user can choose one from three possible IRT models to fit the data: The partial credit model (\code{model="PCM"}), the generalized partial credit model (\code{model="GPCM"}), or the graded response model (\code{model="GRM"}). Ability parameters can be estimated by means of one of three methods: Empirical Bayes (\code{"EB"}), expected a posteriori (\code{"EAP"}), or multiple imputation (\code{"MI"}).

Both item and ability parameters may be provided as function parameters (\code{ip} and \code{ability}, respectively). If \code{ip} is provided then \code{ability} must also be provided. The reason is that the estimation of the ability parameters is done via the function \code{factor.scores} from the \pkg{ltm} package, and this procedure requires an object of a specific class (\code{gpcm} or \code{grm}) containing the estimated item parameters (i.e., providing a matrix of item parameters to \code{ip} is not sufficient).

Aberrant response behavior is (potentially) indicated by small values of lzpoly (i.e., in the left tail of the sampling distribution).
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
\seealso{\code{\link{lz}},\code{\link{lzstar}}}
\examples{
# Load the physical functioning data (polytomous item scores):
data(PhysFuncData);

# Compute the lzpoly scores:
lzpoly(PhysFuncData,Ncat=3);
}
\keyword{univar}
