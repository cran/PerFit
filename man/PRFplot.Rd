\name{PRFplot}
\alias{PRFplot}
\alias{PRFplot}
\title{Person response function (PRF)}
\description{Plot the nonparametric person response function with variability bands.}
\usage{
PRFplot(matrix, respID, h = 0.09, N.FPts = 101, VarBands = FALSE,
 VarBands.area = FALSE, alpha = 0.05,
 Xlabel=NA,Xcex=1.5,Ylabel=NA,Ycex=1.5,title=NA,Tcex=1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{A data matrix of dichotomous item scores: Persons as rows, items as columns, item scores are either 0 or 1, no missing data.}
  \item{respID}{Row index identifiying the respondent for whom the PRF is to be computed.
}
  \item{h}{Bandwidth value. Default is 0.09.}
  \item{N.FPts}{Number of (equidistant) focal points in the [0,1] interval. Default is 101.}
  \item{VarBands}{Logical: Draw the \code{(1-alpha)} variability bands? Default is \code{FALSE}.}
  \item{VarBands.area}{Logical: Draw the area between the \code{(1-alpha)} variability bands? Default is \code{FALSE}.}
  \item{alpha}{Significance level for the variability bands. Default is 0.05.}
  \item{Xlabel}{Define label of x-axis, otherwise a default label is shown.}
  \item{Xcex}{Font size of the label of x-axis.}
  \item{Ylabel}{Define label of y-axis, otherwise a default label is shown.}
  \item{Ycex}{Font size of the label of y-axis.}
  \item{title}{Define the title of the plot, otherwise a default title is shown.}
  \item{Tcex}{Font size of the title of the plot.}
}
\details{This function plots the so-called nonparametric person response function (PRF; Emons, Sijtsma, and Meijer, 2004; Sijtsma and Meijer, 2001) for respondent \code{respID}. The PRF relates item difficulty (0-1 range on the x-axis) with the associated probability of correct response (on the y-axis). The PRF is typically nonincreasing, implying that the probability of answering increasingly difficult items should (at least) not increase.

The code is based on nonparametric kernel smoothing (Emons et al., 2004). The value of the PRF at each focal point (representing a difficulty parameter between 0 and 1) is estimated as a weighted sum score, where scores pertaining to items with difficulty close to the focal point are given the largest weights. 

The weights are functions of the Gaussian kernel function. It is necessary to specify a bandwidth value (\code{h}) in order to compute the weights. The \code{h} value controls the trade-off between bias and sampling variation (Emons et al., 2004). Small \code{h} values reduce bias but increase variance, leading to PRFs that capture too much measurement error. Large \code{h} values, on the other hand, increase bias which renders PRFs that are often too flat, thus missing potentially relevant misfitting response behavior. Therefore, it is important to carefuly specify the value \code{h}. Emons et al. (2004, pp. 10-13), after a simulation study, advised that "\code{h} values between 0.07 and 0.11 are reasonable choices".

Variability bands of level \code{alpha} (0.05 by default) can also be added to the plot. These bands are computed following the jackknife procedure explained in Emons et al. (2004).
}
\references{
Emons, W. M., Sijtsma, K., and Meijer, R. R. (2004) Testing hypotheses about the person-response function in person-fit analysis. \emph{Multivariate Behavioral Research}, \bold{39(1)}, 1--35.

Sijtsma, K., and Meijer, R. R. (2001) The person response function as a tool in person-fit research. \emph{Psychometrika}, \bold{66(2)}, 191--207.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{\code{\link{cutoff}}, \code{\link{plot.PerFit}}, \code{\link{flagged.resp}}}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# As an example, compute the Ht person-fit scores:
Ht.PF <- Ht(InadequacyData);
Ht.PF$PFscores

# Determine which respondents were flagged by Ht at 1%:
flagged.resp(InadequacyData,x=Ht.PF,Qlvl=.01);
# Flagged respondents: 30, 37, 46, 49,...

# Plot the person response function of respondent 30 (flagged as aberrant):
PRFplot(InadequacyData,respID=30);
# Plot the person response function of respondent 35 (not flagged as aberrant):
PRFplot(InadequacyData,respID=35);
}
\keyword{univar}

