\name{plot PerFit}
\alias{plot.PerFit}
\alias{plot}
\title{Plot method for objects of class "PerFit"}
\description{Plot method for objects of class "PerFit".}
\usage{\method{plot}{PerFit}(x, type = "Density", both.scale = TRUE, cutoff = TRUE, method = "Quantile", 
 Qlvl = 0.05, Blvl = 0.05, Breps = 1000, UDlvl = NA,
 Xlabel=NA,Xcex=1.5,title=NA,Tcex=1.5,...)}
\arguments{
  \item{x}{Object of class "PerFit".}
  \item{type}{Type of plot: \code{"Histogram"}, \code{"Density"} (default), \code{"Both"}.}
  \item{both.scale}{Logical: Should the y-axis be adjusted so that both the histogram and the density graphics are completely visible? Default is \code{TRUE}.}
  \item{cutoff}{Logical: Should the estimated statistic cutoff be added to the plot? Default is \code{TRUE}.}
  \item{method}{Method required to compute the cutoff score: \code{"Quantile"}, \code{"Bootstrap"}, \code{"UserDefined"}. Default is \code{"Quantile"}.}
  \item{Qlvl}{Quantile level (value between 0 and 1), only used when \code{method="Quantile"}. Default is 0.05.}
  \item{Blvl}{Bootstrap level (value between 0 and 1), only used when \code{method="Bootstrap"}. Default is 0.05.}
  \item{Breps}{Number of bootstrap resamples, only used when \code{method="Bootstrap"}. Default is 1000.}
  \item{UDlvl}{User-defined level, only used when \code{method="UserDefined"}.}
  \item{Xlabel}{Define label of x-axis, otherwise a default label is shown.}
  \item{Xcex}{Font size of the label of x-axis.}
  \item{title}{Define the title of the plot, otherwise a default title is shown.}
  \item{Tcex}{Font size of the title of the plot.}
  \item{...}{Extra graphical parameters to be passed to \code{plot}.}
}
\details{
This function plots the empirical distribution of the scores of the person-fit statistic specified by the "PerFit" class object \code{x}. A histogram, density, or a combination of both displays is possible. 

The cutoff score is internally computed, for which the user may define \code{method}, \code{Qlvl}, \code{Blvl}, \code{Breps}, and \code{UDlvl} accordingly (see \code{\link{cutoff}}). The value of the cutoff is superimposed to the plot when \code{cutoff=TRUE}. In this case, the adequate "flagging region" is colored, thus indicating the range of values for which the person-fit statistic flags respondents as potentially displaying aberrant behavior. The option \code{both.scale} was introduced to help to better tune the scale of the y-axis.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{
\code{\link{cutoff}}, \code{\link{flagged.resp}}, \code{\link{PRFplot}}
}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# Compute the ZU3 scores:
ZU3.PF <- ZU3(InadequacyData);

# Plot the sample distribution of the ZU3 scores with the quantile 1% cutoff value superimposed:
plot(ZU3.PF,type="Both",Qlvl=.01);

# Plot the sample distribution of the ZU3 scores with a bootstrapped 5% cutoff value superimposed:
plot(ZU3.PF,type="Both",method="Bootstrap");
}
\keyword{univar}

