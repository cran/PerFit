\name{flagged.resp}
\alias{flagged.resp}
\title{Find (potentially) aberrant response patterns}
\description{Find which respondents in the sample were flagged by the specified person-fit statistic.}
\usage{
flagged.resp(matrix, scores = TRUE, ord = TRUE, x, method = "Quantile", 
 Qlvl = 0.05, Blvl = 0.05, Breps = 1000, UDlvl = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{A data matrix of item scores (persons as rows, items as columns). This is the same data matrix used to compute the \code{PerFit} object \code{x}.}
  \item{scores}{Logical: Should item scores of flagged respondents be shown in the output? Default is \code{TRUE}.}
  \item{ord}{Logical: Should items be ordered in increasing order of difficulty (i.e., in decreasing proportion-correct order)? Default is \code{TRUE}. Only used if \code{scores=TRUE}.}
  \item{x}{Object of class "PerFit".}
  \item{method}{Method required to compute the cutoff score: \code{"Quantile"}, \code{"Bootstrap"}, \code{"UserDefined"}. Default is \code{"Quantile"}.}
  \item{Qlvl}{Quantile level (value between 0 and 1), only used when \code{method="Quantile"}. Default is 0.05.}
  \item{Blvl}{Bootstrap level (value between 0 and 1), only used when \code{method="Bootstrap"}. Default is 0.05.}
  \item{Breps}{Number of bootstrap resamples, only used when \code{method="Bootstrap"}. Default is 1000.}
  \item{UDlvl}{User-defined level, only used when \code{method="UserDefined"}.}
}
\details{This function finds the respondents in the dataset that were flagged by the person-fit statistic. This statistic is specified by means of the "PerFit" class object \code{x}. 

The cutoff score is internally computed, for which the user may define \code{method}, \code{Qlvl}, \code{Blvl}, \code{Breps}, and \code{UDlvl} accordingly (see \code{\link{cutoff}}).

If \code{scores=TRUE} then the respondents' item scores will be shown in the output, either in the original item order (\code{ord=FALSE}) or in increasing difficulty order (\code{ord=TRUE}).
}
\value{
If \code{scores=FALSE} the output is a matrix with two columns:
  \item{\code{FlaggedID} }{Row index identifying the flagged respondents.}
  \item{\code{PFscores} }{The values of the person-fit statistic.}
If \code{scores=TRUE} the output is a list with two elements:
  \item{\code{Scores} }{Matrix with columns: \code{FlaggedID}, item scores (\code{It**}), and \code{PFscores}.}
  \item{\code{MeanItemValue} }{The items' mean value (which is nothing more than the proportion-correct for dichotomous items).}
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{
\code{\link{cutoff}}, \code{\link{plot.PerFit}}, \code{\link{PRFplot}}
}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# As an example, compute the Ht person-fit scores:
Ht.PF <- Ht(InadequacyData);
Ht.PF$PFscores

# Determine which respondents were flagged by Ht for the quantile-based 1% cutoff:
flagged.resp(InadequacyData,x=Ht.PF,Qlvl=.01);
}
\keyword{univar}

