\name{cutoff}
\alias{cutoff}
\title{Compute a cutoff value given the scores of a person-fit statistic}
\description{Compute a cutoff value given the scores of a person-fit statistic.}
\usage{cutoff(x, method = "Quantile", Qlvl = 0.05, Blvl = 0.05, Breps = 1000, UDlvl = NA)}
\arguments{
  \item{x}{Object of class "PerFit".}
  \item{method}{Method required to compute the cutoff score: \code{"Quantile"}, \code{"Bootstrap"}, \code{"UserDefined"}. Default is \code{"Quantile"}.}
  \item{Qlvl}{Quantile level (value between 0 and 1), only used when \code{method="Quantile"}. Default is 0.05.}
  \item{Blvl}{Bootstrap level (value between 0 and 1), only used when \code{method="Bootstrap"}. Default is 0.05.}
  \item{Breps}{Number of bootstrap resamples, only used when \code{method="Bootstrap"}. Default is 1000.}
  \item{UDlvl}{User-defined level, only used when \code{method="UserDefined"}.}
}
\details{
This function computes a reference value (referred to as a 'cutoff') associated to the values of a person-fit statistic computed from a sample. The idea is to create a decision rule: Individual person-fit values at or more extreme than the cutoff result in flagging the corresponding respondents as (potentially) displaying aberrant response behavior.

Depending on the person-fit statistic, an "extreme" score might be a very small (e.g., for \code{Ht}) or a very large (e.g., for \code{G}) value. The \code{cutoff} function routinely reports of which type the person-fit statistic being used is (tail="lower" or tail="upper", respectively).

There are three methods available to estimate the cutoff value. When \code{method="Quantile"} the cutoff is the \code{Qlvl} (resp. \code{1-Qlvl}) quantile of the sampling distribution for "lower" (resp. "upper") types of person-fit statistics. When \code{method="Bootstrap"} the cutoff is the median of the bootstrap distribution estimated by computing the \code{Blvl} (resp. \code{1-Blvl}) quantile from each bootstrap  resample (in a total of \code{Breps}) for "lower" (resp. "upper") types of person-fit statistics. Finally, the cutoff can be manually entered by the user (e.g., when it is available from prior data calibration) when \code{method="UserDefined"}.
}
\value{
A list with elements:
  \item{cutoff }{Numeric. Value of the computed cutoff.}
  \item{prop.flagged }{Numeric. Proportion of respondents flagged (that is, with person-fit scores at or more extreme than the cutoff). It is equal to \code{Qlvl} when \code{method="Quantile"} and approximately equal to \code{Blvl} when \code{method="Bootstrap"}.}
  \item{tail}{String with values "lower" or "upper". It indicates the type of person-fit statistic.}
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{
\code{\link{flagged.resp}}, \code{\link{plot.PerFit}}, \code{\link{PRFplot}}
}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# As an example, compute the Ht person-fit scores:
Ht.PF <- Ht(InadequacyData);

# Compute the quantile-based 1% cutoff:
cutoff(Ht.PF,Qlvl=.01);
}
\keyword{univar}

