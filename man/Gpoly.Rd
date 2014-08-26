\name{Gpoly, Gnormed.poly}
\alias{Gpoly}
\alias{Gnormed.poly}
\title{Number of Guttman errors for polytomous items}
\description{Molenaar (1991) and Emons (2008) statistics, based on the number of Guttman errors, for polytomous items.}
\usage{
Gpoly(matrix,Ncat)

Gnormed.poly(matrix,Ncat)
}
\arguments{
  \item{matrix}{A data matrix of polytomous item scores: Persons as rows, items as columns, item scores are integers between 0 and (Ncat-1), no missing data.}
  \item{Ncat}{Number of answer options for each item.}
}
\details{
Molenaar (1991) generalized the \code{G} person-fit statistic to polytomous items, \code{Gpoly}. The idea is based on the so-called item-step difficulty, which is the probability of moving from answer category \eqn{(c)} to answer category \eqn{(c+1)} (\eqn{c=0,\ldots,Ncat-2}{c=0,...,Ncat-2}).

Just like \code{G}, \code{Gpoly} depends on the test length. Emons (2008) developed \code{Gnormed.poly}, which is a normalized version of \code{Gpoly}.

Aberrant response behavior is (potentially) indicated by large values of Gpoly/Gnormed.poly (i.e., in the right tail of the sampling distribution).

The number of answer options, \code{Ncat}, is the same for all items.

\code{Gpoly} reduces to \code{G}, and \code{Gnormed.poly} reduces to \code{Gnormed}, when \code{Ncat=2}.
}
\value{
An object of class "PerFit", which is a list with elements:
  \item{PFscores }{A vector of length \eqn{N} (number of respondents) with the values of the person-fit statistic.}
  \item{PFstatistic}{The name of the person-fit statistic used.}
}
\references{
Emons, W. M. (2008) Nonparametric person-fit analysis of polytomous item scores. \emph{Applied Psychological Measurement}, \bold{32(3)}, 224--247. 

Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.

Meijer, R. R. (1994) The number of Guttman errors as a simple and powerful person-fit statistic. \emph{Applied Psychological Measurement}, \bold{18(4)}, 311--314.

Molenaar, I. W. (1991) A weighted Loevinger H-coefficient extending Mokken scaling to multicategory items. \emph{Kwantitatieve Methoden}, \bold{12(37)}, 97--117.
}
\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{\code{\link{G}}, \code{\link{Gnormed}}}
\examples{
# Load the physical functioning data (polytomous item scores):
data(PhysFuncData);

# Compute the Gpoly scores:
Gpoly(PhysFuncData,Ncat=3);

# Compute the Gnormed.poly scores:
Gnormed.poly(PhysFuncData,Ncat=3);
}
\keyword{univar}
