% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.infinityV.R
\name{check.infinityV}
\alias{check.infinityV}
\title{Check Infinity in Vector}
\usage{
check.infinityV(V, algebraType)
}
\arguments{
\item{V}{is vector.}

\item{algebraType}{is string input that can be \code{minplus} or \code{maxplus}.}
}
\value{
Returns nothing but generates an error if specific conditions met.
}
\description{
\code{check.infinityV} checks infinite value in a vector based on \code{algebraType} input.
}
\details{
The input of this function is a vector and type of tropical algebra. A vector may contain infinite
values that can be positive or negative. Both the positive and negative infinite values works differently on
each algebra type. Due to the difference between \code{minplus} and \code{maxplus} tropical algebra, it is
important to manage them so they can work in their own bounderies. In \code{minplus} -Inf cannot be used while in
\code{maxplus} Inf cannot be used. So the main purpose of this funnction is to check such possibilities that can
cause errors. If this function finds a -Inf in the vector and the type of algebra is \code{minplus} then the
function generates an error. Similarly, if the function finds a Inf in the vector and the type of algebra is
\code{maxplus} then the function also generates an error.
}
\examples{
a <- c(2, Inf, Inf, 0, Inf, Inf, Inf, 10, Inf)
check.infinityV(a, 'minplus')

}
\seealso{
\code{\link{check.infinityM}}
}
