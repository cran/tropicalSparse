#' @title tropicalsparse.axpyi()
#'
#' @description \code{tropicalsparse.axpyi} function multiplies the vector \code{x} by the constant \code{alpha}
#' and adds the result to the vector \code{y}.
#'
#' @param y is a vector.
#' @param alpha is a constant.
#' @param x is a vector.
#' @param algebraType is string input that can be \code{minplus} or \code{maxplus}.
#'
#' @details The input of the function \code{tropicalsparse.axpyi} is two vectors, a constant and type of Tropical
#' Algebra. \code{algebraType} is used to specify type of Tropical Algebra. This can be \code{minplus} or
#' \code{maxplus}. For more details about \code{algebraType}, see \code{detail} section of
#' \code{\link{check.infinityM}} or \code{\link{check.infinityV}}. All inputs of this method are compulsory. The
#' operation is expressed as: y = y + alpha * x where \code{x} and \code{y} must be a vector while \code{alpha}
#' must be a constant.
#'
#' @return Returns a vector.
#'
#' @seealso \code{\link{tropicalsparse.doti}}
#'
#' @examples
#' a <- c(2, Inf, 5, 0, Inf, Inf, Inf, 10, Inf)
#' b <- c(0, 5, Inf, Inf, 12, 2, Inf, Inf, 3)
#' alpha <- 5
#'
#' tropicalsparse.axpyi(a, alpha, b, 'minplus')
#'
#'# [1]   2  10   5   0  17   7 Inf  10   8
#'
#' @export
#'
tropicalsparse.axpyi <- function(y, alpha, x, algebraType){
  if(!(is.vector(x)) || !(is.vector(y))){
    stop("'x' or 'y' is not a vector")
  }

  if(is.array(alpha) || is.character(alpha) || is.complex(alpha) || is.expression(alpha) || is.factor(alpha) || is.function(alpha) || is.list(alpha) || is.matrix(alpha) || is.symbol(alpha) || is.table(alpha) || !(length(alpha) == 1) || is.null(alpha)){
    stop("'alpha' must be a 'single Real' value")
  }

  if(length(x) != length(y)){
    stop("non-comfortable arguments")
  }
  return(tropicalsparse.add(y, tropicalsparse.mul(x, alpha, algebraType = algebraType), algebraType = algebraType))
}
