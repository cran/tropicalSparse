#' @title tropicalsparse.doti()
#'
#' @description \code{tropicalsparse.doti} function multiplies the vector \code{y} with the vector \code{x}.
#'
#' @param y is a vector.
#' @param x is a vector.
#' @param algebraType is string input that can be \code{minplus} or \code{maxplus}.
#'
#' @details The input of this function is \code{x}, \code{y} and \code{algebraType}. If any of the input is missing
#' then the function generates an error. The operation is expressed as: result = yx where \code{x} and \code{y}
#' must be a vector. \code{algebraType} is used to specify type of Tropical Algebra. This can be \code{minplus} or
#' \code{maxplus}. For more details about \code{algebraType}, see \code{detail} section of
#' \code{\link{check.infinityM}} or \code{\link{check.infinityV}}.
#'
#' @return Returns a vector.
#'
#' @seealso \code{\link{tropicalsparse.axpyi}}
#'
#' @examples
#' a <- c(2, Inf, 5, 0, Inf, Inf, Inf, 10, Inf)
#' b <- c(0, 5, Inf, Inf, 12, 2, Inf, Inf, 3)
#'
#' tropicalsparse.doti(a, b, 'minplus')
#'
#'# [1]   2 Inf Inf Inf Inf Inf Inf Inf Inf
#'
#' @export

tropicalsparse.doti <- function(x, y, algebraType){

  if(!(is.vector(x)) || !(is.vector(y))){
    stop("Both the 'vectors' are required")
  }

  if(length(x) != length(y)){
    stop("non-conformable arrays")
  }
  return(tropicalsparse.mul(x,y, algebraType = algebraType))
}
