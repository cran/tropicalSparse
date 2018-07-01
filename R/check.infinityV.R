#' @title Check Infinity in Vector
#'
#' @description \code{check.infinityV} checks infinite value in a vector based on \code{algebraType} input.
#'
#' @param V is vector.
#' @param algebraType is string input that can be \code{minplus} or \code{maxplus}.
#'
#' @details The input of this function is a vector and type of tropical algebra. A vector may contain infinite
#' values that can be positive or negative. Both the positive and negative infinite values works differently on
#' each algebra type. Due to the difference between \code{minplus} and \code{maxplus} tropical algebra, it is
#' important to manage them so they can work in their own bounderies. In \code{minplus} -Inf cannot be used while in
#' \code{maxplus} Inf cannot be used. So the main purpose of this funnction is to check such possibilities that can
#' cause errors. If this function finds a -Inf in the vector and the type of algebra is \code{minplus} then the
#' function generates an error. Similarly, if the function finds a Inf in the vector and the type of algebra is
#' \code{maxplus} then the function also generates an error.
#'
#' @return Returns nothing but generates an error if specific conditions met.
#'
#' @seealso
#' \code{\link{check.infinityM}}
#'
#' @examples
#' a <- c(2, Inf, Inf, 0, Inf, Inf, Inf, 10, Inf)
#' check.infinityV(a, 'minplus')
#'
#' @export
#'
check.infinityV <- function(V, algebraType) {

  i = 1

  if(algebraType == "minplus"){
    while(i<=length(V)){
        if(V[i] == -Inf){
          stop("Invalid input\n-Inf at Vector[", i, "]")
        }
        i = i + 1
      }

  }else if(algebraType == "maxplus"){
    while(i<=length(V)){
      if(V[i] == Inf){
        stop("Invalid input\nInf at Vector[", i, "]")
      }
      i = i + 1
    }
  }
}
