#' @title Count non-infinit values
#'
#' @description \code{counter} is used to get total number of non-infinit values in a matrix.
#'
#' @param M is a matrix.
#'
#' @details The input of this function is a matrix. This function returns the total number of non-infinite values
#' in the given matrix. In order to work properly, \code{M} must be a \code{\link{matrix}} otherwise this method
#' generates an error.
#'
#' @return Returns total number of non-infinit values.
#'
#' @examples
#' a <- matrix(data = c(2, Inf, Inf, 0, Inf, Inf, Inf, 10, Inf), nrow = 3, ncol = 3, byrow = TRUE)
#'
#' counter(a)
#'
#'# [1] 3
#'
#' @export
#'
counter = function(M){
  i = 1
  x = 0
  repeat{
    if(i>nrow(M))
      break
    j = 1
    repeat{
      if(j>ncol(M))
        break
      if(!is.infinite(M[i,j]))
      {
        x <- x+1
      }
      j = j+1
    }
    i = i+1
  }
  return(x)
}
