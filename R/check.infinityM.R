#' @title Check Infinity in Matrix
#'
#' @description \code{check.infinityM} checks infinite value in a matrix based on \code{algebraType} input.
#'
#' @param M is matrix.
#' @param algebraType is string input that can be \code{minplus} or \code{maxplus}.
#'
#' @details The input of this function is a matrix and type of tropical algebra. A matrix may contain infinite
#' values that can be positive or negative. Both the positive and negative infinite values works differently on
#' each algebra type. Due to the difference between \code{minplus} and \code{maxplus} tropical algebra, it is
#' important to manage them so they can work in their own bounderies. In \code{minplus} -Inf cannot be used while in
#' \code{maxplus} Inf cannot be used. So the main purpose of this funnction is to check such possibilities that can
#' cause errors. If this function finds a -Inf in the matrix and the type of algebra is \code{minplus} then the
#' function generates an error. Similarly, if the function finds a Inf in the matrix and the type of algebra is
#' \code{maxplus} then the function also generates an error.
#'
#' @return Returns nothing but generates an error if specific conditions met.
#'
#' @seealso
#' \code{\link{check.infinityV}}
#'
#' @examples
#' a <- matrix(data = c(2, Inf, Inf, 0, Inf, Inf, Inf, 10, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#' check.infinityM(a, 'minplus')
#'
#' @export
#'
check.infinityM <- function(M, algebraType) {

  i = 1

  if(algebraType == "minplus"){
    while(i<=nrow(M)){
      j = 1
      while (j<=ncol(M)) {
        if(M[i, j] == -Inf){
          stop("Invalid input\n-Inf at Matrix[", i,",", j, "]")
        }
        j = j + 1
      }
      i = i + 1
    }

  }else if(algebraType == "maxplus"){
    while(i<=nrow(M)){
      j = 1
      while (j<=ncol(M)) {
        if(M[i, j] == Inf){
          stop("Invalid input\nInf at Matrix[", i,",", j, "]")
        }
        j = j + 1
      }
      i = i + 1
    }
  }
}
