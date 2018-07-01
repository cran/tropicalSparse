#' @title Row/Column Number of a Value
#'
#' @description \code{row.col.Number} method is used to get the row or column number of a specific value in the
#' matrix.
#'
#' @param i is an index of array containing non-infinite values of the matrix.
#' @param x is total number of rows or columns of the matrix.
#' @param arr is an array containing row or column pointer of the matrix.
#'
#' @details The function \code{row.col.Number} recieves three parameters \code{i}, \code{x} and \code{arr}. As
#' mentioned above \code{i} is an index of array containing non-infinite values of the matrix. This array can only
#' be obtained in the \code{CSR} and \code{CSC} storage techniques and has zero sparsity. \code{x} is total number
#' of rows in case of \code{CSR} or total number of columns in case of \code{CSC} of the matrix. \code{arr} is an
#' array containing row pointer in case of \code{CSR} or column pointer in case of \code{CSC} of the matrix. From
#' these inputs \code{row.col.Number} finds row or column number of a specific value in the matrix. This function
#' is used especially for \code{CSR} and \code{CSC} storage techniques.
#'
#' @return Returns the row or column number of a specific value if succeded, otherwise \code{NA}.
#'
#' @seealso
#' \code{\link{tropicalsparse.add}}, \code{\link{tropicalsparse.mul}}.
#'
#' @examples
#' a <- matrix(data = c(2, Inf, Inf, 0, Inf, Inf, Inf, 10, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' List = tropicalsparse.storage(a,'csr','minplus')
#' i = 2
#' row.col.Number(i, nrow(a), List[[1]])
#'# [1] 2
#' @export
#'
row.col.Number = function(i, x, arr){
  k = 1
  while(k<=x){
    if(i<=arr[k+1]){
      return(k)
    }
    k = k+1
  }
  return(NA)
}
