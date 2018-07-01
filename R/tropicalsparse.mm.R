#' @title tropicalsparse.mm()
#'
#' @description \code{tropicalsparse.mm} function performs the matrix-matrix operation on the equation: y = alpha * op(A) * op(B) + beta * op(C).
#'
#' @param alpha is a single real value.
#' @param A is a matrix.
#' @param opA is transpose of \code{A}.
#' @param B is a matrix.
#' @param opB is transpose of \code{B}.
#' @param beta is a single real value.
#' @param C is a matrix.
#' @param opC is transpose of \code{C}.
#' @param store is storage technique.
#' @param algebraType is string input that can be \code{minplus} or \code{maxplus}.
#'
#' @details The input of this function is three matrices, two constants, operation on these matrices, storage
#' technique and  type of Tropical Algebra. Matrices and \code{algebraType} inputs are compulsory while all other
#' inputs are optional. \code{A}, \code{B} and \code{C} must be the matrix of same dimensions and these matrices
#' must be sparse otherwise error occured. \code{alpha} and \code{beta} must be the single real value. \code{opA},
#' \code{opB} and \code{opC} can be set to TRUE to take transpose of \code{A}, \code{B} and \code{C} matrices
#' repectively. \code{store} input can be \code{coo}, \code{csc} and \code{csr} for applying following storage
#' techniques respectively: Coordinate-Wise, Compressed Sparse Row, Compressed Sparse Column. If \code{store} is
#' not specified then functionality is performed without using any storage technique. \code{algebraType} is used
#' to specify type of Tropical Algebra. This can be \code{minplus} or \code{maxplus}. For more details about
#' \code{algebraType}, see \code{detail} section of \code{\link{check.infinityM}} or \code{\link{check.infinityV}}.
#' First of all \code{A} is multiplied with \code{B} and if \code{alpha} is given then the product of \code{A} and
#' \code{B} will be multipied with \code{alpha} otherwise it remais the same. After this, \code{beta} is
#' multiplied with \code{C} only if \code{beta} is given. Finally, both result are added to each other and the
#' resultant matrix is obtained. Same functionailty is applied if any of the \code{store} technique is specified.
#'
#' @return Returns the resultant matrix.
#'
#' @examples
#' a <- matrix(data = c(2, Inf, Inf, 0, Inf, Inf, Inf, 10, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' b <- matrix(data = c(Inf, Inf, 4, Inf, -0.3, Inf, Inf, 2, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' c <- matrix(data = c(1, Inf, Inf, Inf, 0, 6, Inf, Inf, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' tropicalsparse.mm(A = a, alpha = 5, opB = TRUE, B = b, C = c,
#' store = 'csr', algebraType = 'minplus')
#'
#'#      [,1] [,2] [,3]
#'# [1,]    1  Inf  Inf
#'# [2,]  Inf  0.0    6
#'# [3,]  Inf 14.7   17
#'
#' @export
#'
tropicalsparse.mm <- function(alpha = NULL, A, opA = FALSE, B, opB = FALSE, beta = NULL, C, opC = FALSE, store = NULL, algebraType){

  if(!is.matrix(A) || !is.matrix(B) || !is.matrix(C)){
    stop("Invalid Input")
  }

  if(opA == TRUE){
    if(!is.null(A)){
      A = t(A)
    }else{
      stop("matrix 'A' is missing")
    }
  }

  if(opB == TRUE){
    if(!is.null(B)){
      B = t(B)
    }else{
      stop("matrix 'B' is missing")
    }
  }

  if(opC == TRUE){
    if(!is.null(C)){
      C = t(C)
    }else{
      stop("matrix 'C' is missing")
    }
  }

  if(is.array(alpha) || is.character(alpha) || is.complex(alpha) || is.expression(alpha) || is.factor(alpha) || is.function(alpha) || is.list(alpha) || is.matrix(alpha) || is.symbol(alpha) || is.table(alpha) || (!(length(alpha) == 1) && !(is.null(alpha)))){
    stop("'alpha' must be a 'single Real' value")
  }

  if(is.array(beta) || is.character(beta) || is.complex(beta) || is.expression(beta) || is.factor(beta) || is.function(beta) || is.list(beta) || is.matrix(beta) || is.symbol(beta) || is.table(beta) || (!(length(beta) == 1) && !(is.null(beta)))){
    stop("'beta' must be a 'single Real' value")
  }

  if(is.null(store)){

    D = tropicalsparse.mul(A, B, algebraType = algebraType)

    if(!is.null(alpha)){
      D = tropicalsparse.mul(D, alpha, algebraType = algebraType)
    }

    if(!is.null(beta)){
      C = tropicalsparse.mul(C, beta, algebraType = algebraType)
    }

    return(tropicalsparse.add(D, C, algebraType = algebraType))

  }else if(store == 'coo'){

    D = tropicalsparse.mul(A, B, 'coo', algebraType)

    if(!is.null(alpha)){
      D = tropicalsparse.mul(D, alpha, 'coo', algebraType)
    }

    if(!is.null(beta)){
      C = tropicalsparse.mul(C, beta, 'coo', algebraType)
    }

    return(tropicalsparse.add(D, C, 'coo', algebraType))

  }else if(store == 'csr'){

    D = tropicalsparse.mul(A, B, 'csr', algebraType)

    if(!is.null(alpha)){
      D = tropicalsparse.mul(D, alpha, 'csr', algebraType)
    }

    if(!is.null(beta)){
      C = tropicalsparse.mul(C, beta, 'csr', algebraType)
    }

    return(tropicalsparse.add(D, C, 'csr', algebraType))

  }else if(store == 'csc'){

    D = tropicalsparse.mul(A, B, 'csc', algebraType)

    if(!is.null(alpha)){
      D = tropicalsparse.mul(D, alpha, 'csc', algebraType)
    }

    if(!is.null(beta)){
      C = tropicalsparse.mul(C, beta, 'csc', algebraType)
    }

    return(tropicalsparse.add(D, C, 'csc', algebraType))

  }else{
    stop('"store" argument is invalid')
  }
}
