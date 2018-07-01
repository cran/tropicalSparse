#' @title tropicalsparse.mv()
#'
#' @description \code{tropicalsparse.mv} function performs the matrix-vector operation on the equation: y = alpha * op(A) * x + beta * y.
#'
#' @param alpha is a single real value.
#' @param A is a matrix.
#' @param opA is transpose of \code{A}.
#' @param x is a vector.
#' @param beta is a single real value.
#' @param y is a vector.
#' @param store is storage technique.
#' @param algebraType is string input that can be \code{minplus} or \code{maxplus}.
#'
#' @details The input of this function is one matrix, transpose of that matrix, two vectors, two constants,
#' storage technique and  type of Tropical Algebra. The inputs of the matrix, vectors and \code{algebraType} are
#' compulsory while all other inputs are optional. The matrix must be sparse otherwise error occured. \code{alpha}
#' and \code{beta} must be a single real value. \code{opA} can be set to TRUE to take transpose of \code{A}.
#' \code{store} input can be \code{coo}, \code{csc} and \code{csr} for applying following storage techniques
#' respectively: Coordinate-Wise, Compressed Sparse Row, Compressed Sparse Column. If \code{store} is not
#' specified then functionality is performed without using any storage technique. \code{algebraType} is used to
#' specify type of Tropical Algebra. This can be \code{minplus} or \code{maxplus}. For more details about
#' \code{algebraType}, see \code{detail} section of \code{\link{check.infinityM}} or \code{\link{check.infinityV}}.
#' First of all \code{A} is multiplied with \code{x} and if \code{alpha} is given then the product of \code{A}
#' and \code{x} will be multipied with \code{alpha} otherwise it remais the same. After this, \code{beta} is
#' multiplied with \code{y} only if \code{beta} is given. Finally, both result are added to each other and the
#' resultant matrix is obtained. Same functionailty is applied if any of the \code{store} technique is specified.
#'
#' @return Returns the resultant matrix.
#'
#' @examples
#' a <- matrix(data = c(2, Inf, Inf, 0, Inf, Inf, Inf, 10, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' b <- c(2, Inf, 5)
#'
#' c <- c(Inf, 9, Inf)
#'
#' tropicalsparse.mv(A = a, alpha = 5, opA = TRUE, x = b, y = c,
#' store = 'csr', algebraType = 'minplus')
#'
#'#      [,1] [,2] [,3]
#'# [1,]    9  Inf  Inf
#'# [2,]    9    9    9
#'# [3,]  Inf   20  Inf
#'
#' @export
#'
tropicalsparse.mv <- function(alpha = NULL, A, opA = FALSE, x, beta = NULL, y, store = NULL, algebraType){

  if(!is.matrix(A) || !is.vector(x) || !is.vector(y)){
    stop("Invalid Input")
  }

  if(opA == TRUE){
    if(!is.null(A)){
      A = t(A)
    }else{
      stop("matrix 'A' is missing")
    }
  }

  if(is.array(alpha) || is.character(alpha) || is.complex(alpha) || is.expression(alpha) || is.factor(alpha) || is.function(alpha) || is.list(alpha) || is.matrix(alpha) || is.symbol(alpha) || is.table(alpha) || (!(length(alpha) == 1) && !(is.null(alpha)))){
    stop("'alpha' must be a 'single Real' value")
  }

  if(is.array(beta) || is.character(beta) || is.complex(beta) || is.expression(beta) || is.factor(beta) || is.function(beta) || is.list(beta) || is.matrix(beta) || is.symbol(beta) || is.table(beta) || (!(length(beta) == 1) && !(is.null(beta)))){
    stop("'beta' must be a 'single Real' value")
  }

  if(is.null(store)){

    D = tropicalsparse.mul(A, x, algebraType = algebraType)

    if(!is.null(alpha)){
      D = tropicalsparse.mul(D, alpha, algebraType = algebraType)
    }

    if(!is.null(beta)){
      y = tropicalsparse.mul(y, beta, algebraType = algebraType)
    }

    return(tropicalsparse.add(D, y, algebraType = algebraType))

  }else if(store == 'coo'){

    D = tropicalsparse.mul(A, x, 'coo', algebraType)

    if(!is.null(alpha)){
      D = tropicalsparse.mul(D, alpha, 'coo', algebraType)
    }

    if(!is.null(beta)){
      y = tropicalsparse.mul(y, beta, 'coo', algebraType)
    }

    return(tropicalsparse.add(D, y, 'coo', algebraType))

  }else if(store == 'csr'){

    D = tropicalsparse.mul(A, x, 'csr', algebraType)

    if(!is.null(alpha)){
      D = tropicalsparse.mul(D, alpha, 'csr', algebraType)
    }

    if(!is.null(beta)){
      y = tropicalsparse.mul(y, beta, 'csr', algebraType)
    }

    return(tropicalsparse.add(D, y, 'csr', algebraType))

  }else if(store == 'csc'){

    D = tropicalsparse.mul(A, x, 'csc', algebraType)

    if(!is.null(alpha)){
      D = tropicalsparse.mul(D, alpha, 'csc', algebraType)
    }

    if(!is.null(beta)){
      y = tropicalsparse.mul(y, beta, 'csc', algebraType)
    }

    return(tropicalsparse.add(D, y, 'csc', algebraType))

  }else{
    stop('"store" argument is invalid')
  }
}
