#' @title Storage Techniques
#'
#' @description \code{tropicalsparse.storage} function is used to apply \code{coo}, \code{csr} and \code{csc}
#' storage techniques on the sparse matrix in Tropical Algebra.
#'
#' @param M is Matrix
#' @param store is storage technique.
#' @param algebraType is string input that can be \code{minplus} or \code{maxplus}.
#'
#' @details The function \code{tropicalsparse.storage} recieves a matrix as first input, storage technique as
#' second input and the type of Tropical Algebra as third input. All the inputs are compulsory. \code{store} can
#' be \code{coo}, \code{csr} and \code{csc}. \code{algebraType} is used to specify type of Tropical Algebra. This
#' can be \code{minplus} or \code{maxplus}. For more details about \code{algebraType}, see \code{detail} section of
#' \code{\link{check.infinityM}} or \code{\link{check.infinityV}}. If store is equal to \code{coo} then the
#' function returns a list containing three arrays that are \code{row_Indices_COO}, \code{col_Indices_COO} and
#' \code{values_COO}. If store is equal to \code{csc} then the function returns a list containing three arrays
#' that are \code{col_Pointer_CSC}, \code{row_Indices_CSC} and \code{values_CSC}. If store is equal to \code{csr}
#' then the function returns a list containing three arrays that are \code{row_Pointer_CSR},
#' \code{col_Indices_CSR} and \code{values_CSR}. These storage techniques are especially designed for sparse
#' matrices and are very helpful and time saving.
#'
#' @return Returns a list \code{result} that contains three arrays depends upon the \code{store} input.
#'
#' @seealso
#' \code{\link{tropicalsparse.add}}, \code{\link{tropicalsparse.mv}}.
#'
#' @examples
#' a <- matrix(data = c(2, Inf, Inf, 0, Inf, Inf, Inf, 10, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' tropicalsparse.storage(a, 'coo', 'minplus')
#'
#'# $row_Indices_COO
#'# [1] 1 2 3
#'
#'# $col_Indices_COO
#'# [1] 1 1 2
#'
#'# $values_COO
#'# [1]  2  0 10
#'
#' @export
#'
tropicalsparse.storage <- function(M , store, algebraType){

  if(!is.matrix(M)){
    stop("Matrix required")
  }

  if(algebraType == 'minplus'){

    if(store == "csr"){
      row_Pointer_CSR = array(data = NA)
      col_Indices_CSR = array(data = NA)
      values_CSR = array(data = NA)
      x = 1
      i = 1
      k = 0
      row_Pointer_CSR[1] = k
      while(i<=nrow(M)){
        j = 1
        while(j<=ncol(M)){

          if(M[i,j] == -Inf){
            stop("Invalid input\n-Inf at Matrix[", i, ",", j, "]")
          }

          if(!is.infinite(M[i,j])){
            col_Indices_CSR[x] = j
            values_CSR[x] = M[i,j]
            x = x+1
            k = k+1
          }
          j = j+1
        }
        row_Pointer_CSR[i+1] = k
        i = i+1
      }
      k = length(M) - k
      sparcity = (k / length(M)) * 100
      if(sparcity < 60){
        warning("Matrix is not Sparse.")
        return(M)
      }

      result <- list()
      result[["row_Pointer_CSR"]] <- row_Pointer_CSR
      result[["col_Indices_CSR"]] <- col_Indices_CSR
      result[["values_CSR"]] <- values_CSR

      return(result)

    }else if(store == "csc"){
      col_Pointer_CSC = array(data = NA)
      row_Indices_CSC = array(data = NA)
      values_CSC = array(data = NA)
      x = 1
      i = 1
      k = 0
      col_Pointer_CSC[1] = k
      while(i<=ncol(M)){
        j = 1
        while(j<=nrow(M)){

          if(M[i,j] == -Inf){
            stop("Invalid input\n-Inf at Matrix[", i, ",", j, "]")
          }

          if(!is.infinite(M[j,i])){
            row_Indices_CSC[x] = j
            values_CSC[x] = M[j,i]
            x = x+1
            k = k+1
          }
          j = j+1
        }
        col_Pointer_CSC[i+1] = k
        i = i+1
      }
      k = length(M) - k
      sparcity = (k / length(M)) * 100
      if(sparcity < 60){
        warning("Matrix is not Sparse.")
        return(M)
      }

      result <- list()
      result[["col_Pointer_CSC"]] <- col_Pointer_CSC
      result[["row_Indices_CSC"]] <- row_Indices_CSC
      result[["values_CSC"]] <- values_CSC

      return(result)

    }else if(store == "coo"){
      row_Indices_COO = array(data = NA)
      col_Indices_COO = array(data = NA)
      values_COO = array(data = NA)
      i = 1
      x = 1
      k = 0
      repeat{
        if(i>nrow(M))
          break
        j = 1
        repeat{
          if(j>ncol(M))
            break

          if(M[i,j] == -Inf){
            stop("Invalid input\n-Inf at Matrix[", i, ",", j, "]")
          }

          if(!is.infinite(M[i,j])){

            row_Indices_COO[x] <- i
            col_Indices_COO[x] <- j
            values_COO[x] <- M[i, j]
            x <- x + 1
            k = k +1
          }
          j = j+1
        }
        i = i+1
      }
      k = length(M) - k
      sparcity = (k / length(M)) * 100
      if(sparcity < 60){
        warning("Matrix is not Sparse.")
        return(M)
      }

      result <- list()
      result[["row_Indices_COO"]] <- row_Indices_COO
      result[["col_Indices_COO"]] <- col_Indices_COO
      result[["values_COO"]] <- values_COO

      return(result)

    }else{
      stop('"store" argument is invalid')
    }
  }else if(algebraType == 'maxplus'){

    if(store == "csr"){
      row_Pointer_CSR = array(data = NA)
      col_Indices_CSR = array(data = NA)
      values_CSR = array(data = NA)
      x = 1
      i = 1
      k = 0
      row_Pointer_CSR[1] = k
      while(i<=nrow(M)){
        j = 1
        while(j<=ncol(M)){

          if(M[i,j] == Inf){
            stop("Invalid input\nInf at Matrix[", i, ",", j, "]")
          }

          if(!is.infinite(M[i,j])){
            col_Indices_CSR[x] = j
            values_CSR[x] = M[i,j]
            x = x+1
            k = k+1
          }
          j = j+1
        }
        row_Pointer_CSR[i+1] = k
        i = i+1
      }
      k = length(M) - k
      sparcity = (k / length(M)) * 100
      if(sparcity < 60){
        warning("Matrix is not Sparse.")
        return(M)
      }

      result <- list()
      result[["row_Pointer_CSR"]] <- row_Pointer_CSR
      result[["col_Indices_CSR"]] <- col_Indices_CSR
      result[["values_CSR"]] <- values_CSR

      return(result)

    }else if(store == "csc"){
      col_Pointer_CSC = array(data = NA)
      row_Indices_CSC = array(data = NA)
      values_CSC = array(data = NA)
      x = 1
      i = 1
      k = 0
      col_Pointer_CSC[1] = k
      while(i<=ncol(M)){
        j = 1
        while(j<=nrow(M)){

          if(M[i,j] == Inf){
            stop("Invalid input\nInf at Matrix[", i, ",", j, "]")
          }

          if(!is.infinite(M[j,i])){
            row_Indices_CSC[x] = j
            values_CSC[x] = M[j,i]
            x = x+1
            k = k+1
          }
          j = j+1
        }
        col_Pointer_CSC[i+1] = k
        i = i+1
      }
      k = length(M) - k
      sparcity = (k / length(M)) * 100
      if(sparcity < 60){
        warning("Matrix is not Sparse.")
        return(M)
      }

      result <- list()
      result[["col_Pointer_CSC"]] <- col_Pointer_CSC
      result[["row_Indices_CSC"]] <- row_Indices_CSC
      result[["values_CSC"]] <- values_CSC

      return(result)

    }else if(store == "coo"){
      row_Indices_COO = array(data = NA)
      col_Indices_COO = array(data = NA)
      values_COO = array(data = NA)
      i = 1
      x = 1
      k = 0
      repeat{
        if(i>nrow(M))
          break
        j = 1
        repeat{
          if(j>ncol(M))
            break

          if(M[i,j] == Inf){
            stop("Invalid input\nInf at Matrix[", i, ",", j, "]")
          }

          if(!is.infinite(M[i,j])){

            row_Indices_COO[x] <- i
            col_Indices_COO[x] <- j
            values_COO[x] <- M[i, j]
            x <- x+1
            k = k + 1
          }
          j = j+1
        }
        i = i+1
      }
      k = length(M) - k
      sparcity = (k / length(M)) * 100
      if(sparcity < 60){
        warning("Matrix is not Sparse.")
        return(M)
      }

      result <- list()
      result[["row_Indices_COO"]] <- row_Indices_COO
      result[["col_Indices_COO"]] <- col_Indices_COO
      result[["values_COO"]] <- values_COO

      return(result)

    }else{
      stop('"store" argument is invalid')
    }
  }else{
    stop("Invalid algebra type")
  }
}
