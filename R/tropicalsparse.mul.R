#' @title Multiplication With or Without Storage Techniques
#'
#' @description \code{tropicalsparse.mul} function multiplies the provided inputs in Tropical Algebra based on
#' type of Tropical Algebra.
#'
#' @param A is matrix or vector.
#' @param B is matrix or vector.
#' @param store is storage technique.
#' @param algebraType is string input that can be \code{minplus} or \code{maxplus}.
#'
#' @details The compulsory inputs of the function \code{tropicalsparse.mul} are \code{A}, \code{B} and
#' \code{algebraType} while the remaining input is optional that is \code{store}. The inputs \code{A} and
#' \code{B} can be matrix/matrix, matrix/vector, vector/matrix and vector/vector otherwise the function generates
#' an error. For \code{A} and \code{B}, the order of the input does not matter. It can be in any of the following
#' way: the first input of the function is matrix and second input is a vector. Similarly, vise versa.
#' \code{store} can be \code{coo}, \code{csc} and \code{csr} for applying following storage techniques
#' respectively: Coordinate-Wise, Compressed Sparse Row, Compressed Sparse Column. This input is case sensitive.
#' If the \code{store} input is other than the specified storage techniques then the function generates an error.
#' The input \code{algebraType} is used to specify type of Tropical Algebra. This can be \code{minplus} or
#' \code{maxplus}. For more details about \code{algebraType}, see \code{detail} section of
#' \code{\link{check.infinityM}} or \code{\link{check.infinityV}}. \code{\link{tropicalsparse.storage}}
#' function is used to apply storage technique depends upon the input given. If \code{store} input is not
#' specified then the functionality will be performed without using any storage technique.
#'
#' @return Multiplication of \code{A} and \code{B} in Tropical Algebra.
#'
#' @seealso
#' \code{\link{tropicalsparse.add}}, \code{\link{tropicalsparse.storage}}
#'
#' @examples
#' a <- matrix(data = c(2, Inf, Inf, 0, Inf, Inf, Inf, 10, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' b <- matrix(data = c(Inf, Inf, 4, Inf, -0.3, Inf, Inf, 2, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' tropicalsparse.mul(a, b, 'csr', 'minplus')
#'
#'#      [,1] [,2] [,3]
#'# [1,]  Inf  Inf    6
#'# [2,]  Inf  Inf    4
#'# [3,]  Inf  9.7  Inf
#'
#'# also
#'
#' a <- matrix(data = c(5, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 10, 2),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' b <- matrix(data = c(-Inf, -Inf, 3, -Inf, 0.5, -Inf, 1.1, -Inf, -Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' tropicalsparse.mul(a, b, 'coo', 'maxplus')
#'
#'#      [,1] [,2] [,3]
#'# [1,] -Inf -Inf    8
#'# [2,] -Inf -Inf -Inf
#'# [3,]  3.1 10.5 -Inf
#'
#'# also
#'
#' a <- matrix(data = c(2, Inf, Inf, 0, Inf, Inf, Inf, 2, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' b <- c(Inf, 0, 10)
#'
#' tropicalsparse.mul(a, b, algebraType = 'minplus')
#'
#'#      [,1] [,2] [,3]
#'# [1,]  Inf  Inf  Inf
#'# [2,]    0  Inf  Inf
#'# [3,]  Inf   12  Inf
#'
#' @export
#'

tropicalsparse.mul <- function(A, B, store = NULL, algebraType){
  if (algebraType == "minplus") {
    if(is.matrix(A) && is.matrix(B)){
      if(ncol(A) != nrow(B)){
        stop("non-conformable arguments")
      }

      if(is.null(store)){

        check.infinityM(A, algebraType)
        check.infinityM(B, algebraType)

        result = matrix(Inf, nrow = nrow(A), ncol = ncol(B), byrow = TRUE)
        i = 1
        while(i<=nrow(A)){
          j = 1
          while(j<=ncol(B)){
            k = 1
            while(k<=ncol(A)){
              result[i,j] = min(result[i,j], A[i,k] + B[k,j])
              k = k+1
            }
            j = j+1
          }
          i = i+1
        }
        return(result)

      } else if(store == "coo"){

        counterA <- counter(A)
        counterB <- counter(B)
        COOA <- tropicalsparse.storage(A, 'coo', algebraType)
        if(is.matrix(COOA)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        row_IndicesA_COO = COOA[[1]]
        col_IndicesA_COO = COOA[[2]]
        valuesA_COO = COOA[[3]]
        COOB = tropicalsparse.storage(B, 'coo', algebraType)
        if(is.matrix(COOB)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        row_IndicesB_COO = COOB[[1]]
        col_IndicesB_COO = COOB[[2]]
        valuesB_COO = COOB[[3]]
        result = matrix(Inf, nrow = nrow(A), ncol = ncol(B), byrow = TRUE)
        i = 1

        while(i<=counterA){
          j = 1
          while(j<=counterB)
          {
            if(col_IndicesA_COO[i]==row_IndicesB_COO[j]){
              result[row_IndicesA_COO[i],col_IndicesB_COO[j]] = min(result[row_IndicesA_COO[i],col_IndicesB_COO[j]] , (valuesA_COO[i] + valuesB_COO[j]))
            }
            if(col_IndicesA_COO[i]<row_IndicesB_COO[j]){
              break
            }
            j = j+1
          }
          i = i+1
        }
        return(result)

      }else if(store == "csr"){

        counterA <- counter(A)
        counterB <- counter(B)
        CSRA <- tropicalsparse.storage(A, 'csr', algebraType)
        if(is.matrix(CSRA)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        row_PointerA_CSR = CSRA[[1]]
        col_IndicesA_CSR = CSRA[[2]]
        valuesA_CSR = CSRA[[3]]
        CSRB <- tropicalsparse.storage(B, 'csr', algebraType)
        if(is.matrix(CSRB)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        row_PointerB_CSR = CSRB[[1]]
        col_IndicesB_CSR = CSRB[[2]]
        valuesB_CSR = CSRB[[3]]
        result = matrix(Inf, nrow = nrow(A), ncol = ncol(B), byrow = TRUE)
        i = 1

        while(i<=counterA){
          B_Row_Number = col_IndicesA_CSR[i]
          A_Row_Number = row.col.Number(i,nrow(A),row_PointerA_CSR)
          j = row_PointerB_CSR[B_Row_Number] + 1

          while(j <= row_PointerB_CSR[B_Row_Number + 1]){
            result[A_Row_Number, col_IndicesB_CSR[j]] = min(result[A_Row_Number,col_IndicesB_CSR[j]], valuesA_CSR[i] + valuesB_CSR[j])
            j = j+1
          }
          i = i+1
        }
        return(result)

      }else if(store == "csc"){

        counterA <- counter(A)
        counterB <- counter(B)
        CSCA <- tropicalsparse.storage(A, 'csc', algebraType)
        if(is.matrix(CSCA)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        col_PointerA_CSC = CSCA[[1]]
        row_IndicesA_CSC = CSCA[[2]]
        valuesA_CSC = CSCA[[3]]
        CSCB <- tropicalsparse.storage(B, 'csc', algebraType)
        if(is.matrix(CSCB)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        col_PointerB_CSC = CSCB[[1]]
        row_IndicesB_CSC = CSCB[[2]]
        valuesB_CSC = CSCB[[3]]
        result = matrix(Inf, nrow = nrow(A), ncol = ncol(B), byrow = TRUE)
        i = 1

        while(i<=counterA){
          A_Col_Number = row.col.Number(i,ncol(A),col_PointerA_CSC)
          j = 1
          while(j<=counterB){
            if(A_Col_Number==row_IndicesB_CSC[j]){
              B_Col_Number = row.col.Number(j,ncol(B),col_PointerB_CSC)
              result[row_IndicesA_CSC[i], B_Col_Number] = min(result[row_IndicesA_CSC[i], B_Col_Number], valuesA_CSC[i] + valuesB_CSC[j])
            }
            j = j+1
          }
          i = i+1
        }
        return(result)

      }else{
        stop("Invalid Storage Format")
      }

    }else if(is.matrix(A) && is.vector(B) || is.matrix(B) && is.vector(A)){
      if(is.matrix(B) && is.vector(A)){
        temp = A
        A = B
        B = temp
      }
      result = matrix(Inf, nrow(A), ncol(A))
      if(is.null(store)){

        check.infinityM(A, algebraType)
        check.infinityV(B, algebraType)

        if(length(B) == 1){
          if(B != Inf){
            i = 1
            while(i<=nrow(A)){
              j = 1
              while(j<=ncol(A)){
                result[i, j] = A[i, j] + B
                j = j+1
              }
              i = i+1
            }
          }
          return(result)

        }else if(nrow(A) == length(B)){
          i = 1
          while(i<=ncol(A)){
            j = 1
            while(j<=length(B)){
              if(A[j,i] != Inf && B[j] != Inf){
                result[j, i] = A[j, i] + B[j]
              }
              j = j+1
            }
            i = i+1
          }
          return(result)

        }else{
          stop("longer object length is not a multiple of shorter object length")
        }
      }else if(store == "coo"){
        if(length(B) == 1){
          if(B != Inf){
            COO <- tropicalsparse.storage(A, 'coo', algebraType)
            if(is.matrix(COO)){
              return(tropicalsparse.mul(A,B,NULL,algebraType))
            }
            row_Indices_COO = COO[[1]]
            col_Indices_COO = COO[[2]]
            values_COO = COO[[3]]
            x = 1

            while(x <= length(values_COO)){
              i = row_Indices_COO[x]
              j = col_Indices_COO[x]
              result[i,j] = B + values_COO[x]
              x = x+1
            }
          }
          return(result)
        }else if(nrow(A) == length(B)){
          COO <- tropicalsparse.storage(A, 'coo', algebraType)
          if(is.matrix(COO)){
            return(tropicalsparse.mul(A,B,NULL,algebraType))
          }
          row_Indices_COO = COO[[1]]
          col_Indices_COO = COO[[2]]
          values_COO = COO[[3]]
          x1 = 1
          x = 1

          while(x <= nrow(A)){
            y = 1
            while(y <= ncol(A)){
              i = row_Indices_COO[x1]
              j = col_Indices_COO[x1]
              if(!(is.na(i))){
                if(i == x && j == y){
                  if(x1 <= length(values_COO)){
                    result[i,j] = B[x] + values_COO[x1]
                    x1 = x1+1
                  }
                }
              }
              y = y+1
            }
            x = x+1
          }
          return(result)

        }else{
          stop("longer object length is not a multiple of shorter object length")
        }
      }else if(store == "csr"){
        if(length(B) == 1){
          if(B != Inf){
            CSR <- tropicalsparse.storage(A, 'csr', algebraType)
            if(is.matrix(CSR)){
              return(tropicalsparse.mul(A,B,NULL,algebraType))
            }
            row_Pointer_CSR = CSR[[1]]
            col_Indices_CSR = CSR[[2]]
            values_CSR = CSR[[3]]
            x = 1

            while(x <= length(values_CSR)){
              i = row.col.Number(x,nrow(A), row_Pointer_CSR)
              j = col_Indices_CSR[x]
              result[i,j] = B + values_CSR[x]
              x = x+1
            }
          }
          return(result)

        }else if(nrow(A) == length(B)){
          CSR <- tropicalsparse.storage(A, 'csr', algebraType)
          if(is.matrix(CSR)){
            return(tropicalsparse.mul(A,B,NULL,algebraType))
          }
          row_Pointer_CSR = CSR[[1]]
          col_Indices_CSR = CSR[[2]]
          values_CSR = CSR[[3]]
          x1 = 1
          x = 1

          while(x <= nrow(A)){
            y = 1
            while(y <= ncol(A)){
              i = row.col.Number(x1,nrow(A), row_Pointer_CSR)
              j = col_Indices_CSR[x1]
              if(!(is.na(i)) && !(is.na(j))){
                if(i == x && j == y){
                  if(x1 <= length(values_CSR)){
                    result[i,j] = B[x] + values_CSR[x1]
                    x1 = x1+1
                  }
                }
              }
              y = y+1
            }
            x = x+1
          }
          return(result)

        }else{
          stop("longer object length is not a multiple of shorter object length")
        }

      }else if(store == "csc"){
        if(length(B) == 1){
          if(B != Inf){
            CSC <- tropicalsparse.storage(A, 'csc', algebraType)
            if(is.matrix(CSC)){
              return(tropicalsparse.mul(A,B,NULL,algebraType))
            }
            col_Pointer_CSC = CSC[[1]]
            row_Indices_CSC = CSC[[2]]
            values_CSC = CSC[[3]]
            x = 1
            while(x <= length(values_CSC)){
              i = row_Indices_CSC[x]
              j = row.col.Number(x, ncol(A), col_Pointer_CSC)
              result[i,j] = B + values_CSC[x]
              x = x+1
            }
          }
          return(result)

        }else if(nrow(A) == length(B)){
          CSC <- tropicalsparse.storage(A, 'csc', algebraType)
          if(is.matrix(CSC)){
            return(tropicalsparse.mul(A,B,NULL,algebraType))
          }
          col_Pointer_CSC = CSC[[1]]
          row_Indices_CSC = CSC[[2]]
          values_CSC = CSC[[3]]
          x1 = 1
          x = 1

          while(x <= ncol(A)){
            y = 1
            while(y <= nrow(A)){
              i = row_Indices_CSC[x1]
              j = row.col.Number(x1, ncol(A), col_Pointer_CSC)
              if(!(is.na(i)) && !(is.na(j))){
                if(i == y && j == x){
                  if(x1 <= length(values_CSC)){
                    result[i,j] = B[y] + values_CSC[x1]
                    x1 = x1+1
                  }
                }
              }
              y = y+1
            }
            x = x+1
          }
          return(result)

        }else{
          stop("longer object length is not a multiple of shorter object length")
        }
      }else{
        stop("Invalid storage format")
      }
    }else if(is.vector(A) && is.vector(B)){

      check.infinityV(A, algebraType)
      check.infinityV(B, algebraType)

      result = matrix(Inf, length(A), 1)
      result = as.vector(result)

      if((length(A) == 1 && length(B) > 1) || (length(B) == 1 && length(A) > 1)){
        if(length(B) == 1 && length(A) > 1){
          temp = A
          A = B
          B = temp
        }
        if(A != Inf){
          i = 1
          repeat{
            if(i>length(B))
              break
            result[i] = A + B[i]
            i = i+1
          }
        }
        return(result)

      }else{
        if(length(A) == length(B)){
          i = 1
          repeat{
            if(i>length(B))
              break
            if(A != Inf || B != Inf){
              result[i] = A[i] + B[i]
            }
            i = i+1
          }
          return(result)
        }
        stop("non-conformable arguments")
      }
    }else{
      stop("Invalid Input")
    }
  }else if(algebraType == "maxplus") {
    if(is.matrix(A) && is.matrix(B)){
      if(ncol(A) != nrow(B)){
        stop("non-conformable arguments")
      }
      if(is.null(store)){

        check.infinityM(A, algebraType)
        check.infinityM(B, algebraType)

        result = matrix(-Inf, nrow = nrow(A), ncol = ncol(B), byrow = TRUE)
        i = 1

        while(i<=nrow(A)){
          j = 1
          while(j<=ncol(B)){
            k = 1
            while(k<=ncol(A)){
              result[i,j] = max(result[i,j], A[i,k] + B[k,j])
              k = k+1
            }
            j = j+1
          }
          i = i+1
        }
        return(result)

      }else if(store == "coo"){

        counterA <- counter(A)
        counterB <- counter(B)
        COOA <- tropicalsparse.storage(A, 'coo', algebraType)
        if(is.matrix(COOA)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        row_IndicesA_COO = COOA[[1]]
        col_IndicesA_COO = COOA[[2]]
        valuesA_COO = COOA[[3]]
        COOB = tropicalsparse.storage(B, 'coo', algebraType)
        if(is.matrix(COOB)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        row_IndicesB_COO = COOB[[1]]
        col_IndicesB_COO = COOB[[2]]
        valuesB_COO = COOB[[3]]
        result = matrix(-Inf, nrow = nrow(A), ncol = ncol(B), byrow = TRUE)
        i = 1

        while(i<=counterA){
          j = 1
          while(j<=counterB)
          {
            if(col_IndicesA_COO[i]==row_IndicesB_COO[j]){
              result[row_IndicesA_COO[i],col_IndicesB_COO[j]] = max(result[row_IndicesA_COO[i],col_IndicesB_COO[j]] , (valuesA_COO[i] + valuesB_COO[j]))
            }
            if(col_IndicesA_COO[i]<row_IndicesB_COO[j]){
              break
            }
            j = j+1
          }
          i = i+1
        }
        return(result)

      }else if(store == "csr"){

        counterA <- counter(A)
        counterB <- counter(B)
        CSRA <- tropicalsparse.storage(A, 'csr', algebraType)
        if(is.matrix(CSRA)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        row_PointerA_CSR = CSRA[[1]]
        col_IndicesA_CSR = CSRA[[2]]
        valuesA_CSR = CSRA[[3]]
        CSRB <- tropicalsparse.storage(B, 'csr', algebraType)
        if(is.matrix(CSRB)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        row_PointerB_CSR = CSRB[[1]]
        col_IndicesB_CSR = CSRB[[2]]
        valuesB_CSR = CSRB[[3]]
        result = matrix(-Inf, nrow = nrow(A), ncol = ncol(B), byrow = TRUE)
        i = 1

        while(i<=counterA){
          B_Row_Number = col_IndicesA_CSR[i]
          A_Row_Number = row.col.Number(i,nrow(A),row_PointerA_CSR)
          j = row_PointerB_CSR[B_Row_Number] + 1

          while(j <= row_PointerB_CSR[B_Row_Number + 1]){
            result[A_Row_Number, col_IndicesB_CSR[j]] = max(result[A_Row_Number,col_IndicesB_CSR[j]], valuesA_CSR[i] + valuesB_CSR[j])
            j = j+1
          }
          i = i+1
        }
        return(result)

      }else if(store == "csc"){

        counterA <- counter(A)
        counterB <- counter(B)
        CSCA <- tropicalsparse.storage(A, 'csc', algebraType)
        if(is.matrix(CSCA)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        col_PointerA_CSC = CSCA[[1]]
        row_IndicesA_CSC = CSCA[[2]]
        valuesA_CSC = CSCA[[3]]
        CSCB <- tropicalsparse.storage(B, 'csc', algebraType)
        if(is.matrix(CSCB)){
          return(tropicalsparse.mul(A,B,NULL,algebraType))
        }
        col_PointerB_CSC = CSCB[[1]]
        row_IndicesB_CSC = CSCB[[2]]
        valuesB_CSC = CSCB[[3]]
        result = matrix(-Inf, nrow = nrow(A), ncol = ncol(B), byrow = TRUE)
        i = 1

        while(i<=counterA){
          A_Col_Number = row.col.Number(i,ncol(A),col_PointerA_CSC)
          j = 1
          while(j<=counterB){
            if(A_Col_Number==row_IndicesB_CSC[j]){
              B_Col_Number = row.col.Number(j,ncol(B),col_PointerB_CSC)
              result[row_IndicesA_CSC[i], B_Col_Number] = max(result[row_IndicesA_CSC[i], B_Col_Number], valuesA_CSC[i] + valuesB_CSC[j])
            }
            j = j+1
          }
          i = i+1
        }
        return(result)

      }else{
        stop("Invalid Storage Format")
      }

    }else if(is.matrix(A) && is.vector(B) || is.matrix(B) && is.vector(A)){
      if(is.matrix(B) && is.vector(A)){
        temp = A
        A = B
        B = temp
      }

      check.infinityM(A, algebraType)
      check.infinityV(B, algebraType)

      result = matrix(-Inf, nrow(A), ncol(A))
      if(is.null(store)){
        if(length(B) == 1){
          if(B != -Inf){
            i = 1
            while(i<=nrow(A)){
              j = 1
              while(j<=ncol(A)){
                result[i, j] = A[i, j] + B
                j = j+1
              }
              i = i+1
            }
          }
          return(result)

        }else if(nrow(A) == length(B)){
          i = 1
          while(i<=ncol(A)){
            j = 1
            while(j<=length(B)){
              if(A[j,i] != -Inf && B[j] != -Inf){
                result[j, i] = A[j, i] + B[j]
              }
              j = j+1
            }
            i = i+1
          }
          return(result)

        }else{
          stop("longer object length is not a multiple of shorter object length")
        }
      }else if(store == "coo"){
        if(length(B) == 1){
          if(B != -Inf){
            COO <- tropicalsparse.storage(A, 'coo', algebraType)
            if(is.matrix(COO)){
              return(tropicalsparse.mul(A,B,NULL,algebraType))
            }
            row_Indices_COO = COO[[1]]
            col_Indices_COO = COO[[2]]
            values_COO = COO[[3]]
            x = 1

            while(x <= length(values_COO)){
              i = row_Indices_COO[x]
              j = col_Indices_COO[x]
              result[i,j] = B + values_COO[x]
              x = x+1
            }
          }
          return(result)
        }else if(nrow(A) == length(B)){
          COO <- tropicalsparse.storage(A, 'coo', algebraType)
          if(is.matrix(COO)){
            return(tropicalsparse.mul(A,B,NULL,algebraType))
          }
          row_Indices_COO = COO[[1]]
          col_Indices_COO = COO[[2]]
          values_COO = COO[[3]]
          x1 = 1
          x = 1

          while(x <= nrow(A)){
            y = 1
            while(y <= ncol(A)){
              i = row_Indices_COO[x1]
              j = col_Indices_COO[x1]
              if(!(is.na(i))){
                if(i == x && j == y){
                  if(x1 <= length(values_COO)){
                    result[i,j] = B[x] + values_COO[x1]
                    x1 = x1+1
                  }
                }
              }
              y = y+1
            }
            x = x+1
          }
          return(result)

        }else{
          stop("longer object length is not a multiple of shorter object length")
        }
      }else if(store == "csr"){
        if(length(B) == 1){
          if(B != -Inf){
            CSR <- tropicalsparse.storage(A, 'csr', algebraType)
            if(is.matrix(CSR)){
              return(tropicalsparse.mul(A,B,NULL,algebraType))
            }
            row_Pointer_CSR = CSR[[1]]
            col_Indices_CSR = CSR[[2]]
            values_CSR = CSR[[3]]
            x = 1

            while(x <= length(values_CSR)){
              i = row.col.Number(x,nrow(A), row_Pointer_CSR)
              j = col_Indices_CSR[x]
              result[i,j] = B + values_CSR[x]
              x = x+1
            }
          }
          return(result)

        }else if(nrow(A) == length(B)){
          CSR <- tropicalsparse.storage(A, 'csr', algebraType)
          if(is.matrix(CSR)){
            return(tropicalsparse.mul(A,B,NULL,algebraType))
          }
          row_Pointer_CSR = CSR[[1]]
          col_Indices_CSR = CSR[[2]]
          values_CSR = CSR[[3]]
          x1 = 1
          x = 1

          while(x <= nrow(A)){
            y = 1
            while(y <= ncol(A)){
              i = row.col.Number(x1,nrow(A), row_Pointer_CSR)
              j = col_Indices_CSR[x1]
              if(!(is.na(i)) && !(is.na(j))){
                if(i == x && j == y){
                  if(x1 <= length(values_CSR)){
                    result[i,j] = B[x] + values_CSR[x1]
                    x1 = x1+1
                  }
                }
              }
              y = y+1
            }
            x = x+1
          }
          return(result)

        }else{
          stop("longer object length is not a multiple of shorter object length")
        }

      }else if(store == "csc"){
        if(length(B) == 1){
          if(B != -Inf){
            CSC <- tropicalsparse.storage(A, 'csc', algebraType)
            if(is.matrix(CSC)){
              return(tropicalsparse.mul(A,B,NULL,algebraType))
            }
            col_Pointer_CSC = CSC[[1]]
            row_Indices_CSC = CSC[[2]]
            values_CSC = CSC[[3]]
            x = 1
            while(x <= length(values_CSC)){
              i = row_Indices_CSC[x]
              j = row.col.Number(x, ncol(A), col_Pointer_CSC)
              result[i,j] = B + values_CSC[x]
              x = x+1
            }
          }
          return(result)

        }else if(nrow(A) == length(B)){
          CSC <- tropicalsparse.storage(A, 'csc', algebraType)
          if(is.matrix(CSC)){
            return(tropicalsparse.mul(A,B,NULL,algebraType))
          }
          col_Pointer_CSC = CSC[[1]]
          row_Indices_CSC = CSC[[2]]
          values_CSC = CSC[[3]]
          x1 = 1
          x = 1

          while(x <= ncol(A)){
            y = 1
            while(y <= nrow(A)){
              i = row_Indices_CSC[x1]
              j = row.col.Number(x1, ncol(A), col_Pointer_CSC)
              if(!(is.na(i)) && !(is.na(j))){
                if(i == y && j == x){
                  if(x1 <= length(values_CSC)){
                    result[i,j] = B[y] + values_CSC[x1]
                    x1 = x1+1
                  }
                }
              }
              y = y+1
            }
            x = x+1
          }
          return(result)

        }else{
          stop("longer object length is not a multiple of shorter object length")
        }
      }else{
        stop("Invalid storage format")
      }
    }else if(is.vector(A) && is.vector(B)){

      check.infinityV(A, algebraType)
      check.infinityV(B, algebraType)

      result = matrix(-Inf, length(A), 1)
      result = as.vector(result)

      if((length(A) == 1 && length(B) > 1) || (length(B) == 1 && length(A) > 1)){
        if(length(B) == 1 && length(A) > 1){
          temp = A
          A = B
          B = temp
        }
        if(A != -Inf){
          i = 1
          repeat{
            if(i>length(B))
              break
            result[i] = A + B[i]
            i = i+1
          }
        }
        return(result)

      }else{
        if(length(A) == length(B)){
          i = 1
          repeat{
            if(i>length(B))
              break
            if(A != -Inf || B != -Inf){
              result[i] = A[i] + B[i]
            }
            i = i+1
          }
          return(result)
        }
        stop("non-conformable arguments")
      }
    }else{
      stop("Invalid Input")
    }
  }else{
    stop("Invalid algebra type")
  }
}
