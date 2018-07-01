#' @title Addition With or Without Storage Techniques
#'
#' @description \code{tropicalsparse.add} function adds the provided inputs in Tropical Algebra based on type of
#' Tropical Algebra.
#'
#' @param A is matrix or vector.
#' @param B is matrix or vector.
#' @param store is storage technique.
#' @param algebraType is string input that can be \code{minplus} or \code{maxplus}.
#'
#' @details The compulsory inputs of the function \code{tropicalsparse.add} are \code{A}, \code{B} and
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
#' @return Addition of \code{A} and \code{B} in Tropical Algebra.
#'
#' @seealso
#' \code{\link{tropicalsparse.mul}}, \code{\link{tropicalsparse.storage}}
#'
#' @examples
#' a <- matrix(data = c(2, Inf, Inf, 0, Inf, Inf, Inf, 10, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' b <- matrix(data = c(Inf, Inf, 4, Inf, -0.3, Inf, Inf, 2, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' tropicalsparse.add(a, b, 'csr', 'minplus')
#'
#'#     [,1] [,2]  [,3]
#'# [1,]    2  Inf    4
#'# [2,]    0 -0.3  Inf
#'# [3,]  Inf  2.0  Inf
#'
#'# also
#'
#' a <- matrix(data = c(5, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 10, 2),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' b <- matrix(data = c(-Inf, -Inf, 3, -Inf, -0.5, -Inf, 1.1, -Inf, -Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' tropicalsparse.add(a, b, 'coo', 'maxplus')
#'
#'#      [,1] [,2] [,3]
#'# [1,]  5.0 -Inf    3
#'# [2,] -Inf -0.5 -Inf
#'# [3,]  1.1 10.0    2
#'
#'# also
#'
#' a <- matrix(data = c(2, Inf, Inf, 0, Inf, Inf, Inf, 2, Inf),
#' nrow = 3, ncol = 3, byrow = TRUE)
#'
#' b <- c(Inf, 0, 10)
#'
#' tropicalsparse.add(a, b, algebraType = 'minplus')
#'
#'#      [,1] [,2] [,3]
#'# [1,]    2  Inf  Inf
#'# [2,]    0    0    0
#'# [3,]   10    2   10
#'
#' @export
#'

tropicalsparse.add <- function(A, B, store = NULL, algebraType) {
  if (algebraType == "minplus") {
    if (is.matrix(A) && is.matrix(B)) {
      if(nrow(A) != nrow(B) || ncol(A) != ncol(B)){
        stop("non-conformable arguments")
      }

      if (is.null(store)) {
        check.infinityM(A, algebraType)
        check.infinityM(B, algebraType)

        result = matrix(Inf, nrow = nrow(A), ncol = ncol(A), byrow = TRUE)
        i = 1
        repeat {
          if (i > nrow(A))
            break
          j = 1
          repeat {
            if (j > ncol(A))
              break
            result[i, j] = min(A[i, j], B[i, j])
            j = j + 1
          }
          i = i + 1
        }
        return(result)

      } else if (store == "coo") {
        counterA <- counter(A)
        counterB <- counter(B)
        COOA <- tropicalsparse.storage(A, 'coo', algebraType)
        if(is.matrix(COOA)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        row_IndicesA_COO = COOA[[1]]
        col_IndicesA_COO = COOA[[2]]
        valuesA_COO = COOA[[3]]
        COOB = tropicalsparse.storage(B, 'coo', algebraType)
        if(is.matrix(COOB)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        row_IndicesB_COO = COOB[[1]]
        col_IndicesB_COO = COOB[[2]]
        valuesB_COO = COOB[[3]]
        result = matrix(Inf, nrow = nrow(A), ncol = ncol(A), byrow = TRUE)

        if (counterA > counterB) {
          max = counterA
        } else{
          max = counterB
        }
        x = 1
        repeat {
          if (x > max) {
            break
          }
          if (x <= counterA) {
            i = row_IndicesA_COO[x]
            j = col_IndicesA_COO[x]
            y = 1
            fnd = 0
            while (y <= counterB) {
              if (i == row_IndicesB_COO[y] && j == col_IndicesB_COO[y]) {
                result[i, j] = min(valuesA_COO[x], valuesB_COO[y])
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesA_COO[x]
            }
          }

          if (x <= counterB) {
            i = row_IndicesB_COO[x]
            j = col_IndicesB_COO[x]
            y = 1
            fnd = 0

            while (y <= counterA) {
              if (i == row_IndicesA_COO[y] && j == col_IndicesA_COO[y]) {
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesB_COO[x]
            }
          }
          x = x + 1
        }
        return(result)

      } else if (store == "csr") {
        counterA <- counter(A)
        counterB <- counter(B)
        CSRA <- tropicalsparse.storage(A, 'csr', algebraType)
        if(is.matrix(CSRA)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        row_PointerA_CSR = CSRA[[1]]
        col_IndicesA_CSR = CSRA[[2]]
        valuesA_CSR = CSRA[[3]]
        CSRB <- tropicalsparse.storage(B, 'csr', algebraType)
        if(is.matrix(CSRB)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        row_PointerB_CSR = CSRB[[1]]
        col_IndicesB_CSR = CSRB[[2]]
        valuesB_CSR = CSRB[[3]]
        result = matrix(Inf, nrow = nrow(A), ncol = ncol(A), byrow = TRUE)
        if (counterA > counterB) {
          max = counterA
        } else{
          max = counterB
        }
        x = 1
        repeat {
          if (x > max) {
            break
          }
          if (x <= counterA) {
            i = row.col.Number(x, nrow(A), row_PointerA_CSR)
            j = col_IndicesA_CSR[x]
            y = 1
            fnd = 0
            while (y <= counterB) {
              if (i == row.col.Number(y, nrow(B), row_PointerB_CSR) && j == col_IndicesB_CSR[y]) {
                result[i, j] = min(valuesA_CSR[x], valuesB_CSR[y])
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesA_CSR[x]
            }
          }

          if (x <= counterB) {
            i = row.col.Number(x, nrow(B), row_PointerB_CSR)
            j = col_IndicesB_CSR[x]
            y = 1
            fnd = 0

            while (y <= counterA) {
              if (i == row.col.Number(y, nrow(A), row_PointerA_CSR) && j == col_IndicesA_CSR[y]) {
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesB_CSR[x]
            }
          }
          x = x + 1
        }
        return(result)

      } else if (store == "csc") {
        counterA <- counter(A)
        counterB <- counter(B)
        CSCA <- tropicalsparse.storage(A, 'csc', algebraType)
        if(is.matrix(CSCA)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        col_PointerA_CSC = CSCA[[1]]
        row_IndicesA_CSC = CSCA[[2]]
        valuesA_CSC = CSCA[[3]]
        CSCB <- tropicalsparse.storage(B, 'csc', algebraType)
        if(is.matrix(CSCB)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        col_PointerB_CSC = CSCB[[1]]
        row_IndicesB_CSC = CSCB[[2]]
        valuesB_CSC = CSCB[[3]]
        result = matrix(Inf, nrow = nrow(A), ncol = ncol(A), byrow = TRUE)

        if (counterA > counterB) {
          max = counterA
        } else{
          max = counterB
        }
        x = 1
        repeat {
          if (x > max) {
            break
          }
          if (x <= counterA) {
            i = row_IndicesA_CSC[x]
            j = row.col.Number(x, ncol(A), col_PointerA_CSC)
            y = 1
            fnd = 0

            while (y <= counterB) {
              if (i == row_IndicesB_CSC[y] && j == row.col.Number(y, ncol(B), col_PointerB_CSC)) {
                result[i, j] = min(valuesA_CSC[x], valuesB_CSC[y])
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesA_CSC[x]
            }
          }

          if (x <= counterB) {
            i = row_IndicesB_CSC[x]
            j = row.col.Number(x, ncol(B), col_PointerB_CSC)
            y = 1
            fnd = 0
            while (y <= counterA) {
              if (i == row_IndicesA_CSC[y] && j == row.col.Number(y, ncol(A), col_PointerA_CSC)) {
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesB_CSC[x]
            }
          }
          x = x + 1
        }
        return(result)

      } else{
        stop("Invalid Storage Format")
      }
    } else if (is.matrix(A) && is.vector(B) || is.matrix(B) && is.vector(A)) {
      if (is.matrix(B) && is.vector(A)) {
        temp = A
        A = B
        B = temp
      }
      if (is.null(store)) {

        check.infinityM(A, algebraType)
        check.infinityV(B, algebraType)

        if (length(B) == 1) {
          if (!is.infinite(B)) {
            i = 1
            while (i <= nrow(A)) {
              j = 1
              while (j <= ncol(A)) {
                A[i, j] = min(A[i, j], B)
                j = j + 1
              }
              i = i + 1
            }
          }
          return(A)

        } else if (nrow(A) == length(B)) {
          i = 1

          while (i <= ncol(A)) {
            j = 1
            while (j <= length(B)) {
              if (A[j, i] != Inf || B[j] != Inf) {
                A[j, i] = min(A[j, i], B[j])
              }
              j = j + 1
            }
            i = i + 1
          }
          return(A)

        } else{
          stop("longer object length is not a multiple of shorter object length")
        }
      } else if (store == "coo") {
        if (length(B) == 1) {
          if (B != Inf) {
            result = matrix(B, nrow = nrow(A), ncol = ncol(A))
            COO <- tropicalsparse.storage(A, 'coo', algebraType)
            if(is.matrix(COO)){
              return(tropicalsparse.add(A,B,NULL,algebraType))
            }
            row_Indices_COO = COO[[1]]
            col_Indices_COO = COO[[2]]
            values_COO = COO[[3]]
            x = 1

            while (x <= length(values_COO)) {
              i = row_Indices_COO[x]
              j = col_Indices_COO[x]
              result[i, j] = min(B, values_COO[x])
              x = x + 1
            }
            return(result)
          }
          return(A)
        } else if (nrow(A) == length(B)) {
          COO <- tropicalsparse.storage(A, 'coo', algebraType)
          if(is.matrix(COO)){
            return(tropicalsparse.add(A,B,NULL,algebraType))
          }
          row_Indices_COO = COO[[1]]
          col_Indices_COO = COO[[2]]
          values_COO = COO[[3]]
          x1 = 1
          x = 1

          while (x <= nrow(A)) {
            y = 1
            while (y <= ncol(A)) {
              i = row_Indices_COO[x1]
              j = col_Indices_COO[x1]
              if (!(is.na(i))) {
                if (i == x && j == y) {
                  if (x1 <= length(values_COO)) {
                    A[i, j] = min(B[x], values_COO[x1])
                    x1 = x1 + 1
                  }
                } else{
                  A[x, y] = B[x]
                }
              } else{
                A[x, y] = B[x]
              }
              y = y + 1
            }
            x = x + 1
          }
          return(A)

        } else{
          stop("longer object length is not a multiple of shorter object length")
        }
      } else if (store == "csr") {
        if (length(B) == 1) {
          if (B != Inf) {
            result = matrix(B, nrow = nrow(A), ncol = ncol(A))
            CSR <- tropicalsparse.storage(A, 'csr', algebraType)
            if(is.matrix(CSR)){
              return(tropicalsparse.add(A,B,NULL,algebraType))
            }
            row_Pointer_CSR = CSR[[1]]
            col_Indices_CSR = CSR[[2]]
            values_CSR = CSR[[3]]
            x = 1

            while (x <= length(values_CSR)) {
              i = row.col.Number(x, nrow(A), row_Pointer_CSR)
              j = col_Indices_CSR[x]
              result[i, j] = min(B, values_CSR[x])
              x = x + 1
            }
            return(result)
          }
          return(A)

        } else if (nrow(A) == length(B)) {
          CSR <- tropicalsparse.storage(A, 'csr', algebraType)
          if(is.matrix(CSR)){
            return(tropicalsparse.add(A,B,NULL,algebraType))
          }
          row_Pointer_CSR = CSR[[1]]
          col_Indices_CSR = CSR[[2]]
          values_CSR = CSR[[3]]
          x1 = 1
          x = 1

          while (x <= nrow(A)) {
            y = 1
            while (y <= ncol(A)) {
              i = row.col.Number(x1, nrow(A), row_Pointer_CSR)
              j = col_Indices_CSR[x1]
              if (!(is.na(i)) && !(is.na(j))) {
                if (i == x && j == y) {
                  if (x1 <= length(values_CSR)) {
                    A[i, j] = min(B[x], values_CSR[x1])
                    x1 = x1 + 1
                  }
                } else{
                  A[x, y] = B[x]
                }
              } else{
                A[x, y] = B[x]
              }
              y = y + 1
            }
            x = x + 1
          }
          return(A)

        } else{
          stop("longer object length is not a multiple of shorter object length")
        }

      } else if (store == "csc") {
        if (length(B) == 1) {
          if (B != Inf) {
            result = matrix(B, nrow = nrow(A), ncol = ncol(A))
            CSC <- tropicalsparse.storage(A, 'csc', algebraType)
            if(is.matrix(CSC)){
              return(tropicalsparse.add(A,B,NULL,algebraType))
            }
            col_Pointer_CSC = CSC[[1]]
            row_Indices_CSC = CSC[[2]]
            values_CSC = CSC[[3]]
            x = 1
            while (x <= length(values_CSC)) {
              i = row_Indices_CSC[x]
              j = row.col.Number(x, ncol(A), col_Pointer_CSC)
              result[i, j] = min(B, values_CSC[x])
              x = x + 1
            }
            return(result)
          }
          return(A)

        } else if (nrow(A) == length(B)) {
          CSC <- tropicalsparse.storage(A, 'csc', algebraType)
          if(is.matrix(CSC)){
            return(tropicalsparse.add(A,B,NULL,algebraType))
          }
          col_Pointer_CSC = CSC[[1]]
          row_Indices_CSC = CSC[[2]]
          values_CSC = CSC[[3]]
          x1 = 1
          x = 1

          while (x <= ncol(A)) {
            y = 1
            while (y <= nrow(A)) {
              i = row_Indices_CSC[x1]
              j = row.col.Number(x1, ncol(A), col_Pointer_CSC)
              if (!(is.na(i)) && !(is.na(j))) {
                if (i == y && j == x) {
                  if (x1 <= length(values_CSC)) {
                    A[i, j] = min(B[y], values_CSC[x1])
                    x1 = x1 + 1
                  }
                } else{
                  A[y, x] = B[y]
                }
              } else{
                A[y, x] = B[y]
              }
              y = y + 1
            }
            x = x + 1
          }
          return(A)

        } else{
          stop("longer object length is not a multiple of shorter object length")
        }
      } else{
        stop("Invalid storage format")
      }
    } else if (is.vector(A) && is.vector(B)) {

      if(!is.null(store)){
        stop("Store argument must be NULL in Vector-Vector addition.")
      }

      check.infinityV(A, algebraType)
      check.infinityV(B, algebraType)

      if ((length(A) == 1 && length(B) > 1) || (length(B) == 1 && length(A) > 1)) {
        if (length(B) == 1 && length(A) > 1) {
          temp = A
          A = B
          B = temp
        }
        if (A != Inf) {
          i = 1
          repeat {
            if (i > length(B))
              break
            B[i] = min(A, B[i])
            i = i + 1
          }
        }
        return(B)

      } else{
        if (length(A) == length(B)) {
          i = 1
          repeat {
            if (i > length(B))
              break
            if (A != Inf || B != Inf) {
              A[i] = min(A[i], B[i])
            }
            i = i + 1
          }
          return(A)
        }
        stop("non-conformable arguments")
      }
    } else{
      stop("Invalid Input")
    }
  } else if (algebraType == "maxplus") {
    if (is.matrix(A) && is.matrix(B)) {

      if(nrow(A) != nrow(B) || ncol(A) != ncol(B)){
        stop("non-conformable arguments")
      }

      if (is.null(store)) {
        check.infinityM(A, algebraType)
        check.infinityM(B, algebraType)

        result = matrix(-Inf,nrow = nrow(A),ncol = ncol(A),byrow = TRUE)
        i = 1
        repeat {
          if (i > nrow(A))
            break
          j = 1
          repeat {
            if (j > ncol(A))
              break
            result[i, j] = max(A[i, j], B[i, j])
            j = j + 1
          }
          i = i + 1
        }
        return(result)

      } else if (store == "coo") {

        counterA <- counter(A)
        counterB <- counter(B)
        COOA <- tropicalsparse.storage(A, 'coo', algebraType)
        if(is.matrix(COOA)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        row_IndicesA_COO = COOA[[1]]
        col_IndicesA_COO = COOA[[2]]
        valuesA_COO = COOA[[3]]
        COOB = tropicalsparse.storage(B, 'coo', algebraType)
        if(is.matrix(COOB)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        row_IndicesB_COO = COOB[[1]]
        col_IndicesB_COO = COOB[[2]]
        valuesB_COO = COOB[[3]]
        result = matrix(-Inf,nrow = nrow(A),ncol = ncol(A),byrow = TRUE)

        if (counterA > counterB) {
          max = counterA
        } else{
          max = counterB
        }
        x = 1
        repeat {
          if (x > max) {
            break
          }
          if (x <= counterA) {
            i = row_IndicesA_COO[x]
            j = col_IndicesA_COO[x]
            y = 1
            fnd = 0
            while (y <= counterB) {
              if (i == row_IndicesB_COO[y] && j == col_IndicesB_COO[y]) {
                result[i, j] = max(valuesA_COO[x], valuesB_COO[y])
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesA_COO[x]
            }
          }

          if (x <= counterB) {
            i = row_IndicesB_COO[x]
            j = col_IndicesB_COO[x]
            y = 1
            fnd = 0

            while (y <= counterA) {
              if (i == row_IndicesA_COO[y] && j == col_IndicesA_COO[y]) {
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesB_COO[x]
            }
          }
          x = x + 1
        }
        return(result)

      } else if (store == "csr") {

        counterA <- counter(A)
        counterB <- counter(B)
        CSRA <- tropicalsparse.storage(A, 'csr', algebraType)
        if(is.matrix(CSRA)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        row_PointerA_CSR = CSRA[[1]]
        col_IndicesA_CSR = CSRA[[2]]
        valuesA_CSR = CSRA[[3]]
        CSRB <- tropicalsparse.storage(B, 'csr', algebraType)
        if(is.matrix(CSRB)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        row_PointerB_CSR = CSRB[[1]]
        col_IndicesB_CSR = CSRB[[2]]
        valuesB_CSR = CSRB[[3]]
        result = matrix(-Inf,nrow = nrow(A),ncol = ncol(A),byrow = TRUE)
        if (counterA > counterB) {
          max = counterA
        } else{
          max = counterB
        }
        x = 1
        repeat {
          if (x > max) {
            break
          }
          if (x <= counterA) {
            i = row.col.Number(x, nrow(A), row_PointerA_CSR)
            j = col_IndicesA_CSR[x]
            y = 1
            fnd = 0
            while (y <= counterB) {
              if (i == row.col.Number(y, nrow(B), row_PointerB_CSR) && j == col_IndicesB_CSR[y]) {
                result[i, j] = max(valuesA_CSR[x], valuesB_CSR[y])
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesA_CSR[x]
            }
          }

          if (x <= counterB) {
            i = row.col.Number(x, nrow(B), row_PointerB_CSR)
            j = col_IndicesB_CSR[x]
            y = 1
            fnd = 0

            while (y <= counterA) {
              if (i == row.col.Number(y, nrow(A), row_PointerA_CSR) && j == col_IndicesA_CSR[y]) {
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesB_CSR[x]
            }
          }
          x = x + 1
        }
        return(result)

      } else if (store == "csc") {

        counterA <- counter(A)
        counterB <- counter(B)
        CSCA <- tropicalsparse.storage(A, 'csc', algebraType)
        if(is.matrix(CSCA)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        col_PointerA_CSC = CSCA[[1]]
        row_IndicesA_CSC = CSCA[[2]]
        valuesA_CSC = CSCA[[3]]
        CSCB <- tropicalsparse.storage(B, 'csc', algebraType)
        if(is.matrix(CSCB)){
          return(tropicalsparse.add(A,B,NULL,algebraType))
        }
        col_PointerB_CSC = CSCB[[1]]
        row_IndicesB_CSC = CSCB[[2]]
        valuesB_CSC = CSCB[[3]]
        result = matrix(-Inf,nrow = nrow(A),ncol = ncol(A),byrow = TRUE)

        if (counterA > counterB) {
          max = counterA
        } else{
          max = counterB
        }
        x = 1
        repeat {
          if (x > max) {
            break
          }
          if (x <= counterA) {
            i = row_IndicesA_CSC[x]
            j = row.col.Number(x, ncol(A), col_PointerA_CSC)
            y = 1
            fnd = 0

            while (y <= counterB) {
              if (i == row_IndicesB_CSC[y] && j == row.col.Number(y, ncol(B), col_PointerB_CSC)) {
                result[i, j] = max(valuesA_CSC[x], valuesB_CSC[y])
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesA_CSC[x]
            }
          }

          if (x <= counterB) {
            i = row_IndicesB_CSC[x]
            j = row.col.Number(x, ncol(B), col_PointerB_CSC)
            y = 1
            fnd = 0
            while (y <= counterA) {
              if (i == row_IndicesA_CSC[y] && j == row.col.Number(y, ncol(A), col_PointerA_CSC)) {
                fnd = 1
                break
              }
              y = y + 1
            }
            if (fnd != 1) {
              result[i, j] = valuesB_CSC[x]
            }
          }
          x = x + 1
        }
        return(result)

      } else{
        stop("Invalid Storage Format")
      }
    } else if (is.matrix(A) && is.vector(B) || is.matrix(B) && is.vector(A)) {
      if (is.matrix(B) && is.vector(A)) {
        temp = A
        A = B
        B = temp
      }
      if (is.null(store)) {

        check.infinityM(A, algebraType)
        check.infinityV(B, algebraType)

        if (length(B) == 1) {
          if (B != -Inf) {
            i = 1
            while (i <= nrow(A)) {
              j = 1
              while (j <= ncol(A)) {
                A[i, j] = max(A[i, j], B)
                j = j + 1
              }
              i = i + 1
            }
          }
          return(A)

        } else if (nrow(A) == length(B)) {
          i = 1

          while (i <= ncol(A)) {
            j = 1
            while (j <= length(B)) {
              if (A[j, i] != -Inf || B[j] != -Inf) {
                A[j, i] = max(A[j, i], B[j])
              }
              j = j + 1
            }
            i = i + 1
          }
          return(A)

        } else{
          stop("longer object length is not a multiple of shorter object length")
        }
      } else if (store == "coo") {
        if (length(B) == 1) {
          if (B != -Inf) {
            result = matrix(B, nrow = nrow(A), ncol = ncol(A))
            COO <- tropicalsparse.storage(A, 'coo', algebraType)
            if(is.matrix(COO)){
              return(tropicalsparse.add(A,B,NULL,algebraType))
            }
            row_Indices_COO = COO[[1]]
            col_Indices_COO = COO[[2]]
            values_COO = COO[[3]]
            x = 1

            while (x <= length(values_COO)) {
              i = row_Indices_COO[x]
              j = col_Indices_COO[x]
              result[i, j] = max(B, values_COO[x])
              x = x + 1
            }
            return(result)
          }
          return(A)
        } else if (nrow(A) == length(B)) {
          COO <- tropicalsparse.storage(A, 'coo', algebraType)
          if(is.matrix(COO)){
            return(tropicalsparse.add(A,B,NULL,algebraType))
          }
          row_Indices_COO = COO[[1]]
          col_Indices_COO = COO[[2]]
          values_COO = COO[[3]]
          x1 = 1
          x = 1

          while (x <= nrow(A)) {
            y = 1
            while (y <= ncol(A)) {
              i = row_Indices_COO[x1]
              j = col_Indices_COO[x1]
              if (!(is.na(i))) {
                if (i == x && j == y) {
                  if (x1 <= length(values_COO)) {
                    A[i, j] = max(B[x], values_COO[x1])
                    x1 = x1 + 1
                  }
                } else{
                  A[x, y] = B[x]
                }
              } else{
                A[x, y] = B[x]
              }
              y = y + 1
            }
            x = x + 1
          }
          return(A)

        } else{
          stop("longer object length is not a multiple of shorter object length")
        }
      } else if (store == "csr") {
        if (length(B) == 1) {
          if (B != -Inf) {
            result = matrix(B, nrow = nrow(A), ncol = ncol(A))
            CSR <- tropicalsparse.storage(A, 'csr', algebraType)
            if(is.matrix(CSR)){
              return(tropicalsparse.add(A,B,NULL,algebraType))
            }
            row_Pointer_CSR = CSR[[1]]
            col_Indices_CSR = CSR[[2]]
            values_CSR = CSR[[3]]
            x = 1

            while (x <= length(values_CSR)) {
              i = row.col.Number(x, nrow(A), row_Pointer_CSR)
              j = col_Indices_CSR[x]
              result[i, j] = max(B, values_CSR[x])
              x = x + 1
            }
            return(result)
          }
          return(A)

        } else if (nrow(A) == length(B)) {
          CSR <- tropicalsparse.storage(A, 'csr', algebraType)
          if(is.matrix(CSR)){
            return(tropicalsparse.add(A,B,NULL,algebraType))
          }
          row_Pointer_CSR = CSR[[1]]
          col_Indices_CSR = CSR[[2]]
          values_CSR = CSR[[3]]
          x1 = 1
          x = 1

          while (x <= nrow(A)) {
            y = 1
            while (y <= ncol(A)) {
              i = row.col.Number(x1, nrow(A), row_Pointer_CSR)
              j = col_Indices_CSR[x1]
              if (!(is.na(i)) && !(is.na(j))) {
                if (i == x && j == y) {
                  if (x1 <= length(values_CSR)) {
                    A[i, j] = max(B[x], values_CSR[x1])
                    x1 = x1 + 1
                  }
                } else{
                  A[x, y] = B[x]
                }
              } else{
                A[x, y] = B[x]
              }
              y = y + 1
            }
            x = x + 1
          }
          return(A)

        } else{
          stop("longer object length is not a multiple of shorter object length")
        }

      } else if (store == "csc") {
        if (length(B) == 1) {
          if (B != -Inf) {
            result = matrix(B, nrow = nrow(A), ncol = ncol(A))
            CSC <- tropicalsparse.storage(A, 'csc', algebraType)
            if(is.matrix(CSC)){
              return(tropicalsparse.add(A,B,NULL,algebraType))
            }
            col_Pointer_CSC = CSC[[1]]
            row_Indices_CSC = CSC[[2]]
            values_CSC = CSC[[3]]
            x = 1
            while (x <= length(values_CSC)) {
              i = row_Indices_CSC[x]
              j = row.col.Number(x, ncol(A), col_Pointer_CSC)
              result[i, j] = max(B, values_CSC[x])
              x = x + 1
            }
            return(result)
          }
          return(A)

        } else if (nrow(A) == length(B)) {
          CSC <- tropicalsparse.storage(A, 'csc', algebraType)
          if(is.matrix(CSC)){
            return(tropicalsparse.add(A,B,NULL,algebraType))
          }
          col_Pointer_CSC = CSC[[1]]
          row_Indices_CSC = CSC[[2]]
          values_CSC = CSC[[3]]
          x1 = 1
          x = 1

          while (x <= ncol(A)) {
            y = 1
            while (y <= nrow(A)) {
              i = row_Indices_CSC[x1]
              j = row.col.Number(x1, ncol(A), col_Pointer_CSC)
              if (!(is.na(i)) && !(is.na(j))) {
                if (i == y && j == x) {
                  if (x1 <= length(values_CSC)) {
                    A[i, j] = max(B[y], values_CSC[x1])
                    x1 = x1 + 1
                  }
                } else{
                  A[y, x] = B[y]
                }
              } else{
                A[y, x] = B[y]
              }
              y = y + 1
            }
            x = x + 1
          }
          return(A)

        } else{
          stop("longer object length is not a multiple of shorter object length")
        }
      } else{
        stop("Invalid storage format")
      }
    } else if (is.vector(A) && is.vector(B)) {

      check.infinityV(A, algebraType)
      check.infinityV(B, algebraType)

      if ((length(A) == 1 && length(B) > 1) || (length(B) == 1 && length(A) > 1)) {
        if (length(B) == 1 && length(A) > 1) {
          temp = A
          A = B
          B = temp
        }

        if (A != -Inf) {
          i = 1
          repeat {
            if (i > length(B))
              break
            B[i] = max(A, B[i])
            i = i + 1
          }
        }
        return(B)

      } else{
        if (length(A) == length(B)) {
          i = 1
          repeat {
            if (i > length(B))
              break
            if (A != -Inf || B != -Inf) {
              A[i] = max(A[i], B[i])
            }
            i = i + 1
          }
          return(A)
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
