#error calculation

#' Root Mean Square Error
#'
#' This function calculates the root mean square error and checks that the idenifiers match. Note that if n=1, then the root mean square error is equivalent to the average error.
#' @param predicted the predictions including the first column with a character vector with idenifiers to check that the rows of predicted and actual match 
#' @param actual the real result including the first column with a character vector
#'
#' @return the Root Mean Square Error
#' @export
#'
#' @examples
#' predict1 = data.frame(state = c("TX", "GA", "OR"), dem = c(.45, .49, .63), rep = c(.55, .51, .37))
#' actual1 = data.frame(state = c("TX", "GA", "OR"), dem = c(.42, .48, .62), rep = c(.58, .52, .38))
#' rmse(predict1, actual1)
#' predict2 = data.frame(state = c("OK", "MS", "WY"), dem = c(.42, .44, .3), rep = c(.56, .55, .7), other = c(0.02, 0.01, 0))
#' actual2 = data.frame(state = c("OK", "MS", "WY"), dem = c(.41, .38, .3), rep = c(.55, .58, .69), other = c(0.04, 0.04, 0.01))
#' rmse(predict2, actual2)

rmse = function(predicted, actual){
  if(!identical(predicted[, 1], actual[ , 1])){
    warning("Predicted names do not match names of actual data")
  }
  if(!is.factor(predicted[ , 1])){
    stop("First column of predicted is not a character vector")
  }
  if(!is.factor(actual[ , 1])){
    stop("First column of actual is not a character vector")
  }
  if(nrow(predicted) != nrow(actual) | ncol(predicted) != ncol(actual)){
    stop("Incompatiable dimensions of predicted and actual")
  }
  correction = ifelse(ncol(actual)>2, 2, 1)
  return(sum((predicted[, -1] - actual[, -1])^2)/correction)
}

#' Average Error by State
#' This function calculates the average error and checks that the idenifiers match.
#' @param predicted the predictions including the first column with a character vector with idenifiers to check that the rows of predicted and actual match 
#' @param actual the real result including the first column with a character vector
#'
#' @return This function returns the average error for each state.
#' @export
#'
#' @examples
#' predict1 = data.frame(state = c("TX", "GA", "OR"), dem = c(.45, .49, .63), rep = c(.55, .51, .37))
#' actual1 = data.frame(state = c("TX", "GA", "OR"), dem = c(.42, .48, .62), rep = c(.58, .52, .38))
#' average_error(predict1, actual1)
#' predict2 = data.frame(state = c("OK", "MS", "WY"), dem = c(.42, .44, .3), rep = c(.56, .55, .7), other = c(0.02, 0.01, 0))
#' actual2 = data.frame(state = c("OK", "MS", "WY"), dem = c(.41, .38, .3), rep = c(.55, .58, .69), other = c(0.04, 0.04, 0.01))
#' average_error(predict2, actual2)
average_error = function(predicted, actual){
  if(!identical(predicted[, 1], actual[ , 1])){
    warning("Predicted names do not match names of actual data")
  }
  if(!is.factor(predicted[ , 1])){
    stop("First column of predicted is not a character vector")
  }
  if(!is.factor(actual[ , 1])){
    stop("First column of actual is not a character vector")
  }
  if(nrow(predicted) != nrow(actual) | ncol(predicted) != ncol(actual)){
    stop("Incompatiable dimensions of predicted and actual")
  }
  return(sum(predicted[, -1] - actual[, -1]))
}

#' Average Error by State
#'
#' This function calculates the root mean square error for each state over a single or multiple years and checks that the idenifiers match.
#' @param predicted the predictions including the first column with a character vector with idenifiers to check that the rows of predicted and actual match 
#' @param actual the real result including the first column with a character vector
#'
#' @return a data frame with the state and the average error for that state 
#' @export
#'
#' @examples
#' predict1 = data.frame(state = c("TX", "GA", "OR"), dem1 = c(.45, .49, .63), dem2 = c(.55, .51, .37))
#' actual1 = data.frame(state = c("TX", "GA", "OR"), dem1 = c(.42, .48, .62), dem2 = c(.58, .52, .38))
#' ae_by_state(predict1, actual1)
#' predict2 = data.frame(state = c("OK", "MS", "WY"), dem1 = c(.42, .44, .3))
#' actual2 = data.frame(state = c("OK", "MS", "WY"), dem1 = c(.41, .38, .3))
#' ae_by_state(predict2, actual2)
ae_by_state = function(predicted, actual){
  if(!identical(predicted[, 1], actual[ , 1])){
    warning("Predicted names do not match names of actual data")
  }
  if(!is.factor(predicted[ , 1])){
    stop("First column of predicted is not a character vector")
  }
  if(!is.factor(actual[ , 1])){
    stop("First column of actual is not a character vector")
  }
  if(nrow(predicted) != nrow(actual) | ncol(predicted) != ncol(actual)){
    stop("Incompatiable dimensions of predicted and actual")
  }
  ae = matrix(0, nrow = nrow(predict), ncol = (ncol(predict)-1))
  for(i in 1:nrow(predicted)){
    for(j in 1:(ncol(predict) - 1)){
      ae[i, j] = average_error(predicted[i,j ], actual[i, j])
    }
    tae = rowMeans(ae)
  }
  return(data.frame("state" = predicted[, 1]), "Average Error by predict" = ae, "Total Average Error"  = tae)
  
}
