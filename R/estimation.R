#file to get the 

#' MLE estimates of a Normal Distribution
#' Fits a normal 
#' @param data - a vector of data
#'
#' @return a list with the elements mean with the mean of the distribution and 
#' @export
#'
#' @examples
#' set.seed(1)
#' data1  = rnorm(100)
#' fitnormalvecMLE(data1)
#' data2 = rnorm(100, mean  = .5, sd = 0.03)
#' fitnormalvecMLE(data2, logit = T)
fitnormalvecMLE <- function(data){
  if(!is.vector(data)){
    stop("data must be a vector")
  }
  mean = mean(data)
  var = (length(data) - 1) * var(data)/length(data) # var uses n-1 so we must multiply by (n-1)/n to get MLE
  return(list(mean = mean, var = var))
}



#' Beta Method of Moments
#' Calculates the method of monements estimators of a beta distribution fitted to the data
#' @param data - a vector of data 
#'
#' @return The fitted a,b of a beta(a,b) distribtion
#' @export
#'
#' @examples
#' set.seed(2)
#' data1 = rbeta(100, 1, 1)
#' fitbetavecMOM(data1)
#' data2 = rbeta(100, 43, 27)
#' fitbetavecMOM(data2)
fitbetavecMOM <-function(data){
  
  datamean =  mean(data)
  datavar = var(data)
  a = (datamean^2*(1-datamean)-datavar*datamean)/datavar
  b = a * (1-datamean)/datamean
  return(list(a = a, b = b))
}