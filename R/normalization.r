#this file will contain the normalization methods



#' Proportional Normalization of data frame
#' 
#' Takes a data frame and normalizes it such that the rowSum is one
#' @param data the numeric data to be normailized
#'
#' @return A data frame that is proportionally normalized
#' @export
#'
#' @examples 
#' df1 = data.frame(a = c(0.4, 0.3, .35, .45), b = c(.45, .40, .40, .3))
#' propnorm(df1)
#' 
#' df2  = data.frame(a = c(0.3, .34, .38), b = c(0.2, 0.21, 0.25), c = c(.45, .42, .41))
#' propnorm(df2)
propnorm = function(data){
  if(!is.data.frame(data) && !is.matrix(data)){
    stop("data is not matrix or data frame")
  }
  normdata = data / rowSums(data)
  
  return(normdata)
}

#' Proportional Normalization of Data Frame
#'
#' @param df data frame to be  
#' @param columns a vector of columns to be normalized
#'
#' @return a new data frame with both original and the normalized data with column names ending in normalized
#' @export
#'
#' @examples
#' df1 = data.frame(a = c(0.4, 0.3, .35, .45), b = c(.45, .40, .40, .3))
#' propnormdf(df1)
#' 
#' df2  = data.frame(a = c(0.3, .34, .38), b = c(0.2, 0.21, 0.25), c = c(.45, .42, .41))
#' propnormdf(df2)
propnormdf = function(df, columns){
  data = df[columns]
  normdata = propnorm(data)
  colnames(normdata) = paste(colnames(data), "normalized")
  return(cbind(df,normdata))
}


#' Title
#' @param df data frame to be  
#' @param columns a vector of columns to be normalized
#'
#' @return a new data frame with both original and the normalized data with column names ending in normalized
#' @export
#'
#' @examples
#' df1 = data.frame(a = c(0.4, 0.3, .35, .45), b = c(.45, .40, .40, .3))
#' propnormdfreplace(df1)
#' 
#' df2  = data.frame(a = c(0.3, .34, .38), b = c(0.2, 0.21, 0.25), c = c(.45, .42, .41))
#' propnormdfreplace(df2)
propnormdfreplace = function(df, columns){
  normdata = propnorm(df[, columns])
  df[, columns] = normdata
  return(df)
}  
  

