#this file will contain the normalization methods



#proportional normalization of data frames
#takes a data 
#' Proportional Normalization of data frame
#' 
#' Takes a data frame and normalizes it such that the rowSum is one
#' @param data 
#'
#' @return A data frame that is proportionally normalized
#' @export
#'
#' @examples
propnorm <- function(data){
  if(!is.data.frame(data) && !is.matrix(data)){
    stop("data is not matrix or data frame")
  }
  normdata = data / rowSums(data)
  
  return(normdata)
}

#' Title
#'
#' @param df data frame to be  
#' @param columns a vector of columns to be normalized
#'
#' @return a new data frame with the normalized data added with column names ending in normalized
#' @export
#'
#' @examples
propnormdf <- function(df, columns){
  data = df[columns]
  normdata = propnorm(data)
  colnames(normdata) = paste(colnames(data), "normalized")
  return(cbind(df,normdata))
}


  
  

