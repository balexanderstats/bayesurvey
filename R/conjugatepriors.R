#conjugate prior functions

#' Title
#'
#' @param data a vector  of the data or a single point
#' @param priormean mean of prior distribution
#' @param priorvar variance of prior distribution
#' @param datavar datavariance,  if null var is called to estimate it
#' @param n sample size
#' @param invgamma a logical indicator of whether or not a inverse gamma prior for sigma is (T) or isn't included
#' @param a0 prior a if gamma is T
#' @param b0 prior b if gamma is T
#' @param V0 prior V0 where the likelihood is N(priormean)
#'
#' @return
#' @export
#'
#' @examples
unigausscp = function(data, priormean, priorvar, datavar = NULL,  n = NULL, invgamma = F, V0 = NULL, a0 = NULL, b0 = NULL){
  if(is.null(n)){
    n = length(data)
  }
  if(is.null(datavar)){
    datavar = var(data)
  }
  datamean =  mean(data)
  if(invgamma = F){
    #this does just a gaussian - gaussian conjugate prior assuming sigma is known or fixed
    datavar = var(data)
    postmean = priormean * (datavar / (n * priorvar + datavar)) + datamean * (n * priorvar)/(n * priorvar + datavar)
    postvar = (n/datavar + 1/priorvar)^-1
    posta= NULL
    postb = NULL
    dataweight = (n * priorvar)/(n * priorvar + datavar)
  
    }
  if(invgamma = T){
    #adds a prior distribution and posterior of sigma
    postvar = 1/(1/v0 + n)
    postmean = (priormean/priorvar + n * datamean)/postV
    posta = a0 + n/2
    postb = b0 + 0.5*(priormean^2/priorvar+sum(data^2)-postmean^2/postvar)
    dataweight = (n/postV)/(n/postV + (priorvar^-1)/postvar)
  }
  sdpost = sqrt(postvar)
  return(list(postmean = postmean, postvar = postvar, postsd = postsd, posta = posta,  postb = postb, dataweight = dataweight))
}


#' Iterative Gaussian Conjugate Prior
#' This function iteratives over a set of data.  The posterior becomes the new prior for the next data point
#'
#' @param data- data for analysis
#' @param priormean - intial prior mean
#' @param priorvar initial prior variance
#' @param datavar a vector of the variance for each data point
#' @param n a vector of sample sizes
#' @param invgamma a logical indicator of whether or not a inverse gamma prior for sigma is (T) or isn't included
#' @param V0 initial v0
#' @param a0 initial a0
#' @param b0 inital b0
#'
#' @return
#' @export
#'
#' @examples
unigausscpiterative = function(data, priormean, priorvar, datavar = NULL, n, invgamma = F, V0 = NULL, a0 = NULL, b0 = NULL){
  dataweight = rep(0, length(data))
  postmeans = rep(0, length(data))
  postvar = rep(0, length(data))
  postsd = rep(0, length(data))
  posta = rep(a0, length(data))
  postb = rep(b0, length(data))
  if(is.null(datavar)){
    datavar = (data*(1-data))/n
  }
  if(invgamma = F){
    for(i in 1:length(data)){
      out = unigausscp(data[i], priormean, priorvar, datavar = datavar[i], n[i])
      postmeans[i] = out$priormean
      postvar[i] = out$priorvar
      priorsd[i] = out$priorsd
      dataweight[i] = out$dataweight
      priormean = out$priormean
      priorvar = out$priorvar
      
    }
    
  }
  if(invgamma = T){
    for(u in 1:length(data)){
      out = unigausscp(data[i], priormean, priorvar, datavar, n, invgamma, V0, a0, b0)
      postmeans[i] = out$priormean
      postvar[i] = out$priorvar
      priorsd[i] = out$priorsd
      posta[i] = out$posta
      postb[i] = out$postb
      dataweight[i] = out$dataweight
      priormean = out$priormean
      priorvar = out$priorvar
      
      
    }
  }
}


