#conjugate prior functions

#' Title
#'
#' @param data a vector of the data
#' @param priormean mean of prior distribution
#' @param priorvar variance of prior distribution
#' @param n sample size
#' @param invgamma a logical indicator of whether or not a gamma prior for sigma is (T) or isn't included
#' @param a0 prior a if gamma is T
#' @param b0 prior b if gamma is T
#' @param V0 prior V0 where the likelihood is N(priormean)
#'
#' @return
#' @export
#'
#' @examples
unigausscp = function(data, priormean, priorvar, n = NULL, invgamma = F, V0 = NULL, a0 = NULL, b0 = NULL){
  if(is.null(n)){
    n = length(data)
  }
  datamean =  mean(data)
  if(invgamma = F){
    #this does just a gaussian - gaussian conjugate prior assuming sigma is known or fixed
    datavar = var(data)
    postmean = priormean * (datavar / (n * priorvar + datavar)) + datamean * (n * priorvar)/(n * priorvar + datavar)
    postvar = (n/datavar + 1/priorvar)^-1
    apost = NULL
    bpost = NULL
    dataweight = (n * priorvar)/(n * priorvar + datavar)
  
    }
  if(invgamma = T){
    #adds a prior distribution and posterior of sigma
    postvar = 1/(1/v0 + n)
    postmean = (priormean/priorvar + n * datamean)/postV
    apost = a0 + n/2
    bpost = b0 + 0.5*(priormean^2/priorvar+sum(data^2)-postmean^2/postvar)
    dataweight = (n/postV)/(n/postV + (priorvar^-1)/postvar)
  }
  sdpost = sqrt(postvar)
  return(list(postmean = postmean, postvar = postvar, postsd = postsd, apost = apost,  bpost = bpost, dataweight = dataweight))
}


unigausscpiterative = function(data, priormean, priorvar, n, invgamma = F, V0 = NULL, a0 = NULL, b0 = NULL){
  dataweight = rep(0, length(data))
  postmeans = rep(0, length(data))
  postvar = rep(0, length(data))
  postsd = rep(0, length(data))
  if(invgamma = F){
    for(i in 1:length(data)){
      out = unigausscp(data[i], priormean, priorvar, n[i])
      postmeans[i] = out$priormean
      postvar[i] = out$priorvar
      priorsd[i] = out$priorsd
      dataweight[i] = out$dataweight
      priormean = out$priormean
      priorvar = out$priorvar
      
    }
  }
  if(invgamma = T){
    #complete later
    return(1)
  }
}


