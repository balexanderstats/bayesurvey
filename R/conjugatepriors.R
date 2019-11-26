#conjugate prior functions

#' Title
#'
#' @param data a vector  of the data or a single point
#' @param priormean mean of prior distribution
#' @param priorvar variance of prior distribution if invgamma is not true, shift parameter for variance if invgamma is true
#' @param datavar data variance,  if null var is called to estimate it
#' @param singlepoll if true it calculates the variance based on the sample size
#' @param n approriate sample size
#' @param invgamma a logical indicator of whether or not a inverse gamma prior for sigma is (T) or isn't included
#' @param a0 prior a if invgamma is T
#' @param b0 prior b if invgamma is T
#'
#' @return A list with the following components:  priormean: the prior mean, priorvar: the prior var, a0: the prior of the inverse gamma distribution if invgamma is true, b0: the prior of the inverse gamma distribution if invgamma is true, n: which is the sample size if supplied,  datavar: which is the variance of the data if it is supplied, invgamma: which states if the inverse gamma prior was used, postmean: which is the posterior mean, postvar: which is the posterior variance, postsd: which is the posterior standard deviation, posta: which is the posterior a of a inverse gamma (a,b), and postb: which is the posterior of b of a inverse gamma(a, b), and dataweight: a scalar of the weight of the data over the prior. 
#' @export
#'
#' @examples
#' set.seed(10)
#' data1 = rnorm(30, mean = .48, sd = 0.05)
#' n1 = floor(runif(1, 200, 800))
#' datavar1 = mean(sqrt(data1[1]*(1-data1[1])/n1))
#' unigausscp(data1, 0.5, 0.05)
#' unigausscp(data1[1], 0.5, 0.05, singlepoll = T, n = n1)
#' unigausscp(data1, 0.5, 0.05, invgamma = T, a0 = 0.0001, b0 = 0.0001)
#' unigausscp(data1[1], 0.5, 0.05, singlepoll = T, n = n1, invgamma = T, a0 = 0.0001, b0 = 0.0001)
#' 
#' 
unigausscp = function(data, priormean, priorvar, datavar = NULL, singlepoll = F,  n = NULL, invgamma = F, a0 = NULL, b0 = NULL){
  if(is.null(n)){
    n = length(data)
  }
  if(is.null(datavar)){
    if(singlepoll == T){
      datavar = sqrt(data*(1-data)/n)
    }
    if(singlepoll == F){
      datavar = var(data)
    }
  }
  datamean =  mean(data)
  if(invgamma == F){
    #this does just a gaussian - gaussian conjugate prior assuming sigma is known or fixed
    postmean = priormean * (datavar / (n * priorvar + datavar)) + datamean * (n * priorvar)/(n * priorvar + datavar)
    postvar = (n/datavar + 1/priorvar)^-1
    posta= NULL
    postb = NULL
    dataweight = (n * priorvar)/(n * priorvar + datavar)
  
    }
  if(invgamma == T){
    if(singlepoll == T){
      ssqdata = round(datamean * n, 0)
    }
    #adds a prior distribution and posterior of sigma
    if(singlepoll == F){
      ssqdata = sum(data^2)
      }
    postv = 1/(1/priorvar + n)
    postmean = (priormean/priorvar + n * datamean) * postv
    posta = a0 + n/2
    postb = b0 + 0.5 * (priormean^2 /priorvar + ssqdata - postmean^2/postv)
    postvar = postb / ((posta - 1)*postv)
    dataweight = (n/postv)/(n/postv + (priorvar^-1)/postvar)
  }
  postsd = sqrt(postvar)
  return(list(priormean  = priormean,  priorvar = priorvar, a0 = a0, b0 = b0,  n = n, datavar = datavar, invgamma = invgamma, postmean = postmean, postvar = postvar, postsd = postsd, posta = posta,  postb = postb, dataweight = dataweight))
}


#' Iterative Gaussian Conjugate Prior
#' It iterates over the data and calls the function unigausscp for each point, setting the new prior as the posterior from the last data point.
#'
#' @param data data for analysis
#' @param priormean - intial prior mean
#' @param priorvar initial prior variance
#' @param datavar a vector of the variance for each data point
#' @param n a vector of sample sizes
#' @param invgamma a logical indicator of whether or not a inverse gamma prior for sigma is (T) or isn't included
#' @param a0 initial a0
#' @param b0 inital b0
#'
#' @return A list with the following elements: finalpostmean: a scalar of the final mean of the posterior, finalpostvar: the final posterior variance, finalpostsd: the final posterior standard deviation, dataweights: a vector of the weight of the data in the posterior, postmeans: a vector of the posterior means for each data point, postvar:  a vector of the posterior variance for each data point, postsds: a vector of the posterior standard deviations for each data point, posta: a vector of the posterior values of a if invgamma = T, and postb:  a vector of the posterior values of b in invgamma = T.
#' @export
#'
#' @examples
#' set.seed(12)
#' data1 = c(rnorm(10, mean = .48, sd = 0.05), rnorm(20, mean = .49, sd = 0.05),  rnorm(20, mean = .51, sd = 0.05))
#' n = floor(runif(50, 200, 800))
#' unigausscpiterative(data1, 0.5, 0.05, n = n) 
#' unigausscpiterative(data1, 0.5, 0.05, n = n, invgamma = T, a0 = 0.0001, b0 = 0.0001)
unigausscpiterative = function(data, priormean, priorvar, datavar = NULL, n, invgamma = F,  a0 = NULL, b0 = NULL){
  dataweights = rep(0, length(data))
  postmeans = rep(0, length(data))
  postvars = rep(0, length(data))
  postsds = rep(0, length(data))
  posta = rep(a0, length(data))
  postb = rep(b0, length(data))
  if(is.null(datavar)){
    datavar = (data*(1-data))/n
  }
  if(invgamma == F){
    for(i in 1:length(data)){
      out = unigausscp(data[i], priormean, priorvar, datavar = datavar[i], n[i])
      postmeans[i] = out$postmean
      postvars[i] = out$postvar
      postsds[i] = out$postsd
      dataweights[i] = out$dataweight
      priormean = out$postmean
      priorvar = out$postvar
      
    }
    
  }
  if(invgamma == T){
    for(i in 1:length(data)){
      out = unigausscp(data[i], priormean, priorvar, datavar = datavar[i], singlepoll = T, n = n[i], invgamma, a0, b0)
      postmeans[i] = out$postmean
      postvars[i] = out$postvar
      postsds[i] = out$postsd
      posta[i] = out$posta
      postb[i] = out$postb
      dataweights[i] = out$dataweight
      priormean = out$postmean
      priorvar = out$postvar
      a0 = out$posta
      b0 = out$postb
      
      
      
    }

  }
  finalpostmean = priormean
  finalpostvar = priorvar
  finalpostsd = sqrt(priorvar)
  return(list(finalpostmean = finalpostmean,  finalpostvar = finalpostvar, finalpostsd = finalpostsd, dataweights = dataweights, postmeans = postmeans, postvars = postvars, postsds = postsds, posta = posta, postb = postb))
}


