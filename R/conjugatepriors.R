#conjugate prior functions

#' Univariate Iterative Gaussian Conjugate Prior
#' This function takes a vector of proportions and iteratives over it applying a gaussian conjugate prior calculation at each point, and using the posterior as the new prior in the next iteration.
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
#' unigausscp(data1, 0.5, 1, invgamma = T, a0 = 0.0001, b0 = 0.0001)
#' 
#' 
unigausscp = function(data, priormean, priorvar, datavar = NULL, singlepoll = F,  n = NULL, invgamma = F, a0 = NULL, b0 = NULL){
  #checks that n is the number of polls if singlepoll = F
  if(!is.null(n)){
    if(singlepoll == F  & n != length(data)){
      n = length(data)
    }
  }
  
  #checks if data is a vector
  if(!is.vector(data)){
    stop("Data is not a vector")
  }
  
  #generates n if not provided
  if(is.null(n)){
    n = length(data)
  }
  #gets datavar for a single poll or multiple polls if datavar is not provided
  if(is.null(datavar)){
    if(singlepoll == T){
      datavar = sqrt(data*(1-data)/n)
    }
    if(singlepoll == F){
      datavar = var(data)
    }
  }
  #calcualtes data mean
  datamean =  mean(data)
  
  #No inverse gamma prior
  if(invgamma == F){
    #this does just a gaussian - gaussian conjugate prior assuming sigma is known or fixed
    #calculates posterior mean
    postmean = priormean * (datavar / (n * priorvar + datavar)) + datamean * (n * priorvar)/(n * priorvar + datavar)
    #calculates posterior variance
    postvar = (n/datavar + 1/priorvar)^-1
    #returns null for the parameters relating to the invgamma distribution
    posta = NULL
    postb = NULL
    postv = NULL
    #calculates the percent that the data plays in the weighted average of the mean
    dataweight = (n * priorvar)/(n * priorvar + datavar)
    
  }
  # Includes inverse gamma prior
  #adds a prior distribution and posterior of sigma
  
  if(invgamma == T){
    if(is.null(a0) | is.null(b0)){
      stop("a0 and b0 must be specified for inverse gamma model")
    }
    if(a0 < 0 | b0 < 0){
      stop("a0 and b0 must be positive")
    }
    #gets sum(x_i^2) where x_i is each respondent of the polls
    if(singlepoll == T){
      ssqdata = round(datamean * n, 0)
    }

    if(singlepoll == F){
      #gets sum(x_i^2) where x_i is the result of a poll
      ssqdata = sum(data^2)
    }
    #calculates posterior for v
    postv = 1/(1/priorvar + n)
    #calculates posterior for the mean
    postmean = (priormean/priorvar + n * datamean) * postv
    #calculates posterior a and b
    posta = a0 + n/2
    postb = b0 + 0.5 * (priormean^2 /priorvar + ssqdata - postmean^2/postv)
    #calculates the posterior variance
    postvar = postb / ((posta - 1)*postv)
    #calculates the percent that the data plays in the weighted average of the mean
    dataweight = (n/postv)/(n/postv + (priorvar^-1)/postvar)
  }
  postsd = sqrt(postvar)
  return(list(priormean  = priormean,  priorvar = priorvar, postv = postv, a0 = a0, b0 = b0,  n = n, datavar = datavar, invgamma = invgamma, postmean = postmean, postvar = postvar, postsd = postsd, posta = posta,  postb = postb, dataweight = dataweight))
}


#' Iterative Gaussian Conjugate Prior
#' It iterates over the data and calls the function unigausscp for each point, setting the new prior as the posterior from the last data point.
#' 
#' @param data data for analysis
#' @param priormean - intial prior mean
#' @param priorvar initial prior variance, if inverse gamma is desired this is the initial V0 parameter of the 
#' @param datavar a vector of the variance for each data point
#' @param n a vector of sample sizes

#'
#' @return A list with the following elements: finalpostmean: a scalar of the final mean of the posterior, finalpostvar: the final posterior variance, finalpostsd: the final posterior standard deviation, dataweights: a vector of the weight of the data in the posterior, postmeans: a vector of the posterior means for each data point, postvar:  a vector of the posterior variance for each data point, postsds: a vector of the posterior standard deviations for each data point, posta: a vector of the posterior values of a if invgamma = T, and postb:  a vector of the posterior values of b in invgamma = T.
#' @export
#'
#' @examples
#' set.seed(12)
#' data1 = c(rnorm(10, mean = .48, sd = 0.05), rnorm(20, mean = .49, sd = 0.05),  rnorm(20, mean = .51, sd = 0.05))
#' n = floor(runif(50, 200, 800))
#' unigausscpiterative(data1, 0.5, 0.05, n = n)
#' data2 = c(rnorm(50, mean = 0.5, 0.05))
#' unigausscpiterative(data2, 0.5, 0.05, n = n)

unigausscpiterative = function(data, priormean, priorvar, datavar = NULL, n){
  #initialize the vector that stores the dataweight
  dataweights = rep(0, length(data))
  #intialize the vector that stores the posterior mean
  postmeans = rep(0, length(data))
  #initalize the vector that stores the posterior variance
  postvars = rep(0, length(data))
  #initialize the vector that stores the posterior standard deviation
  postsds = rep(0, length(data))
  #initalize the vectors for the parameters v, a,b, of a normal-gamma(mean, v, a, b) distribution 
  #calculates the variance of the poll if it is not provided
  if(is.null(datavar)){
    datavar = (data*(1-data))/n
  }
  
  if(sum(data > 1) > 1 | sum(data < 0) > 1){
    stop("The data contains numbers outside of 0 and 1.")
  }
  #iterates over the polls and uses the posterior as the prior in the next iteration
  for(i in 1:length(data)){
    #call function
    out = unigausscp(data[i], priormean, priorvar, datavar = datavar[i], n[i])
    #store posterior mean, variance, standard deviation, the data weight into their vectors
    postmeans[i] = out$postmean
    postvars[i] = out$postvar
    postsds[i] = out$postsd
    dataweights[i] = out$dataweight
    #resets priormean and priorvar to be the posterior values of this loop
    priormean = out$postmean
    priorvar = out$postvar
    
  }

  finalpostmean = out$postmean
  finalpostvar = out$postvar
  finalpostsd = sqrt(finalpostvar)
  return(list(finalpostmean = finalpostmean,  finalpostvar = finalpostvar, finalpostsd = finalpostsd, dataweights = dataweights, postmeans = postmeans, postvars = postvars, postsds = postsds))
}


