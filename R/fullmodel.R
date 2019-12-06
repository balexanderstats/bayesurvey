#big model program

#' Iterative Gaussian Conjugate Prior Model for American Election
#'
#' This function completely fits the Iterative Gaussian Proportional model in Alexander and Ellingson (2019) given raw poll data, election data, and basic parameters.
#' Given polling data on congressional districts or senate seats this could make similar predictions.
#' 
#' @param poll_data a data frame containing all necessary datat
#' @param stateloc the column number that contains the states for the poll
#' @param proploc the set of columns that contain poll results that need to be normalized
#' @param candidateloc the location of the candidate of interest
#' @param varloc the location of the variance if it is provided
#' @param nloc the location of the sample sizes
#' @param election_data the election data used to make the categorization decisions
#' @param cutoffs the cutoffs to assign the prior categories
#' @param groupnames  the names of the prior categories
#' @return a dataframe with the following rows: State with the state, Posterior Mean with the posterior mean, Posterior Variance with the posterior variance, and Posterior Standard Deviation with the posterior standard deviation.
#' @export

#' @examples
#' require(politicaldata)
#' elect2008  = subset(pres_results , year == 2008)
#' elect2008$margin = elect2008$dem - elect2008$rep
#' elect2012 = subset(pres_results , year == 2012)
#' elect2012$margin = elect2012$dem - elect2012$rep
#' electdata = data.frame("state" = elect2008$state, "2008" =  elect2008$margin, "2012" = elect2012$margin)
#' data(polls2016)
#' polls = polls2016[ complete.cases(polls2016[, c(2, 3, 12)]), ]
#' head(polls)
#' nloc = which(colnames(polls2016) == "observations")
#' stateloc  = which(colnames(polls) == "State")
#' iterativegaussianmodelprop(polls, stateloc, c(2,3),  3, nloc = nloc, election_data = electdata)
#' iterativegaussianmodelprop(polls, stateloc, c(2,3), 3, nloc = nloc, cutoffs = c(-.15, -.1, -0.05, 0.05, .1, .15), election_data = electdata)
iterativegaussianmodelprop = function(poll_data, stateloc, proploc, candidateloc,  varloc = NULL, nloc, election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue")){
  #step 1 get prior assignments
  
  #normalize data
  newdf = propnormdfreplace(poll_data, proploc)
  #call add category to polls to get prior information
  priorout = addcategorytopolls(newdf, candidateloc,  stateloc, election_data, cutoffs, groupnames)
  #saves new data frame with categories
  finaldf = priorout$new_poll_data
  #gets the names of state and the number of states
  statenames = unique(poll_data[ , stateloc])
  statenum = length(statenames)
  #saves prior mean and variance from the output
  priormeans = priorout$priormean
  priorvars =  priorout$priorvar

  #initializes posterior means
  postmeans = rep(NA, statenum)
  postvars = rep(NA, statenum)
  
  #loops through the states
  for(i in 1:statenum){
    #subsets the polls to get just one state
    poll_temp = subset(finaldf, finaldf[, stateloc] == statenames[i])
    
    #get the name of the group that state is in
    groupnametemp = poll_temp$priorcat[1]
    #get the location of that group in groupnames to use to find mean and variance
    groupnameloc = which(groupnames == groupnametemp, arr.ind = T)
    #gets priormean and variance for the state
    priormeantemp = priormeans[groupnameloc]
    priorvartemp = priorvars[groupnameloc]
    #calls iterative function
    postout = unigausscpiterative(poll_temp[, candidateloc], priormeantemp, priorvartemp, datavar = NULL, poll_temp[, nloc])
    #saves posterior
    postmeans[i] = postout$finalpostmean
    postvars[i] = postout$finalpostvar
  }
  postsd = sqrt(postvars)
  return(data.frame("State" = statenames, "Posterior Mean" = postmeans, "Posterior Variance" = postvars, "Posterior Standard Deviation" = postsd))
}

#' Gaussian Conjugate Prior Model for American Elections
#' This function completely fits the Iterative Gaussian Proportional model in Alexander and Ellingson (2019) given raw poll data, election data, and basic parameters.
#' Given polling data on congressional districts or senate seats this could make similar predictions.
#'
#' @param poll_data a data frame containing all necessary datat
#' @param stateloc the column number that contains the states for the poll
#' @param proploc the set of columns that contain poll results that need to be normalized
#' @param candidateloc the location of the candidate of interest
#' @param varloc the location of the variance if it is provided
#' @param nloc the location of the sample sizes
#' @param invgamma an indicator if an inverse gamma model is fit
#' @param v0 the hyperparameter for the shift parameter of a normal gamma distribution (mean, v0, a, b)
#' @param a0 the hyperparameter for the a paramater of a normal gamma distribution (mean, v0, a, b)
#' @param b0 the hyperparameter for the b parameter of a normal gamma distribution (mean, v0, a, b)
#' @param election_data the election data used to make the categorization decisions
#' @param cutoffs the cutoffs to assign the prior categories
#' @param groupnames  the names of the prior categories
#' @return a dataframe with the following rows: State with the state, Posterior Mean with the posterior mean, Posterior Variance with the posterior variance, and Posterior Standard Deviation with the posterior standard deviation.
#' @export

#' @examples
#' require(politicaldata)
#' elect2008  = subset(pres_results , year == 2008)
#' elect2008$margin = elect2008$dem - elect2008$rep
#' elect2012 = subset(pres_results , year == 2012)
#' elect2012$margin = elect2012$dem - elect2012$rep
#' electdata = data.frame("state" = elect2008$state, "2008" =  elect2008$margin, "2012" = elect2012$margin)
#' data(polls2016)
#' polls = polls2016[complete.cases(polls2016[, c(2, 3, 12)]), ]
#' head(polls)
#' nloc = which(colnames(polls2016) == "observations")
#' stateloc  = which(colnames(polls) == "State")
#' noniterativegaussianmodelprop(polls, stateloc, c(2,3),  3, nloc = nloc, election_data = electdata)
#' noniterativegaussianmodelprop(polls, stateloc, c(2,3), 3, nloc = nloc, invgamma = TRUE, v0 = 1, a0 = 0.0001, b0=0.0001, election_data = electdata)
noniterativegaussianmodelprop = function(poll_data, stateloc, proploc, candidateloc,  varloc = NULL, nloc, invgamma = F, v0 = NULL,  a0 = NULL, b0 = NULL, election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue")){
  #step 1 get prior assignments
  
  #normalizes the polling data
  newdf = propnormdfreplace(poll_data, proploc)
  #gets the list of state names
  statenames = unique(poll_data[ , stateloc])
  statenum = length(statenames)
  #calls function to add the prior category to polls
  priorout = addcategorytopolls(newdf, candidateloc,  stateloc, election_data, cutoffs, groupnames)
  
  #saves data frame with the categories
  finaldf = priorout$new_poll_data
  #gets the prior mean from the add categories output
  priormeans = priorout$priormean
  #sets prior variance to v0 for inverse gamma model and to the add category output for other functions
  if(invgamma == T){
    priorvars = rep(v0, length(groupnames))
  }
  if(invgamma == F){
    priorvars = priorout$priorvar
  }
  
  #initialize posterior mean and variance
  postmeans = rep(NA, statenum)
  postvars = rep(NA, statenum)
  
  for(i in 1:statenum){
    #subset the data to get only polls of one state
    poll_temp = subset(finaldf, finaldf[, stateloc] == statenames[i])
    #gets the groupname for the state
    groupnametemp = poll_temp$priorcat[1]
    #finds the location of the group name
    groupnameloc = which(groupnames == groupnametemp, arr.ind = T)
    #calls prior mean and variance for that group
    priormeantemp = priormeans[groupnameloc]
    priorvartemp = priorvars[groupnameloc]
    if(priormeantemp == 0 | is.na(priormeantemp)){
      print(i)
    }
    #runs model
    postout = unigausscp(poll_temp[, candidateloc], priormeantemp, priorvartemp, invgamma = invgamma, a0 = a0, b0 = b0)
    postmeans[i] = postout$postmean
    postvars[i] = postout$postvar
  }
  postsd = sqrt(postvars)
  return(data.frame("State" = statenames, "Posterior Mean" = postmeans, "Posterior Variance" = postvars, "Posterior Standard Deviation" = postsd))
}
