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
#' @param npolls the number of polls to include, default is all polls
#' @param dateloc the location of the date used to determine the last npolls
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
iterativegaussianmodelprop = function(poll_data, stateloc, proploc, candidateloc,  varloc = NULL, nloc,  npolls = NULL, dateloc = NULL, election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue")){
  #step 1 get prior assignments
  #gets the names of state and the number of states
  statenames = unique(election_data[ , 1])
  statenum = length(statenames)
  if(!is.null(varloc)){
    datavar = poll_data[, varloc]}
  if(is.null(varloc)){
    datavar = NULL
  }
  #normalize data
  newdf = propnormdfreplace(poll_data, proploc)
  if(!is.null(npolls)){
    lastpolls = numeric(0)
    for(i in 1:length(statenames)){
      #Step 2: For each state get index of last n polls and add to array of last polls
      poll_temp = subset(newdf, newdf[, stateloc] == statenames[i])
      poll_temp = poll_temp[order(poll_temp[, dateloc], decreasing = T), ]
      if(nrow(poll_temp) > npolls){
        lastpolls = c(lastpolls, poll_temp$X[1:npolls])
      }
      else{
        lastpolls = c(lastpolls, poll_temp$X)
      }
      
      #Step 3: redefine poll data
      #go through original 
      
    }
    newdf = newdf[newdf$X %in% lastpolls, ]
  }
  
  #call add category to polls to get prior information
  priorout = addcategorytopolls(newdf, candidateloc,  stateloc, election_data, cutoffs, groupnames)
  #saves new data frame with categories
  finaldf = priorout$new_poll_data

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
    if(nrow(poll_temp) == 0){
      priorassign = getpriorassign(election_data = election_data, cutoffs = cutoffs, groupnames = groupnames)
      groupnametemp = priorassign[priorassign[, 1] == statenames[i] , 2] 
      groupnameloc = which(groupnames == groupnametemp, arr.ind = T)
      postmeans[i] = priormeans[groupnameloc]
      postvars[i] = priorvars[groupnameloc]
    }
    
    
    #calls iterative function
    else{
      #get the name of the group that state is in
      groupnametemp = poll_temp$priorcat[1]
      #get the location of that group in groupnames to use to find mean and variance
      groupnameloc = which(groupnames == groupnametemp, arr.ind = T)
      #gets priormean and variance for the state
      priormeantemp = priormeans[groupnameloc]
      priorvartemp = priorvars[groupnameloc]
      postout = unigausscpiterative(poll_temp[, candidateloc], priormeantemp, priorvartemp, datavar = datavar, poll_temp[, nloc])
      #saves posterior
      postmeans[i] = postout$finalpostmean
      postvars[i] = postout$finalpostvar
    }
    
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
#' @param npolls the number of polls to include, default is all polls
#' @param dateloc the location of the date used to determine the last npolls
#' @param idloc the location of an id variable (necessary when npoll isn't null)
#' @param invgamma an indicator if an inverse gamma model is fit
#' @param v0 the hyperparameter for the shift parameter of a normal gamma distribution (mean, v0, a, b)
#' @param a0 the hyperparameter for the a paramater of a normal gamma distribution (mean, v0, a, b)
#' @param b0 the hyperparameter for the b parameter of a normal gamma distribution (mean, v0, a, b)
#' @param logit if a logit transformation is done.
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
#' dateloc =  which(colnames(polls) == "end_date")
#' noniterativegaussianmodelprop(polls, stateloc, c(2,3),  3, nloc = nloc, election_data = electdata)
#' noniterativegaussianmodelprop(polls, stateloc, c(2,3),  3, dateloc = 9, nloc = nloc, election_data = electdata, npolls = 10)
#' noniterativegaussianmodelprop(polls, stateloc, c(2,3), 3, nloc = nloc, invgamma = TRUE, v0 = 1, a0 = 0.0001, b0=0.0001, election_data = electdata)
#' noniterativegaussianmodelprop(polls, stateloc, c(2,3), 3, nloc = nloc, invgamma = TRUE, v0 = 1, a0 = 0.0001, b0=0.0001, logit = T, election_data = electdata)


noniterativegaussianmodelprop = function(poll_data, stateloc, proploc, candidateloc,  varloc = NULL, nloc, npolls = NULL, dateloc = NULL, idloc = NULL, invgamma = F, v0 = NULL,  a0 = NULL, b0 = NULL, logit = F, election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue")){
  #step 1 get prior assignments
  
  #normalizes the polling data
  newdf = propnormdfreplace(poll_data, proploc)
  #gets the list of state names
  statenames = unique(election_data[ , 1])
  statenum = length(statenames)

  
  #implementing last n polls
  #Step 1: Loop over states
  if(!is.null(npolls)){
    lastpolls = numeric(0)
    for(i in 1:length(statenames)){
      #Step 2: For each state get index of last n polls and add to array of last polls
      poll_temp = subset(newdf, newdf[, stateloc] == statenames[i])
      poll_temp = poll_temp[order(poll_temp[, dateloc], decreasing = T), ]
      if(nrow(poll_temp) > npolls){
        lastpolls = c(lastpolls, poll_temp[1:npolls, idloc])
      }
      else{
        lastpolls = c(lastpolls, poll_temp[, idloc])
      }
      
      #Step 3: redefine poll data
      #go through original 
      
    }
    newdf = newdf[newdf[, idloc] %in% lastpolls, ]
  }
  #calls function to add the prior category to polls
  priorout = addcategorytopolls(newdf, candidateloc,  stateloc, election_data, cutoffs, groupnames, logit = logit)
  
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
  cilow = rep(NA, statenum)
  cihigh = rep(NA, statenum)
  winprob = rep(NA, statenum)
  
  for(i in 1:statenum){
    #subset the data to get only polls of one state
    poll_temp = subset(finaldf, finaldf[, stateloc] == statenames[i])
    if(nrow(poll_temp) <  2){
      priorassign = getpriorassign(election_data = election_data, cutoffs = cutoffs, groupnames = groupnames)
      groupnametemp = priorassign[priorassign[, 1] == statenames[i] , 2] 
      groupnameloc = which(groupnames == groupnametemp, arr.ind = T)
      if(logit == F){
        postmeans[i] = priormeans[groupnameloc]
        postvars[i] = priorout$priorvar[groupnameloc]
        cilow[i] = postmeans[i] - 1.96*sqrt(postvars[i])
        cihigh[i] = postmeans[i] + 1.96*sqrt(postvars[i])
        winprob[i] = pnorm(.5, mean  = postmeans[i], sd = sqrt(postvars[i]))
      }
      if(logit == T){
        polls_tempgroup = finaldf[finaldf$priorcat == groupnametemp, ]
        postmeans[i] = mean(polls_tempgroup[ ,candidateloc])
        postvars[i] = var(polls_tempgroup[, candidateloc])
        cilow[i] = postmeans[i] - 1.96*sqrt(postvars[i])
        cihigh[i] = postmeans[i] + 1.96*sqrt(postvars[i])
        winprob[i] = pnorm(.5, mean  = postmeans[i], sd = sqrt(postvars[i]))
        
      }
      
    }
    
    else{
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
    postout = unigausscp(poll_temp[, candidateloc], priormeantemp, priorvartemp,logit = logit, invgamma = invgamma, a0 = a0, b0 = b0)
    postmeans[i] = postout$postmean
    postvars[i] = postout$postvar
    cilow[i] = postout$cilow
    cihigh[i] = postout$cihigh
    winprob[i] = postout$winprob
    }
    #
    ciwidth = (cihigh -cilow)/2
  }

  postsd = sqrt(postvars)
  return(data.frame("State" = statenames, "Posterior Mean" = postmeans, "Posterior Variance" = postvars, "Posterior Standard Deviation" = postsd, "Credible Interval Low" = cilow, "Credible Interval High" = cihigh, "Margin of Error" = ciwidth, "Win Probability" = winprob))
}




#' Gaussian Rolling Average
#' @param poll_data 
#' @param stateloc 
#' @param proploc 
#' @param candidateloc 
#' @param npolls 
#' @param varloc 
#' @param nloc 
#' @param invgamma 
#' @param v0 
#' @param a0 
#' @param b0 
#' @param election_data 
#' @param cutoffs 
#' @param groupnames 
#'
#' @return
#' @export
#'
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
#' gaussianrollingaverage(polls, stateloc = stateloc, proploc = c(2,3), 3, nloc = nloc, election_data = electdata)
gaussianrollingaverage = function(poll_data, stateloc, proploc, candidateloc, npolls = 10,  varloc = NULL, nloc, invgamma = F, v0 = NULL,  a0 = NULL, b0 = NULL, election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue")){
  
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
    if(nrow(poll_temp) > npolls){
      poll_temp = poll_temp[order(poll_temp$end_date, decreasing = T),]
      poll_temp = poll_temp[1:npolls,]
    }
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


#' Rolling Average Iterative Gaussian Bayesian Model
#'
#' @param poll_data 
#' @param stateloc 
#' @param proploc 
#' @param candidateloc 
#' @param npolls 
#' @param varloc 
#' @param nloc 
#' @param invgamma 
#' @param v0 
#' @param a0 
#' @param b0 
#' @param election_data 
#' @param cutoffs 
#' @param groupnames 
#'
#' @return
#' @export
#'
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
#' gaussianrollingaverage(polls, stateloc = stateloc, proploc = c(2,3), 3, nloc = nloc, election_data = electdata)
gaussianrollingaverageiterative = function(poll_data, stateloc, proploc, candidateloc, npolls = 10,  varloc = NULL, nloc, invgamma = F, v0 = NULL,  a0 = NULL, b0 = NULL, election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue")){
  
  #normalizes the polling data
  newdf = propnormdfreplace(poll_data, proploc)
  #gets the list of state names
  statenames = unique(poll_data[ , stateloc])
  statenum = length(statenames)
  #calls function to add the prior category to polls
  priorout = addcategorytopolls(newdf, candidateloc,  stateloc, election_data, cutoffs, groupnames)
  if(!is.null(varloc)){
    datavar = poll_data[, varloc]}
  if(is.null(varloc)){
    datavar = NULL
  }
  
  
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
    if(nrow(poll_temp) > npolls){
      poll_temp = poll_temp[order(poll_temp$end_date, decreasing = T), ]
      poll_temp = poll_temp[1:npolls, ]
    }
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
    postout = unigausscpiterative(poll_temp[, candidateloc], priormeantemp, priorvartemp, datavar = datavar, n = nloc)
    postmeans[i] = postout$finalpostmean
    postvars[i] = postout$finalpostvar
  }
  postsd = sqrt(postvars)
  return(data.frame("State" = statenames, "Posterior Mean" = postmeans, "Posterior Variance" = postvars, "Posterior Standard Deviation" = postsd))
}




#' Title
#'
#' @param poll_data 
#' @param stateloc 
#' @param proploc 
#' @param candidateloc 
#' @param npolls 
#'
#' @return
#' @export
#'
#' @examples
#' require(politicaldata)
#' elect2008  = subset(pres_results , year == 2008)
#' elect2008$margin = elect2008$dem - elect2008$rep
#' elect2012 = subset(pres_results , year == 2012)
#' elect2012$margin = elect2012$dem - elect2012$rep
#' electdata = data.frame("state" = elect2008$state, "2008" =  elect2008$margin, "2012" = elect2012$margin)
#' data(polls2016)
#' polls = polls2016[ complete.cases(polls2016[, c(2, 3, 12)]), ]
#' stateloc  = which(colnames(polls) == "State")
#' pollrollingaverage(polls, stateloc = stateloc, proploc = c(2,3), 3)
pollrollingaverage=function(poll_data, stateloc, proploc, candidateloc, npolls = 10){
  #normalizes the polling data
  newdf = propnormdfreplace(poll_data, proploc)
  #gets the list of state names
  statenames = unique(poll_data[ , stateloc])
  statenum = length(statenames)
  
  #initialize state mean and variance
  statemeans = rep(NA, statenum)
  statevars = rep(NA, statenum)
  
  for(i in 1:statenum){
    #subset the data to get only polls of one state
    poll_temp = subset(newdf, newdf[, stateloc] == statenames[i])
    if(nrow(poll_temp) > npolls){
      poll_temp = poll_temp[order(poll_temp$end_date, decreasing = T),]
      poll_temp = poll_temp[1:npolls,]
      statemeans[i] = mean(poll_temp[, candidateloc])
      statevars[i] = var(poll_temp[, candidateloc])
    }
    else{
      statemeans[i] = mean(poll_temp[, candidateloc])
      statevars[i] = var(poll_temp[, candidateloc])
      
    }
    statesd = sqrt(statevars)
    
  }
  return(data.frame("State" = statenames, "Predicted Mean" = statemeans, "Predicted Variance" = statevars, "Predicted Standard Deviation" = sqrt(statevars)))
}




#margin based models
#prior function is the same
#have input for normalize vs not normalize
# no iterative model
# rolling average


#' Title
#'
#' @param poll_data 
#' @param stateloc 
#' @param marginloc 
#' @param npolls the number of polls to include, default is all polls
#' @param dateloc the location of the date used to determine the last npolls
#' @param normalize whether or not to normalize poll date before margin is calculated
#' @param invgamma 
#' @param v0 
#' @param a0 
#' @param b0 
#' @param election_data 
#' @param cutoffs 
#' @param groupnames 
#'
#' @return
#' @export
#'
#' @examples
#' require(politicaldata)
#' elect2008  = subset(pres_results , year == 2008)
#' elect2008$margin = elect2008$dem - elect2008$rep
#' elect2012 = subset(pres_results , year == 2012)
#' elect2012$margin = elect2012$dem - elect2012$rep
#' electdata = data.frame("state" = elect2008$state, "2008" =  elect2008$margin, "2012" = elect2012$margin)
#' data(polls2016)
#' polls = polls2016[ complete.cases(polls2016[, c(2, 3, 12)]), ]
#' stateloc  = which(colnames(polls) == "State")
#' polls2016$margin = polls2016$Trump - polls2016$Clinton
#' marginloc = which(colnames(polls2016) == "margin")
#' margingaussianmodelprop(polls2016, stateloc = stateloc, marginloc = marginloc, election_data = electdata)
margingaussianmodelprop = function(poll_data, stateloc, marginloc, nloc, proploc = NULL,  npolls = NULL, dateloc = NULL, normalize = F, invgamma = F, v0 = NULL,  a0 = NULL, b0 = NULL, election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue")){
  #step 1 get prior assignments
  
  #normalizes the polling data
  
  if(normalize == T){
    newdf = propnormdfreplace(poll_data, proploc)
    newdf[, marginloc ] = 100*(newdf[, proploc[1]] - newdf[, proploc[2]])}
  if(normalize == F){
    newdf = poll_data}
  #gets the list of state names
  statenames = unique(election_data[ , 1])
  statenum = length(statenames)
  if(!is.null(npolls)){
    lastpolls = numeric(0)
    for(i in 1:length(statenames)){
      #Step 2: For each state get index of last n polls and add to array of last polls
      poll_temp = subset(newdf, newdf[, stateloc] == statenames[i])
      poll_temp = poll_temp[order(poll_temp[, dateloc], decreasing = T), ]
      if(nrow(poll_temp) > npolls){
        lastpolls = c(lastpolls, poll_temp$X[1:npolls])
      }
      else{
        lastpolls = c(lastpolls, poll_temp$X)
      }
      
      #Step 3: redefine poll data
      #go through original 
      
    }
    newdf = newdf[newdf$X %in% lastpolls, ]
  }
  
  #calls function to add the prior category to polls
  priorout = addcategorytopollsmargin(newdf, marginloc,  stateloc, election_data, cutoffs, groupnames)
  
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
    if(nrow(poll_temp) < 2){
      priorassign = getpriorassign(election_data = election_data, cutoffs = cutoffs, groupnames = groupnames)
      groupnametemp = priorassign[priorassign[, 1] == statenames[i] , 2] 
      groupnameloc = which(groupnames == groupnametemp, arr.ind = T)
      postmeans[i] = priormeans[groupnameloc]
      postvars[i] = priorvars[groupnameloc]
    }
    else{
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
    if(nrow(poll_temp) == 1){
      n = poll_temp[, nloc]
    }
    if(nrow(poll_temp) > 1){
      n = NULL
    }
    #runs model
    postout = unigausscp(poll_temp[, marginloc], priormeantemp, priorvartemp, n = n,  invgamma = invgamma, a0 = a0, b0 = b0)
    postmeans[i] = postout$postmean
    postvars[i] = postout$postvar
    }
  }
  postsd = sqrt(postvars)
  return(data.frame("State" = statenames, "Posterior Mean" = postmeans, "Posterior Variance" = postvars, "Posterior Standard Deviation" = postsd))
}
