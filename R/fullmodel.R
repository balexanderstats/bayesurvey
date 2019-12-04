#big model program

#' Title
#'
#' @param poll_data 
#' @param stateloc 
#' @param proploc 
#' @param varloc 
#' @param nloc 
#' @param invgamma 
#' @param VO 
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
#' polls = polls2016[, complete.cases(polls2016[, c(2, 3, 12)])]
#' head(polls)
#' nloc = which(colnames(polls2016) == "observations")
#' stateloc  = which(colnames(polls) == "State")
#' iterativeelectionmodelprop(polls, stateloc, c(2,3),  3, nloc = nloc, election_data = electdata)
iterativeelectionmodelprop = function(poll_data, stateloc, proploc, candidateloc,  varloc = NULL, nloc, invgamma = F, VO = NULL, a0 = NULL, b0 = NULL, election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue")){
  #step 1 get prior assignments
  newdf = propnormdfreplace(poll_data, proploc)
  priorout = addcategorytopolls(newdf, candidateloc,  stateloc, election_data, cutoffs, groupnames)
  finaldf = priorout$new_poll_data
  priormeans = priorout$priormean
  priorvars = priorout$priorvar
  statenames = unique(poll_data[ , stateloc])
  statenum = length(statenames)
  postmeans = rep(NA, statenum)
  postvars = rep(NA, statenum)
  for(i in 1:length(statenames)){
    poll_temp = subset(finaldf, finaldf[, stateloc] == statenames[i])
    groupnametemp = poll_temp$priorcat[1]
    groupnameloc = which(groupnames == groupnametemp, arr.ind = T)
    priormeantemp = priormeans[groupnameloc]
    priorvartemp = priorvars[groupnameloc]
    postout = unigausscpiterative(poll_temp[, candidateloc], priormeantemp, priorvartemp, datavar = NULL, poll_temp[, nloc], invgamma,  a0, b0)
    postmeans[i] = postout$finalpostmean
    postvars[i] = postout$finalpostvar
  }
  return(data.frame("State" = statenames, "Posterior Mean" = postmeans, "Posterior Variance" = postvars))
}