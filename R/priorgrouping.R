#data splitting



#' Title
#'
#' @param election_data the election data in df form. The first column must contain the name of that state
#' @param cutoffs the cutoffs used to split the data into the categories
#' @param groupnames 
#' @param weights optional weights for a weighted average of the columns
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
#' data = data.frame(state = elect2008$state, 2008 =  elect2008$margin, 2012 = elect2012$margin)
#' getpriorassign(data)
#' weight = c(0.25,0.75)
#' getpriorassign(data , weights = weight)
getpriorassign = function(election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue"), weights =  NULL){
  statenames = as.character(election_data[,1])
  if(!is.null(weights)){
    weightmat = matrix(weights, ncol = length(weights), nrow = nrow(data))
    weighteddata = weightmat * election_data[, -1]
    avgdata = rowSums(weighteddata)/(ncol(election_data)-1)
  }
  else{
    avgdata = rowMeans(election_data[, -1])
  }
  assignments = rep("", nrow(election_data))
  assignments[which(avgdata < cutoffs[1], arr.ind = T)] = groupnames[1]
  categories = length(groupnames)
  assignments[which(avgdata > cutoffs[categories - 1])] = groupnames[categories]
  if(categories > 2){
   for(i in 2:(categories-1)){
     assignments[which((avgdata > cutoffs[i-1]) & (avgdata < cutoffs[i]), arr.ind = T)] = groupnames[i]
   } 
  }
  
  return(cbind(statenames, assignments))
}



#' Adds Category assignment to poll data.
#' This function takes the poll data and adds the assigned prior region to the end of the data frame.  It also returns the mean of the prior distributions and the variance of the prior distributions, along with the vector appended to the data frame.
#' @param poll_data the data to fit the prior distribution , the first column must be the state name
#' @param proploc the column the adjusted poll proportion is located
#' @param stateloc
#' @param election_data 
#' @param cutoffs 
#' @param groupnames 
#' @param yearweight weights the columns differently 
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
#' set.seed(16)
#' poll1 = polls2016[sample(1:nrow(polls2016), 500) ,]
#' poll2 = polls2016[sample(1:nrow(polls2016), 500), ]
#' propnormpolls1 = propnormdf(poll1, c(2,3))
#' propnormpolls2 = propnormdf(poll2, c(2,3))
#' addcategorytopolls(propnormpolls1, 30, 18, electdata)
#' addcategorytopolls(propnormpolls2, 30, 18, electdata, cutoffs = c(-.15, -.1, -0.05, 0.05, .1, .15))

addcategorytopolls = function(poll_data, proploc, stateloc, election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue"), yearweight = NULL){
  assignment = getpriorassign(election_data, cutoffs, groupnames)
  categories  = length(groupnames)
  priormean = rep(NA, categories)
  priorvar = rep(NA, categories)
  priorcat = rep("", nrow(poll_data))
  for(i in 1:categories){
    statenames = assignment[which(assignment[ , 2] == groupnames[i]) , 1]
    priorcat[which(poll_data[, stateloc] %in% statenames)] = groupnames[i]
    polls_temp = poll_data[priorcat == groupnames[i], ]
    priormean[i] = mean(polls_temp[ ,proploc])
    priorvar[i] = var(polls_temp[, proploc])
  }
  newdf = cbind(poll_data, priorcat)
  return(list(new_poll_data = newdf, priorcat = priorcat, priormean = priormean, priorvar = priorvar))  
  }
  