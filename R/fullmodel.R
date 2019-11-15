#big model program

iterativeelectionmodelprop = function(poll_data, stateloc, proploc, varloc, nloc, invgamma, VO = NULL, a0 = NULL, b0 = NULL, election_data, cutoffs = c(-.2,-.1, -0.025, 0.025, .1, .2), groupnames = c("Strong Red", "Red", "Lean Red", "Competitive", "Lean Blue", "Blue", "Strong Blue"), priorweights =  NULL){
  #step 1 get prior assignments
  assignments = getpriorassign(election_data, cutoffs, groupnames, priorweights)
  categories = length(groupnames)
  priormean = rep(NA, categories)
  priorvar = rep(NA, categories)
  #step 2 loop through groups and get posteriors
  for(i in 1:categories){
    polls_temp = subset(new_poll_data, assignments == groupnames[i])
    statesingroup = unique(polls)
    priormean[i] = mean(polls_temp[,proploc])
    priorvar[i] = var(polls_temp[, proploc])
    
  }
}