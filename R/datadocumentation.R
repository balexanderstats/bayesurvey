#' 2016 Pollster Polls
#'
#' A dataset containing state level poll data from Huffington's Post Pollster.  The data also contains new variables added 
#'
#' \itemize{
#'   \item X. An ID Variable
#'   \item Trump. The support in a poll for Trump (on a scale of 0 to 100 with 100 being 100%)
#'   \item Clinton. The support in a poll for Clinton (on a scale of 0 to 100 with 100 being 100%)
#'   \item Other. The support in a poll for Other (on a scale of 0 to 100 with 100 being 100%). This is not present in all data values.
#'   \item Undecided. The percent of undecided in a poll (on a scale of 0 to 100 with 100 being 100%). This is not present in all data values.
#'   \item poll-slug another ID variable
#'   \item survey_house.  The survey house who conducted the poll
#'   \item start_date. the start date of the poll
#'   \item end_date. the end date of the poll
#'   \item question_text. the question used in the survey
#'   \item sample_subpopulation. whether the poll was Registered Voters, Likely Voters or Adults
#'   \item observation. The sample size of the poll
#'   \item margin_of_error The provided margin of error of the poll
#'   \item mode.  How the poll was conducted
#'   \item partisanship.  If the poll was conducted by a partisan group
#'   \item partisian_affiliation.  The affiliation of the partisan group
#'   \item year.  The year the poll was conducted (2016 for all values)
#'   \item state.  The state in which the poll was conducted
#'   \item enddaysuntil. The days until the election at the end date of the poll. This variable was added to the original pollster data set.
#'   \item startdaysuntil. The days until the election at the start date of the poll.  This variable was added to the original pollster data set.
#'   \item middaysuntil. The midpoint between endddaysuntil and startdaysuntil.  This variable was added to the dataset.
#'   \item Dem2P. the proportionally normalized support for Clinton so that Dem2P+Rep2P=1
#'   \item Rep2P. the proportionally normalized support for Trump so that Dem2P+Rep2P = 1
#'   \item DemVote. the proportionally normalized election result for Clinton so that DemVote+RepVote = 1
#'   \item RepVote. the proportionally normalized election results for Trump so that DemVote + RepVote = 1
#'   \item Error. the difference between Dem2P and DemVote or the error of the poll
#'   \item Johnson. The support for Gary Johnson if provided
#'   \item McMullin.  The support for McMUllin if provied
#' }
#'
#' @docType data
#' @keywords datasets
#' @name polls2016
#' @usage data(polls2016)
#' @format A data frame with 2650 rows and 28 variables
NULL
