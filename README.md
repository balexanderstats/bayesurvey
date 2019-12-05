# bayesurvey
This package implements Bayesian conjugate prior models with a focus on handling survey data.  The models involve a empirical bayes approach that pools polls from similar states to set the hyperparameters. It is designed for predicting American Presidential Election The methodology comes from: "Poll-Based Bayesian Models to Predict United States Presidential Elections". B. Alexander, L. Ellingson, Joint Statistical Meetings Proceedings 2019. It is available [here.](http://www.balexanderstatistics.com/wp-content/uploads/2019/12/JSMProceedingsAlexander19.pdf). This project focuses on implementing 2 models from that paper that showed promise:  The Gaussian Iterative and Gaussian Polls model.  

While it was designed for American election modelling it contains various functions useful for other applications. Many of these functions rely on each other and allow flexibility to adapt the method with a slightly different setup. A different prior can be specified and a inverse gamma prior can be added for the variance of the estimate.

This package also includes polling data taken from [Huffington Post Pollster](https://elections.huffingtonpost.com/pollster).  There are data frames containing all polls for all states on the site for the 2008, 2012, and 2016 Presidential elections. National polls were not included to simplify the process of making state level prediction.  

This package is currently located on GitHub. It can be installed via the devtools package.

``` r
# install.packages("devtools")  #if devtools is not installed on your machine
devtools::install_github("balexanderstats/bayesurvey")
```
This package has the following functions:

average_error:  This function finds the average error of predictions and works for any predictions even multivariate predictions.

addcategorytopolls:  This function applies a cutoff based cluster assignment method like implemented in the paper and adds the assigned to the dataframe of poll data. This makes estimating prior parameters easier. The number of groups and the cutoffs are changable.

fitbetavecMOM:  While the beta based model in the paper is currently not implemented this function calculates the method of moment estimators of a beta distribution given a data set.

fitnormalvecMLE:  This function fits the MLE of a normal distribution given data.

getpriorassign: This function takes election data and assigns the states (or other geographical area) to categories of similar states, and returns the category each state is assigned.

iterativegaussianmodelprop:  This function is a complete implementation of the Iterative Gaussian Proportional model in the paper given poll and election data.

noniterativegaussianmodelprop: This function is a complete implementation of the Gaussian Polls Proportional model in the paper given poll and election data.

propnorm: This function takes a dataframe and a vector of locations of numeric data and proportionally normalizes it so that the rows sum to 1. This is designed for proportions.

propnormdf: This function takes a dataframe and adds columns with the normalized data at the end of the data frame.

propnormreplace:  This function replaces the data with the normalized version. This is not recommended for use by itself.

rmse:  This function calculates the root mean square error of a prediction.  This function works on any data.

unigausscp:  This function is a univariate gaussian conjugate prior with or without an inverse gamma prior on the variance.  It works with any data.

unigaussiterative:  This function is a univariate gaussian conjugate prior for proportion data.  It estimates the data variance with the variance of a sample proportion.



