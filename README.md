This package has been created to perform a multivariate 
ordinal regression using a Bayesian approach. In cases 
that the response variables of a system are ordinal in nature,
 (e.g. Likert items), it returns the Bayesian posterior density 
 distributions of the coefficients of explanatory variables in
 a multivariate regression. After obtaining the distributions in R,
 a user can simply calculate 95% Highest Density Interval (HDI) for 
 each of the coefficients and if the interval does not include zero, 
 the respective explanatory variable is statistically significant.
 There are two functions available in this package as follows: “mor” and “rtruncnorm”.