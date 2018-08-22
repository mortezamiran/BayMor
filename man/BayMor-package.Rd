\name{BayMor-package}
\alias{BayMor-package}
\alias{BayMor}
\docType{package}
\title{
\packageTitle{BayMor}
}
\description{
\packageDescription{BayMor}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{BayMor}
\packageIndices{BayMor}
#@mor \cr
   Step 1: Create a N by m matrix of responses, Y, where N is the number of observations and m is the number of dimension. 
           Please note that the first column of the matrix must correspond to the first set of observations, 
		   equation 1 of your multivariate regression, and the last column to the last set of observations, equation m. \cr Step 2: Create a vector “l” for the number of levels response variables in equations 1 through m. 
           The first element of the resulting vector must correspond to equation 1 and the last element to equation m. \cr		   
   Step 3: For explanatory variables of each set of observations, create a matrix,  X_1,X_2,\ldots,X_m. 
           X_1 refers to the matrix of explanatory variables for the first equation and X_m is the matrix of explanatory variables for the m equation.\cr  	   
   Step 4: Create a list with X_1,X_2,\ldots,X_m as X=list(X_1,X_2,\ldots,X_m).\cr
   Step 5: Create a vector for number of explanatory variables in each of the equations 1 through m. (i.e. nx=c(r1,r2,\ldots,rm)).\cr
   Step 6: Determine Nit, the number of iteration of the MCMC.\cr 
   Step 7: run the function, determine the burn-in steps, and then either find HDI, or plot the density.\cr
}


