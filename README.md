# pprRFA
R package providing basic functionalities to use projection pursuit regression (PPR) in context of regional frequency analysis.  

###Motivation

Regression method can be used to predict at ungauged location a hydrological variable Y based on the avalaible site characteristics

A Ridge function has the form g(a'x) where x is vector of explanatory variable,a is a direction, or a unitary vector of coefficients and g is a nonparametric function. A PPR model is regression model defines by several Ridge terms

E(Y|x) = u + g_1(a_1'x) + g_2(a_2'x) + ... + error
 
and can be seen as a generalized additive model with estimated linear predictor a'x. 

This packages provide wrapping functions and graphical tools to create a gam object of the well known mgcv R packages, with additional direction.

For further information please consult the R documentation. 

###Installation

The package is not in the CRAN repository. It can be download from the present github page and built locally.

###Version
- version v0.2: correct issues

- version v0.1: first complete version

###Note 
This package still needs further testing and an improved documentation.
