## Multivariate multinomial logistic regression as applied in LaBianca et al. to profiles of PGS
## multivariate to note multiple variables of interest (PGS) being fit jointly, 
##    in addition to jointly fit covariates
## Currently coded for 4 groups - the ADHD +/- ASD scenario in thee paper


## Required packages
library( nnet )
library( data.table )

## Some functions

HET.v <- function( betas. hessian, variable, outcomes ) {

  N <- length( outcomes )
  h.inds <- paste( outcomes, variable, sep=":" )
  
  if ( N == 2 ) {
    Tmat <- matrix( c( -1,1 ), 1, 2, byrow=T )
    Tbeta <- Tmat %*% betas
    V_Hess <- Tmat %*% solve( hessian )[ h.inds, h.inds ] %*% t( Tmat )
    chi2 <- Tbeta %*% solve( V_Hess ) %*% Tbeta
    nlog10pval <- pchisq( chi2, df=1, lower.tail=F, log.p=T ) / log(10) * (-1)        
  }
  
  if ( N == 3 ) {  
    Tmat <- matrix( c( -1,1,0,0,-1,1 ), 2, 3, byrow=T )
    Tbeta <- Tmat %*% betas
    V_Hess <- Tmat %*% solve( hessian )[ h.inds, h.inds ] %*% t( Tmat )
    chi2 <- Tbeta %*% solve( V_Hess ) %*% Tbeta
    nlog10pval <- pchisq( chi2, df=2, lower.tail=F, log.p=T ) / log(10) * (-1)        
  }
}

## Set up input

# Assumes all data is in one data.table, with informative column names

data <- data.table( data )

# Define outcome / column name in data
#   Outcome should be numbered groups
#   e.g., 0,1,2,3 for controls, ADHD only, ASD only, ADHD and ASD, respectiveely

vo <- "z"

# Define Variables of interest / column names in data

voi <- c( "PRS19", "PRS22", "PRS23", "PRS10" ) 

# Define Variables of no interest (e.g., covariates) / column names in data

voni <- c( "C1", "C2", "C3", "C4", "C5" )

## Run analysis

# fit global null model

null.expr <- as.formula( paste( vo, 
                               paste( voni, collapse=" + " ), 
                               sep=" ~ " ) )
null.model <- multinom( null.expr, data=data, Hess=T, maxit=500 )

# fit full model

out.expr <- as.formula( paste( vo, 
                              paste( paste( voi, collapse=" + " ), paste( voni, collapse=" + " ), sep=" + " ), 
                              sep=" ~ " ) )
out.model <- multinom( out.expr, data=data, Hess=T, maxit=500 )

## Basic summary stats
out.betas <- summary( out.model )$coefficient[ ,voi ]
out.se <- summary( out.model )$standard.errors[ ,voi ]
out.z <- out.betas / out.se
out.p <- pnorm( -abs( out.z ) ) * 2

## Heterogeneity test summary stats

out.hetp.12 <- NULL
out.hetp.13 <- NULL
out.hetp.22 <- NULL
out.hetp.123 <- NULL

for ( var in voi ) {

  print( var )
  temp.hetp.12 <- 10^-HET.v( out.betas[ c(1,2), var ], out.model$Hessian, var, c("1","2") )
  temp.hetp.13 <- 10^-HET.v( out.betas[ c(1,3), var ], out.model$Hessian, var, c("1","3") )
  temp.hetp.23 <- 10^-HET.v( out.betas[ c(2,3), var ], out.model$Hessian, var, c("2","3") )
  temp.hetp.123 <- 10^-HET.v( out.betas[ c(1,2,3), var ], out.model$Hessian, var, c("1","2","3") )























