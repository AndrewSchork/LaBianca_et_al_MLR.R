## Multivariate multinomial logistic regression as applied in LaBianca et al. to profiles of PGS
## Multivariate to note multiple variables of interest (PGS) being fit jointly, 
##    in addition to jointly fit covariates
## Coded for 3 groups - the adult vs. child diagnosis scenario from the paper


## Required packages
library( nnet )
library( data.table )

## Some functions

HET.v <- function( betas, hessian, variable, outcomes ) {

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
  return(nlog10pval)
}

## Set up input / edit this section

# output file / .txt

outFilePath <- "path/to/an/outfile.txt"

# Assumes all data is in one data.table, with informative column names

data <- data.table( yourData )

# Define outcome / column name in data
#   Outcome should be numbered groups
#   e.g., 0,1,2 for controls, adult ADHD, or child ADHD, respectiveely

vo <- "outcome"

# Define Variables of interest / column names in data

voi <- c( "voi1", "voi2", "voi3", ... ) 

# Define Variables of no interest (e.g., covariates) / column names in data

voni <- c( "cov1", "cov2", "cov3", "cov4", "cov5", ... )

## Run analysis / don't edit for 4 groups analysis

# Scale variables
for ( var in c( voi,voni ) ) {
  data[ ,var ] <- scale( data[ ,var ] )
}

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

for ( var in voi ) {

  print( var )
  temp.hetp.12 <- 10^-HET.v( out.betas[ c(1,2), var ], out.model$Hessian, var, c("1","2") )
  
  out.hetp.12 <- c( out.hetp.12, temp.hetp.12 )

}

## Global p: joint significance of all voi for variability simultaneously across all groups

out.anova <- anova( null.model, out.model )
out.p.global.all <- out.anova[[ 7 ]][ 2 ]

## Global p: marginal significance of each voi testing for variability simultaneously across all groups

out.p.global.each <- NULL

for ( var in voi ) {
  
  temp.voi <- voi[ voi != var ]
  temp.null.expr <- as.formula( paste( vo, 
                              paste( paste( temp.voi, collapse=" + " ), paste( voni, collapse=" + " ), sep=" + " ), 
                              sep=" ~ " ) )
  temp.null.model <- multinom( temp.null.expr, data=data, Hess=T, maxit=500 )
  temp.anova <- anova( temp.null.model, out.model )
  temp.p.global.each <- temp.anova[[ 7 ]][ 2 ]

  out.p.global.each <- c( out.p.global.each, temp.p.global.each )

}

## Build output

output.voi <- rbind( 
                    cbind( c( "Beta.1", "Beta.2" ), signif( out.betas,4 ) ),
                    cbind( c( "se.1", "se.2" ), signif( out.se,4 ) ),
                    cbind( c( "z.1", "z.2" ), signif( out.z,4 ) ),
                    cbind( c( "p.1", "p.2" ), signif( out.p,4 ) ),
                    c( "het.p.12", signif( out.hetp.12,4 ) ),
                    c( "p.marginal.global", signif( out.p.global.each,4 ) ),
                    c( "p.joint.global", rep( signif( out.p.global.all,4 ), length(out.p.global.each) ) )
                   )

## Save output to outFilePath

write.table( output.voi, file=outFilePath, row.names=F, col.names=T, quote=F )










