#' Decomposing the subclonal structure of tumors with two-way mixture models on copy number aberrations.
#'
#' The core function of CloneDeMix.
#'
#' @name CloneDeMix_core
#' @author An-Shun Tai \email{daansh13@gmail.com}
#' @param dataX Gene expression of a tumor sample.
#' @param threshold The threhold of convergence. Default is 10^-5.
#' @param CNVstate A vector of considered copy number state.
#' @param numbG The number of subclones.
#' @param baseline A vector of length s. It's expression profile from the paired normal sample. If it's empty, the baseline is calculated from sample mean.
#' @param iterationC Maximum number of iterations. Default is 1000.
#' @return A list containing MCP vector, the weight of copy number state, error rate, AIC, and BIC.
#' @export
#' @examples
CloneDeMix_core <- function(dataX,threshold=10^-5,CNVstate=c(0:4),numbG=4, baseline=NULL ,iterationC=10^3){
 if(is.null(baseline)){baseline <- rep(round(mean(dataX)/2),length(dataX))
                   warning("Baseline is calculated by sample mean")}

# State <- union(CNVstate,c(2))
  State <- CNVstate

#initial value
 phi_init <- matrix(runif(length(State)*numbG,0,1),length(State),numbG) ; phi_init <- phi_init/sum(phi_init)
 r_init <- seq(0.1,0.9,0.8/(numbG-1))

 MatrixL <- matrix(1:(numbG*length(CNVstate)),length(CNVstate),numbG)


# iteration
 stop <- FALSE; count <- 1; error <- c()
 while(stop==FALSE){

  Prop_a <- matrix(NA,numbG*length(CNVstate),length(dataX))
  for(j in 1:length(CNVstate)){
  for(k in 1:numbG){
   Pr <- dpois(dataX,baseline*(2*(1-r_init[k])+CNVstate[j]*r_init[k] )) *  phi_init[j,k]
   Prop_a[MatrixL[j,k],] <- Pr +10^-100
  }}

  Wijk <- t(t(Prop_a)/colSums(Prop_a))


#estimation
  phi_hat <- matrix(rowSums(Wijk)/length(dataX),length(CNVstate),numbG)

  r_hat <- numeric(numbG)
  for(k in 1:numbG){
   r0 <- r_init[k]
   stopRR <- FALSE
   while(stopRR==FALSE){
    F1 <-  -sum( -baseline%*%t(Wijk[MatrixL[,k],]*(CNVstate-2))+dataX%*%t( Wijk[MatrixL[,k],]*(CNVstate-2)/( (CNVstate-2)*r0+2 ) ) )
    F2 <-  sum( dataX%*%t( Wijk[MatrixL[,k],]*(CNVstate-2)^2/( (CNVstate-2)*r0+2 )^2 ) )
    r1 <- r0-F1/F2
    r1[r1<0] <- 10^-5; r1[r1>1] <- 1-10^-5
    ifelse(abs(r1-r0)<threshold,stopRR <- TRUE,r0 <- r1)
   }
   r_hat[k] <- r1
  }

#evalute performance
  error[count] <- mean(abs( c(phi_init-phi_hat,
                              r_init - r_hat) ))
  if(error[count] < threshold){stop <- TRUE
  }else{
   phi_init <- phi_hat
   r_init <- r_hat
   count <- count +1
  }#if

  if(count > iterationC){stop <- TRUE
                         warning("NOT converge")}
 }#iteration

 LikH <- c()
 for(j in 1:length(CNVstate)){
 for(k in 1:numbG){

  LikH <- cbind(LikH,phi_hat[j,k]*dpois(dataX,baseline*(2*(1-r_init[k])+CNVstate[j]*r_init[k] )))

 }}
 LikH_s <- rowSums(LikH)

 Loglik <- sum(log(LikH_s+ min(LikH_s[LikH_s!=0])*10^-10))
 aic <- -2*Loglik+2*(numbG+numbG*length(CNVstate))
 bic <- -2*Loglik+log(length(dataX))*(numbG+numbG*length(CNVstate))

structure(list(phi=phi_hat, r=r_hat, Weight=Wijk, Error=error,AIC=aic,BIC=bic))
}

