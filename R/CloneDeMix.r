#' Decomposing the subclonal structure of tumors with two-way mixture models on copy number aberrations.
#'
#' A R package for deconvoluting subclonal copy number and mutational cellular prevalence.
#'
#' @name CloneDeMix
#' @author An-Shun Tai \email{daansh13@gmail.com}
#' @param tumor A G by N expression matrix of tumor samples, where G is the gene number and N is the sample size.
#' @param threshold The threhold of convergence. Default is 10^-5.
#' @param CNVstate A vector of considered copy number state.
#' @param normal A vector of length s. It's expression profile from the paired normal sample. If it's empty, the baseline is calculated from sample mean.
#' @param iterationC Maximum number of iterations. Default is 1000.
#' @param method The criteria for the final clone number.
#' @return A list is containing an estimated CNV matrix, MCP matrix, and the number of predicted clone number.
#' @export
#' @examples
#' data("ESCC_chr1")
#' res <- CloneDeMix(tumor=ESCC_chr1$tumor, normal=ESCC_chr1$normal)
#' head(res$CNV); head(res$MCP)
CloneDeMix <- function(tumor, normal=NULL, threshold=10^-5, iterC=10^3,CNVstate=c(0:10),method="aic"){
  if(is.null(normal)){base <- round(rowMeans(tumor)/2)
  warning("Baseline is calculated by sample mean")}else{
  base <- round(normal/2)
  }

  CNV.est <- matrix(NA,nrow(tumor),ncol(tumor))
  prop.est <- matrix(NA,nrow(tumor),ncol(tumor))

  for(j in 1:ncol(tumor)){
    X <- round(tumor[,j])

    filterU <- X/base < max(CNVstate)  ; filterL <- X/base > 0.05
    filter <- filterU & filterL

    ic <- c()
    for(n.clone in 2:max(CNVstate)){
      f <- CloneDeMix_core(dataX=X[filter],threshold=threshold,CNVstate=CNVstate,numbG=n.clone, baseline=base[filter],iterationC=iterC)
      if(method=="aic"){
        ic[n.clone-1] <- f$AIC
      }
      if(method=="bic"){
        ic[n.clone-1] <- f$BIC
      }
    }

    f.n.clone <- which.min(ic)+1
    f <- CloneDeMix_core(dataX=X[filter],threshold=threshold,CNVstate=CNVstate,numbG=f.n.clone, baseline=base[filter],iterationC=iterC)

    MatrixL <- matrix(1:(f.n.clone*length(CNVstate) ),length(CNVstate),f.n.clone)
    detCNV <- function(x){
      output <- matrix(NA,length(CNVstate),f.n.clone)
      output[MatrixL] <- x
      which.max(rowSums(output))
    }
    detR <- function(x){
      output <- matrix(NA,length(CNVstate),f.n.clone)
      output[MatrixL] <- x
      which.max(colSums(output))
    }

    CNV.est[filter,j] <- apply(f$Weigh,2,detCNV)-1
    prop.est[filter,j] <- f$r[apply(f$Weigh,2,detR)]

    CNV.est[which(filterU==F),j] <- max(CNVstate)
    CNV.est[which(filterL==F),j] <- 0
  }


  #### correction ####
  prop.est[is.na(prop.est)]<-1

  #test for non-normal
  Pvalue_n <- matrix(NA,nrow(tumor),ncol(tumor))

  for(j in 1:ncol(tumor)){
    X <- round(tumor[,j])
    Pvalue_n[,j] <- ppois(X,2*base)
  }

  Pvalue_n[Pvalue_n>0.5] <- 1-Pvalue_n[Pvalue_n>0.5]
  CNV.est.f <- CNV.est
  CNV.est.f[Pvalue_n>0.001] <- 2
  prop.est.f <- prop.est
  prop.est.f[CNV.est.f==2] <- 0

  colnames(CNV.est.f) <- colnames(prop.est.f) <- colnames(tumor)
  rownames(CNV.est.f) <- rownames(prop.est.f) <- rownames(tumor)
  n.clone <- matrix(apply(prop.est.f,2,function(x)length(unique(x))-1),ncol(tumor),1)
  colnames(n.clone) <- "clone number"
  rownames(n.clone) <- colnames(tumor)

  structure(list(CNV=CNV.est.f, MCP=prop.est.f,n.clone=n.clone))
}
