#' Group Assignment Function
#'
#' This function creates the various group poolings associated with each set partition
#'
#' @param npg vector of the number of observations in each group
#' @param part.vec vector of the partition of interest, as specified by the setparts() function
#' @export
#' @examples
#' grp.asn()
grp.asn <- function(npg, part.vec){
  nobs    <- sum(npg)
  numcol  <- length(unique(part.vec))
  grps    <- factor(unique(part.vec))
  maxgrps <- length(part.vec) 
  desmat  <- matrix(0, nrow = nobs, ncol = numcol) 
  
  predesmat <- matrix(0, nrow = nobs, ncol = maxgrps) 
  
  for(i in 1:maxgrps){
    if (i == 1){
      predesmat[,i] <- rep(c(1,0), times = c(npg[i], nobs - npg[i]))
    } else {
      predesmat[,i] <- rep(c(0,1,0), times = c(sum(npg[1:(i-1)]), npg[i], nobs - npg[i] - sum(npg[1:(i-1)]))) 
    }
  }
  
  for(i in 1:numcol){
    if (is.vector(predesmat[ ,part.vec == grps[i]])) {
      desmat[,i] <- predesmat[ ,part.vec == grps[i]]
      desmat[,i] <- as.numeric(levels(grps))[i]*desmat[,i]
    } else {
      desmat[,i] <- rowSums(predesmat[ ,part.vec == grps[i]])
      desmat[,i] <- as.numeric(levels(grps))[i]*desmat[,i]
    }
  }
  grpvec <- rowSums(desmat)
  
  #grpinfo <- list(Vec = grpvec)
  return(grpvec)
}