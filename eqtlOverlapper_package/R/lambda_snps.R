#' Calculate lambda statistic (median and mean) for a vector of SNPs
#' 
#' @param snps Dataframe, with a column titled "p_snp" with SNP p-values
#'
#' @return Vector of median and mean lambdas for that vector of p-values
#'
#' @author Chloe X Yap
#'
#' @export

#Calculates mean and median lambda
lambda_snps <- function(snps) {

  # Calculate lambda
  pval <- as.numeric(snps$p_snp)
    Z <- qnorm(pval/2, lower.tail = FALSE)
    #if p<1e-100, will get "Inf" --> replace with what you get if you have 1e-100
    Z[which(Z == "Inf")] <- qnorm(1e-100/2, lower.tail = FALSE) 

  # Z score to lambda
  lambda_med <- median(na.omit(Z)^2)/0.456
  lambda_mean <- mean(na.omit(Z)^2)/1

  lambda <- c(lambda_med, lambda_mean)

return(lambda)
}