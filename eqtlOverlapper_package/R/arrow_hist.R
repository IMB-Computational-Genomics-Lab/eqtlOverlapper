#' Generate histograms of lambda distributions with red line indicating the
#' lambda statistic for the discovery analysis transprobe SNPs
#' 
#' @param discsnps_p Vector of discovery analysis transprobe SNPs
#' @param lambda Dataframe with columns of mean and median lambdas, per resamp
#'  Colnames as "median_lambda", and "mean_lambda"
#' @param resamp_n Integer for number of resampling performed
#' @param out_dir Output directory
#'
#' @return 1) Lambda for the discovery analysis transprobe SNPs AND
#'  2) A median and mean lambda distribution histogram with red lines
#'  for the lambda statistic from the discovery analysis transprobe SNPs,
#'  written to the current working directory.
#'
#' @author Chloe X Yap
#'
#' @export

#
arrow_hist <- function(discsnps_p, lambda, resamp_n = 1000, out_dir = ".") { 

all_lambdas.med <- lambda$median_lambda
all_lambdas.mean <- lambda$mean_lambda

print("... Plotting lambda histograms")  
  #Find lambda of actual analysis #####
  pval <- as.numeric(discsnps_p)
  Z <- qnorm(pval/2, lower.tail = FALSE)
  #if p<1e-100, will get "Inf" -->all_lambdas.mean replace with what you get with 1e-100
  Z[which(Z == "Inf")] <- qnorm(1e-100/2, lower.tail = FALSE) 

  med <- median(na.omit(Z)^2)/0.456
  mea <- mean(na.omit(Z)^2)/1

print("... Calculating ranks")
  medrank <- resamp_n-length(which(all_lambdas.med > med))
  mearank <- resamp_n-length(which(all_lambdas.mean > mea))
  medrank_leg <- paste("rank = ", medrank, "/", resamp_n, sep = "")
  mearank_leg <- paste("rank = ", mearank, "/", resamp_n, sep = "")

  actual_lambda <- data.frame(med, mea)

  p_med <- hist(all_lambdas.med)
  dev.off()
  p_mean <- hist(all_lambdas.mean)
  dev.off()

print("... Plotting median lambda") 
  #Plot histogram
  pdf(paste(out, "lambda_med_hist.pdf", sep = ""))
  hist(all_lambdas.med, xlab = expression(lambda), 
     xlim = c(0,(round((max(med+1,max(p_med[[1]]))),0))), 
     ylim = c(0,max(p_med[[2]])), 
     main = paste(disease_name, "lambda (median) histogram", sep = " "))
  legend("topright", medrank_leg , bty="n", text.col = "red")
  abline(v = med, col = "red")
  dev.off()

print("... Plotting mean lambda") 
  #Plot histogram
  pdf(paste(out, "lambda_mean_hist.pdf", sep = ""))
  hist(all_lambdas.mean, xlab = expression(lambda), 
     xlim = c(0,(round((max(mea+1,max(p_mean[[1]]))),0))), 
     ylim = c(0,max(p_mean[[2]])), 
     main = paste(disease_name, "lambda (mean) histogram", sep = " "))
  legend("topright", mearank_leg , bty="n", text.col = "red")    
  abline(v = mea, col = "red")
  dev.off()

  return(actual_lambda)
} 