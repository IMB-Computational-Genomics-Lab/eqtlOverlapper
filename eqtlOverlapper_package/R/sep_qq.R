#' Plot an individual QQ-plot for each trans-gene region, and produce array.
#'
#' BEFORE using this function, pQQ2.fun.lite() and points_pQQ2.fun()
#' need to be loaded. These functions are available in this package:
#' \code{\link{pQQ2.fun.lite}}
#' \code{\link{points_pQQ2.fun}}
#' In turn, these QQ plot function are slightly-modified versions 
#' of Haplin::pQQ
#' \url{http://www.inside-r.org/packages/cran/Haplin/docs/pQQ}
#' 
#' @param snps Dataframe of trans-gene SNPs;
#'    Require columns with names:
#'    1) "p_snp" - for p-value   (colClass: numeric)
#'    2) "qqid" - for identifier (colClass: character)
#'    3) "rsid" - for SNP ID     (colCLass: character)
#' @param id_ordered Vector of qqids. Required to match transprobe SNP p-values
#'    to the probe. Default: "qq_out[["id_ordered"]]$qqid".
#' @param colouryay Vector of colour IDs corresponding to each 
#'    transprobe region. Length must match length(id_ordered).
#'    Default: use directly from qq_plot(resamp = F) output so colours match 
#'    ie. "qq_out[["colouryay"]]"
#'
#' @return Arrays of QQ plots. Colours used will match those representing each
#'  trans-gene region from the qq_plot function. QQ plots into automatically 
#'  be grouped into pages of 6 or 4 trans-eQTL associations, depending on 
#'  if the total number of associations is better divided by 6 or 4.
#'
#' @author Chloe X Yap
#'
#' @export

sep_qq <- function(snps,
                  id_ordered = qq_out[["id_ordered"]]$qqid,
                  colouryay = qq_out[["colouryay"]]) {

print("... Plotting compartmentalised QQ plot")

#1.0 PREP FOR GRPAH

if (length(colouryay) > 4) {

    #3.1 Figure out structure of panels
      image_quant.pre <- unique((1:(length(colouryay) %/% 6))+1)
      image_quant = NULL
      for (t in (1:length(colouryay) %/% 6)+1) {
        image_quant[t] <- ((t-1)*6)+1 #get index of each group eg.1,7,13,25
      }
      remainder <- length(colouryay) %% 6

    #If there is a remainder of 3 or more, put into groups of 6
    if (remainder > 2 | remainder == 0) {

      #Make images
        for (s in 1:(length(image_quant)-1)) {

          pdf(file = paste(name, "_QQ_colour_split_", image_quant.pre[s]-1, 
            "of", ceiling(length(colouryay)/6), ".pdf", sep = ""), 
            height = 8.27, width = 11.69)
          par(mfrow = c(2,3)) #no plot.new margin error

          for (r in image_quant[s]:(image_quant[s]+5)) { #separates images into 6s
            pQQ2.fun.lite(snps[,"p_snp"], main = paste(id_ordered[r], sep = " "))
            points_pQQ2.fun(snps[which(snps$qqid == id_ordered[r]),"p_snp"], resamp = FALSE,
              colour = colouryay[r], nlabs.2 = 0)
          }

          dev.off()  
        }

      #Make for remainder
        pdf(file = paste(name, "_QQ_colour_split_", image_quant.pre[s], 
          "of", ceiling(length(colouryay)/6), ".pdf", sep = ""), 
          height = 8.27, width = 11.69)
        par(mfrow = c(2,3)) #no plot.new margin error

        for (r in (image_quant[s]+1):(image_quant[s]+remainder)) {
          pQQ2.fun.lite(snps[,"p_snp"], main = paste(id_ordered[r], sep = " "))
          points_pQQ2.fun(snps[which(snps$qqid == id_ordered[r]),"p_snp"], resamp = FALSE,
            folour = colouryay[r], nlabs.2 = 0)         
        }

        dev.off()

    } else { #ELSE put into groups of 4

      #Reorganise structure
        image_quant.pre <- unique((1:length(colouryay) %/% 4)+1)
        image_quant = NULL
        for (t in (1:length(colouryay) %/% 4)+1) {
          image_quant[t] <- ((t-1)*4)+1 #get index of each group eg.1,5,9,13
        }
        remainder <- length(colouryay) %% 4

      #Make images
        for (s in 1:(length(image_quant)-1)) {
          pdf(file = paste(name, "_QQ_colour_split_", image_quant.pre[s]-1, 
            "of", ceiling(length(colouryay)/4), ".pdf", sep = ""), 
            height = 8.27, width = 11.69)
          par(mfrow = c(2,2))

          for (r in image_quant[s]:(image_quant[s]+3)) { #separates images into 4s
            pQQ2.fun.lite(snps[,"p_snp"], main = paste(id_ordered[r], sep = " "))
            points_pQQ2.fun(snps[which(snps$qqid == id_ordered[r]),"p_snp"], resamp = FALSE,
              colour = colouryay[r], nlabs.2 = 0)
          }

          dev.off()  
        }

      #Make for remainder
        pdf(file = paste(name, "_QQ_colour_split_", image_quant.pre[s], 
          "of", ceiling(length(colouryay)/4), ".pdf", sep = ""), 
          height = 8.27, width = 11.69)
        par(mfrow = c(2,2)) #no plot.new margin error

        for (r in (image_quant[s]+1):(image_quant[s]+remainder)) {
          pQQ2.fun.lite(snps[,"p_snp"], main = paste(id_ordered[r], sep = " "))
          points_pQQ2.fun(snps[which(snps$qqid == id_ordered[r]),"p_snp"], resamp = FALSE,
            colour = colouryay[r], nlabs.2 = 0)
        }

        dev.off() 
    }

} else if (length(colouryay) <= 4) {

      #Reorganise structure
        image_quant.pre <- unique((1:length(colouryay) %/% 4)+1)
        image_quant = NULL
        for (t in (1:length(colouryay) %/% 4)+1) {
          image_quant[t] <- ((t-1)*4)+1
        }
        remainder <- length(colouryay) %% 4

      #Make images
        for (s in 1:length(image_quant)) {
          pdf(file = paste(name, "_QQ_colour_split_", image_quant.pre[s], 
            "of", ceiling(length(colouryay)/4), ".pdf", sep = ""), 
            height = 8.27, width = 11.69)
          par(mfrow = c(2,2))

          for (r in 1:length(colouryay)) {
            pQQ2.fun.lite(snps[,"p_snp"], main = paste(id_ordered[r], sep = " "))
            points_pQQ2.fun(snps[which(snps$qqid == id_ordered[r]),"p_snp"], resamp = FALSE,
              colour = colouryay[r], nlabs.2 = 0)
          }
          
          dev.off()  
        }
  }
  return(invisible())  
}
