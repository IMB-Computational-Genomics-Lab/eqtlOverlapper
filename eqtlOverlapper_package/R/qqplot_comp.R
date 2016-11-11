#' Plot QQ-plot with GWAS SNP p-values compartmentalised by qqid.
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
#'    1) "p_snp" - for p-value       (colClass: numeric)
#'    2) "qqid" - for identifier (colClass: character)
#'    3) "rsid" - for SNP ID     (colCLass: character)
#' @param name Character: phenotype name.
#' @param Directory for output. Default is current working directory
#'
#' @importFrom colorspace rainbow_hcl
#'
#' @return QQ plot, with SNPs compartmentalised by trans-gene region.
#'    The colours (from red to purple on the colour wheel) are ordered 
#'    by ascending -log10 p-value of the sentinel SNP in the region. 
#'    The sentinel SNPs are also labelled by their rsid on the plot.
#'    A list is also returned as output containing:
#'      1) Vector of the colour IDs used.
#'      2) Dataframe of the qqids, ordered by -log10(p-value) of the 
#'         sentinel trans-gene SNP.
#'      3) Vector of row indices, ordered by -log10(p-value) of the 
#'         sentinel trans-gene SNP. Used to construct circle_plot.
#'
#' @author Chloe X Yap
#'
#' @export

qqplot_comp <- function(snps, 
                    name = "out", 
                    out_dir = ".") {

    #1.0 PREP FOR GRAPH
    print("... Prepping for QQ plot")
      
      #1.1 Get vector of colour IDs
      colouryay <- rainbow_hcl(length(unique(snps$qqid)), start = 0, end = 280, 
        c = 100, l = 65)
      
      #1.2 Prepare legend - vector of both SNP and gene name
      id_order.id <- unique(snps$qqid) #holds qqid names
      id_order.min = NULL
      subset_qq = NULL

      #1.3 Get points for each unique qqid
      for (p in 1:length(id_order.id)) {
      
        subset <- which(snps$qqid == id_order.id[p])

        id_qq <- snps[subset,] #take out rows for each unique qqid
        id_order.min <- rbind(id_order.min, 
        id_qq[which(id_qq[,"p_snp"] == min(id_qq[,"p_snp"])),c("p_snp", "rsid")][1,]) 
            #if there are multiple minimums, only take the first in the table
      }

      #1.4 Order qqid by significance (get out the rownames to use as indices)
      id_order <- as.data.frame(cbind(id_order.id, id_order.min))
      colnames(id_order) <- c("qqid", "p_snp", "rsid")
      id_order$qqid <- as.character(as.matrix(id_order$qqid)) #qqid
      id_order$p_snp <- as.numeric(as.matrix(id_order$p_snp)) #sentinel p-value
      id_order$rsid <- as.character(as.matrix(id_order$rsid)) #rsid
      rownames(id_order) <- 1:nrow(id_order)
      
      id_ordered <- id_order[order(as.numeric(id_order$p_snp), decreasing = FALSE),] 
      #orders rows by descending p-value ie. smallest first
        #data frame of [,c(qqid, p_snp, SNP_id)]
      id_ordered.row <- rownames(id_ordered)

  #2.0 START GRAPH
    print("... Plotting whole QQ plot")

    #2.1 Make graph area
    pQQ2.fun.lite(snps[,"p_snp"], colour = "gray20", pin = c(8.27, 11.69), 
      main = paste(name, "QQ Plot", sep = " "))

    #2.2 Plot p-values onto graph area
      #Get p-values from each probe
      print("...... Adding points from each qqid")
      for (p in nrow(id_ordered):1) { #so that bigger associations = top layer
        #Get p-values for each qqid
        subset_qq <- snps[which(snps$qqid == id_ordered[p,"qqid"]),]
        subset_qq[,"rsid"] <- as.character(as.matrix(subset_qq[,"rsid"]))
        p_qq <- subset_qq[,"p_snp"] 
          #take out matching rows, p-value column

        #Add points
        points_pQQ2.fun(p_qq, colouryay[p], subset.qq = subset_qq, resamp = FALSE,
          nlabs.2 = length(which(subset_qq[,"p_snp"] == min(subset_qq[,"p_snp"]))))
    }

  #3.0 STORE VARIABLES FOR OUTPUT
    out <- list()
    out[["colouryay"]]      <- colouryay
    out[["id_ordered"]]     <- id_ordered
    out[["id_ordered.row"]] <- id_ordered.row

return(out)
}
