#' QQ plot points to overlay onto the "base" plot from pQQ2.fun.lite()
#'
#' N.B. This code is adapted from Haplin::pQQ
#' \url{http://www.inside-r.org/packages/cran/Haplin/docs/pQQ}
#' It has only been slightly adapted to suit the purposes of this pipeline.
#' It is placed within this package for convenience of access.
#'
#' This function is also used within the following functions:
#' \code{\link{qqplot_comp}}
#' \code{\link{sep_qq}}
#'
#' @param pvals Vector of p-values
#' @param resamp TRUE or FALSE; 
#'    If resamp = TRUE, get smaller points
#'    If resamp = FALSE, this function is meant to be called within the 
#'    qqplot_comp() function. This will make points bigger, and will assign
#'    GWAS SNPs corresponding to different transprobe regions different colours
#' @param colour Colour for points. Default is "orange", though if
#'    resamp = TRUE, different colours will be assigned per transprobe region.
#' @param subset.qq Only relevant when resamp = FALSE. Labels GWAS SNP with the 
#'    lowest p-value in each transprobe region with its rsid.
#' @param nlabs.2 Only relevant when resamp = FALSE. Determines how many rsid 
#'    labels to have (see subset.qq)
#' @param pch_point Point size. Default is '.', - tiny point for resamp = TRUE
#'
#' @return QQ plot points, overlaid on "base" QQ plot from pQQ2.fun.lite()
#'
#' @author adapted from Haplin::pQQ
#'
#' @export

  points_pQQ2.fun <- function (pvals,
                              resamp = TRUE,
                              colour = "orange", 
                              subset.qq = NA,
                              nlabs.2 = NA,
                              pch_point = '.', ...) {
    
  if (any(is.na(pvals))) 
    stop("Missing p-values not allowed", call. = F)
  sim <- F
  .pvals <- sort(pvals)
  .logpvals <- -log10(.pvals)
  .order <- seq(along = .pvals)
  .n <- length(.pvals)
  .aux <- .order/(.n + 1)
  .logaux <- -log10(.aux)
  .xt <- 1/2
  .fix.order <- c(.xt, .order)
  .fix.logaux <- -log10(.fix.order/(.n + 1))
  .b.median <- qbeta(p = 1/2, shape1 = .fix.order, shape2 = .n - 
                       .fix.order + 1)
  .f.median <- approxfun(x = .fix.logaux, y = -log10(.b.median))
  .adj.logaux <- .f.median(.logaux)
  .adj.fix.logaux <- .f.median(.fix.logaux)

  if (resamp == FALSE) {

    points(.adj.logaux, .logpvals, pch = 16, #larger point
    cex = 0.7, col = as.character(colour))
    #Under construction
    if (nlabs.2 > 0) {
      .ind <- seq(length.out = nlabs.2)
      text(.adj.logaux[.ind], .logpvals[.ind] - 0.05, 
        labels = subset.qq[which(pvals == min(pvals)),"rsid"], 
           srt = 0, font = 1, cex = 0.7, adj = 0, col = "gray50", 
           xaxs = "i", yaxs = "i")
    }

  } else if (resamp == TRUE) {

    points(.adj.logaux, .logpvals, pch = pch_point, #smaller point
    cex = 0.7, col = colour)
    
    }

  return(invisible())
  }