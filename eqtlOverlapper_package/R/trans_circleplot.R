#' Produce a circle plot showing trans-eQTL SNP/probe associations.
#' 
#' @param df Dataframe of 
#'  Column names must include:
#'    1) "rsid"
#'    2) "snpchr"
#'    3) "snpbp"
#'    4) "probe"
#'    5) "probechr"
#'    6) "probebp"
#'    7) "qqid"
#'    8) "p_eqtl" N.B. this is the p-value for the trans-eQTL association
#'    9) "p_snp"
#' @param gws_snps Dataframe with GWS SNPs
#'  Column names must include:
#'    1) "rsid"
#'    2) "p_snp"
#'    3) "snpchr"
#'    4) "snpbp"
#'  If input is NA, then GWS SNPs will not be shown on the circle plot.
#' @param id_order If qq_plot() function was run, input out[["id-ordered.row"]]
#'  Vector of row numbers corresponding to the original "df" order, but ordered 
#'  by descending -log10 p-value of the sentinel trans-gene SNP. 
#'  Eg. c(3, 1, 2) if the trans-eQTL association represents that the trans-eQTL
#'  association in the 3rd row of df had the strongest sentinel SNP (-log10(p))
#' @param hg_chr Use for GRanges object seqname
#'  Extract using code:
#'  hg19chr <- fetchExtendedChromInfoFromUCSC("hg19")
#'  seq_name <- hg19chr[1:22,"NCBI_seqlevel"]
#' @param hg_chr Use for GRanges object seqlength
#'  Extract using code:
#'  hg19chr <- fetchExtendedChromInfoFromUCSC("hg19")
#'  seq_length <- hg19chr[1:22, "UCSC_seqlength"]
#' @param gen_name Use for GRanges genome name. Default = "hg19"
#' @param legend Whether to put legend on the circle plot
#'  One from: c(FALSE, TRUE, "both")
#'
#' @importFrom ggbio circle
#' @importFrom colorspace rainbow_hcl
#'
#' @return Circle plot showing for each trans-eQTL association the positions of 
#'  the SNP (red point) and probe (green line). The links connecting these
#'  associations are colour-coded to match the transprobe regions on the
#'  compartmentalised QQ plot. Their size is proportional to the strength of the 
#'  eQTL association. If an additional dataframe of GWS SNPs and their  
#'  positions are input into the gws_snps argument, they will be plotted as 
#'  a Manhattan plot on a track.
#'
#' @author Chloe X Yap
#'
#' @export

circle_plot <- function(df,
                        gws_snps = NA, #if gws_snps = NA ... else ... 
                                       #use variable (gws_snps) passed to gws_snps. 
                                       #gws_snps needs positional coordinates. 
                                       #Can extract from FDb. Colnames "chr", "bp"
                        id_order = qq_out[["id_ordered.row"]],
                        hg_chr = seq_name, 
                        hg_len = seq_length,
                        gen_name = "hg19") {

#------------------------------------------------------------------------------
# 1. Check libraries
#------------------------------------------------------------------------------

  require(ggbio)
  require(colorspace)
  require(GenomicRanges)

#------------------------------------------------------------------------------
# 2. trans-eQTL orders, colours need to be correct
#------------------------------------------------------------------------------

  #Get colour vector
  colouryay <- rainbow_hcl(nrow(df), start = 0, end = 280, c = 100, l = 65)
  
  #Order df by trans-gene sentinel SNP strength (-log10(p))
  df <- df[as.numeric(id_order),]

#------------------------------------------------------------------------------
# 3. Convert positional coordinates to GRanges objects
#------------------------------------------------------------------------------

  # Make base build to use as base of circle plot
  build.ir <- IRanges(start = 1, end = hg_len)
  build.gr <- GRanges(seqnames = as.integer(hg_chr), strand = NA, ranges = build.ir)
  seqlevels(build.gr) <- hg_chr
  seqlengths(build.gr) <- hg_len
  genome(build.gr) <- gen_name

  # Make gene.gr table - probe coordinates
  ir1 <- IRanges(start = as.numeric(df[,"probebp"]), end = as.numeric(df[,"probebp"]))
  gene.gr <- GRanges(seqnames = as.integer(df[,"probechr"]), strand = NA, ranges = ir1)
    seqlevels(gene.gr) <- hg_chr
    seqlengths(gene.gr) <- hg_len
    genome(gene.gr) <- gen_name
  
  # Make mut.gr table - SNP coordinates
  ir2 <- IRanges(start = as.numeric(df[,"snpbp"]), end = as.numeric(df[,"snpbp"]))
  mut.gr <- GRanges(seqnames = df[,"snpchr"], strand = NA, ranges = ir2)
    seqlevels(mut.gr) <- hg_chr
    seqlengths(mut.gr) <- hg_len
    genome(mut.gr) <- gen_name
  
  # Make to.gr - makes links on circle plot
  ir3 <- IRanges(start = as.numeric(df[,"probebp"]), end = as.numeric(df[,"probebp"]))
  to.gr <- GRanges(seqnames = df[,"probechr"], strand = NA, ranges = ir3)
    seqlevels(to.gr) <- hg_chr
    seqlengths(to.gr) <- hg_len
    genome(to.gr) <- gen_name
    
   #Add metadata #N.B. neglog10 p-values are for the **TRANS-EQTL ASSOCIATION*
  elementMetadata(mut.gr) <- DataFrame(assoc.pval.neglog10 = -log10(df[,"p_eqtl"]), 
    to.gr = to.gr, 
    trans_eQTL.ID = c(as.character(df[,"qqid"])), 
    leadSNP.pval = -log10(df[,"p_snp"]))

  mut.gr$trans_eQTL.ID <- factor(mut.gr$trans_eQTL.ID, 
    levels = unique(rev(mut.gr$trans_eQTL.ID)))


#------------------------------------------------------------------------------
# 4. Generate circle plots
#------------------------------------------------------------------------------

# 4.1 Plots
# 4.1.1: Plot if gws_snps variable is input (want Manhattan plot on track)
if (is.na(gws_snps) != TRUE) {

    # Only take gws_snps rows on autosomal chromosomes
    gws_snps <- gws_snps[which(gws_snps$snpchr %in% 1:22),]
    #Taking p-values for GWS SNPs
    gws_snp_logp <- -log10(gws_snps[,"p_snp"])

    #Take GWS SNP positional information
    ir4 <- IRanges(start = gws_snps[,"snpbp"], end = gws_snps[,"snpbp"])
    snps.gr <- GRanges(seqnames = as.integer(gws_snps[,"snpchr"]), strand = NA, ranges = ir4)
      seqlevels(snps.gr, force = TRUE) <- hg_chr #force drop chr23 results
      seqlengths(snps.gr) <- hg_len
      genome(snps.gr) <- gen_name
    
    elementMetadata(snps.gr) <- DataFrame(gwssnp.pval = gws_snp_logp)

    # Circle plot
    p <- ggbio() +
     scale_colour_manual(values = colouryay) + 
     circle(mut.gr, geom = "link", linked.to = "to.gr", 
        aes(color = trans_eQTL.ID, size = assoc.pval.neglog10)) +
     circle(snps.gr, geom = "point", aes(y = gwssnp.pval), size = 1, alpha = 0.5, grid = TRUE, 
        grid.background = "lightblue", grid.n = 10, trackWidth = 15, radius = 40,
        ylim = c(min(gws_snp_logp), max(gws_snp_logp))) + 
     circle(mut.gr, geom = "point", aes(y = leadSNP.pval), fill = "red",
        size = 3, grid = FALSE, trackWidth = 10, radius = 40, color = "red",
        ylim = c(min(gws_snp_logp), max(gws_snp_logp))) +         
     circle(gene.gr, geom = "rect", color = "darkgreen", trackWidth = 10, radius = 40) +
     circle(build.gr, geom = "ideo", fill = "orange") +
     circle(build.gr, geom = "scale", size = 2) +
     circle(build.gr, geom = "text", aes(label = seqnames), vjust = 10, size = 3)

# 4.1.2: Plot if gws_snps variable is not input

  } else if (is.na(gws_snps) == TRUE) {

    p <- ggbio() +
    scale_colour_manual(values = colouryay) + 
    circle(mut.gr, geom = "link", linked.to = "to.gr", 
        aes(color = trans_eQTL.ID, size = assoc.pval.neglog10)) + 
    circle(mut.gr, geom = "point", aes(y = leadSNP.pval), grid = TRUE, 
        trackWidth = 5, radius = 40, color = "red", ylim = c(0,max(mut.gr$leadSNP.pval))) + 
    circle(gene.gr, geom = "rect", color = "darkgreen", trackWidth = 5, radius = 40) +
    circle(build.gr, geom = "ideo", fill = "orange") +
    circle(build.gr, geom = "scale", size = 2) +
    circle(build.gr, geom = "text", aes(label = seqnames), vjust = 10, size = 3)

  }

    
return(p)
}