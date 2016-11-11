#' Extract trans-probe region SNPs
#' 
#' @param trans_eqtl_id Dataframe with eqtl probe and qqid;
#'    Eg. list_out[["trans_eQTLs.id"]] from trans_cister_eqtls() output;
#'    Else, colnames as:
#'    1) "Probe"   (colClass: character)
#'      And if resamp = FALSE: ...
#'    2) "qqid"    (colClass: character) - a unique trans-eQTL identifier.
#' @param gwas Dataframe of gwas summary statistics.
#'    gwas column *indices* in order: 
#'    1) rsid        (colClass: character)
#'    2) gwas p-val  (colClass: numeric)
#'    Recommend also having chr and bp columns for output completeness
#'    as columns after rsid and gwas p-val, though this is not essential. 
#'    However, if planning on LD clumping, chr columns are required.
#' @param snp_gr GRanges object with SNP position coordinates;
#'    names(snp_gr) corresponds with SNP rsIDs;
#'    Ensure snp_gr genome build matches with probe_gr file.
#' @param probe_gr GRangesList object with probe position coordinates; 
#'    names(probe_gr) corresponds with Illumina probe IDs;
#'    Where there are multiple probe coordinates, the first will be used;
#'    Ensure probe_gr genome build matches with snp_gr file.
#' @param trans_gene Numeric: extract SNPs around probe start +/- this region; 
#'    Default = 5e4.
#' @param out_dir Character: output folder, followed by /; Default is "."
#' @param name Character: phenotype name.
#' @param resamp TRUE or FALSE
#'
#' @importFrom GenomicRanges GRangesList flank findOverlaps mcols
#' @importFrom plyr join
#'
#' @return Dataframe of trans-gene SNPs with columns titled:
#'  If resamp = FALSE:
#'    1) "rsid", 2) "p_snp", 3) "probe", 4) "qqid";
#'    If gwas input had other columns, these will also be incorporated.
#'  If resamp = TRUE:
#'    1) "snpchr", 2) "p_snp", 3) "rsid"
#'    Only actually need "p_snp" and "rsid"; "snpchr" is to 
#'    make clumping faster, if desired
#'
#' @author Chloe X Yap
#'
#' @export

transprobe <- function(trans_eqtl_id,
                        gwas = gwas_pos,
                        snp_gr = snp137common,
                        probe_gr = probe_coord,
                        trans_gene = 5e4,
                        write = TRUE,
                        out_dir = ".",
                        name = "out",
                        resamp = FALSE) {
  
    # Ensure probe and gene columns are characters

    if ((nrow(trans_eqtl_id) == 0) | (length(which(is.na(trans_eqtl_id[1,]))) == ncol(trans_eqtl_id))) {

      qqplotSNPs <- data.frame("NA", 0.5)
      colnames(qqplotSNPs) <- c("rsid", "p_snp")
      
    } else {

    trans_eqtl_id$probe <- as.character(as.matrix(trans_eqtl_id$probe))

    if (resamp == FALSE) {
    trans_eqtl_id$qqid <- as.character(as.matrix(trans_eqtl_id$qqid))      
    }

      print("Taking the vector of probes, checking if any are missing")
      #Take vector of probes
      trans_probes <- trans_eqtl_id$probe[which(trans_eqtl_id$probe %in% names(probe_gr))]

      if (length(trans_probes) < length(trans_eqtl_id$probe)) {
          miss <- trans_eqtl_id$probe[-which((trans_eqtl_id$probe %in% names(probe_gr)))]
          warning(paste("WARNING:", miss, "- trans-eQTL probe/s missing from probe file."))
      }
      
      if (resamp == TRUE) {

        probesngenes <- data.frame(trans_eqtl_id$probe, stringsAsFactors = F)
        colnames(probesngenes) <- c("probe")

      } else if (resamp == FALSE) {

        probesngenes <- data.frame(trans_eqtl_id$probe, trans_eqtl_id$qqid, stringsAsFactors = F)
        colnames(probesngenes) <- c("probe", "qqid")        

      }
      probes_bioc <- probe_gr[trans_probes]
  
      #Multiple probes - only take first entry
      probes_bioc <- lapply(1:length(probes_bioc), function(x) probes_bioc[[x]][1])
      probes_bioc <- GRangesList(probes_bioc)
      names(probes_bioc) <- trans_probes
  
      #Get +/-trans_gene region
      print("Using flank() to extend the trans-gene region")      
      probes_bioc <- flank(probes_bioc, as.integer(trans_gene), both = TRUE)
      
      #findOverlaps to get SNPs within trans_gene region
      print("Finding SNPs within the trans-gene region")      
      overlapsnps_ind <- findOverlaps(snp_gr, probes_bioc, ignore.strand = TRUE, select = "all")
      
      #"Self-subset" matching SNPs. Better than subsetByOverlaps() as that only takes unique. 
        #This will take all matching
      print("Subsetting matching SNPs within the trans-gene region")
      overlapsnps <- snp_gr[queryHits(overlapsnps_ind),]
    
      #Add probe IDs as metatdata columns. Once again, use overlapsnps_ind as will get matching length
      mcols(overlapsnps) <- names(probes_bioc[subjectHits(overlapsnps_ind)])
      
      #Make dataframe of SNPid, probeid
      if (resamp == FALSE) {

      qqplotSNPs2join <- data.frame(names(overlapsnps), overlapsnps$value, stringsAsFactors = F)
      colnames(qqplotSNPs2join) <- c("rsid", "probe")

      print("Joining bioconductor SNP and gwas SNP output")      
      bioc_in_gwas <- qqplotSNPs2join[which(qqplotSNPs2join$rsid %in% gwas$rsid),]
      gwas_in_bioc <- gwas[which(gwas$rsid %in% bioc_in_gwas$rsid),]
      qqplotSNPs <- join(gwas_in_bioc, bioc_in_gwas, by = "rsid", type = "inner", match = "all")
    
      #Join by probe ID --> get Gene
      qqplotSNPs <- join(qqplotSNPs, probesngenes, by = "probe", type = "inner", match = "all")

      } else if (resamp == TRUE) {

      qqplotSNPs2join <- data.frame(names(overlapsnps), stringsAsFactors = F)
      colnames(qqplotSNPs2join) <- c("rsid")

      print("Joining bioconductor SNP and gwas SNP output")      
      bioc_in_gwas <- data.frame(qqplotSNPs2join[which(qqplotSNPs2join$rsid %in% gwas$rsid),])
      colnames(bioc_in_gwas) <- "rsid"
      gwas_in_bioc <- gwas[which(gwas$rsid %in% bioc_in_gwas$rsid),]
      qqplotSNPs <- join(gwas_in_bioc, bioc_in_gwas, by = "rsid", type = "inner", match = "all")
    
      }
    
      #Get list of unique surrounding SNPs
      if (resamp == FALSE & write == TRUE) {
        setwd(out_dir)
        write.table(qqplotSNPs, file = paste("qqplotSNPs_for_", name, ".txt", sep = ""), 
          sep = "\t", row.names = F, col.names = T, quote = F)
      }

    }

    # If there are no SNPs in the transprobe region, and are resampling, assign p=0.5
    if (resamp == TRUE & nrow(qqplotSNPs) == 0) {

      qqplotSNPs <- data.frame("NA", 0.5)
      colnames(qqplotSNPs) <- c("rsid", "p_snp")

    }

    # If there are no SNPs in the transprobe region, and not resampling, give warning
    if (resamp == FALSE & nrow(qqplotSNPs) == 0) {

      warning(paste("WARNING: no SNPs found in transprobe region."))      

    }


  return(qqplotSNPs)

  }

