#' Find trans-eQTLs and cister-eQTLs
#' 
#' @param eqtlfile Dataframe of eqtls. Cols and classes described in eqtl_cols.
#' @param eqtl_cols Vector of eqtlfile column indices in order: 
#'    1) SNP rsid     (colClass: character)
#'    2) eQTL p-value (colClass: numeric)
#'    3) SNP chr      (colClass: integer)
#'    4) probe chr    (colClass: integer)
#'    5) probe ID     (colClass: character)
#'    6) gene name    (colClass: character).
#' @param query_snps Vector of query SNPs (class: character); 
#'    if already extracted, input "all" (default).
#' @param sig Numeric: eQTL p-value significance threshold. Default = 1e-6
#' @param resamp TRUE or FALSE;
#'    FALSE value creates tables of sig eQTLs, trans-eQTLs, cister-eQTLs and 
#'    trans-eQTL IDs in out_dir folder. Default is TRUE;
#'    TRUE does not create any tables; Use FALSE with resampling analysis.
#' @param write TRUE or FALSE;
#'    Variable only required if resamp = FALSE;
#'    TRUE value writes all sig eQTLs, trans-eQTLs, cister-eQTLs and 
#'    trans-eQTL IDs to the specified out_dir. 
#'    FALSE value means that out_dir and name are not required variables.
#' @param out_dir Character: output directory. Default = "./"
#' @param name Character: Used for naming output files 
#'    Eg. phenotype. Default = "out"
#' @param keep_allsig TRUE or FALSE;
#'    TRUE value stores all eQTLs with p-value < "sig" in the output list.
#'    FALSE recommended for resampling analysis; Default is FALSE.
#'
#' @return If resamp = FALSE; return list with entries for
#'    all sig eQTLs, trans-eQTLs, cister-eQTLs, and trans-eQTL IDs 
#'    (N.B. trans-eQTL IDs only show trans-eQTLs with autosomal probes);
#'    If resamp = TRUE, return list with trans-eQTLs associations, 
#'    number of cis-eQTLs,
#'    number of trans-eQTL SNPs.
#'    In addition, keep_allsig = TRUE returns a data frame of all significant 
#'    eQTLs as a list entry.
#'
#' @examples
#' trans_cister_eqtls(eqtl_df, 
#'                    eqtl_cols = c(1:6), 
#'                    query_snps = "all", 
#'                    sig = 1e-6, 
#'                    out_dir = "Output",                  
#'                    name = "out",
#'                    resamp = FALSE,
#'                    keep_allsig = FALSE)
#'
#' @author Chloe X Yap
#'
#' @export

trans_cister_eqtls <- function(eqtlfile,
                                eqtl_cols = c(2,15,1,10,7,9), #rsid_p_snpc_probec_probe_gene
                                query_snps = "all",
                                sig = 1e-6,
                                resamp = FALSE,                                 
                                write = TRUE,
                                out_dir = ".",
                                name = "out",
                                keep_allsig = FALSE) {

  list_out <- list()

      rsid      <- eqtl_cols[1]
      p         <- eqtl_cols[2]
      snpchr    <- eqtl_cols[3]
      probechr  <- eqtl_cols[4]
      probe     <- eqtl_cols[5]
      gene      <- eqtl_cols[6]

      colnames(eqtlfile)[rsid]      <- "rsid"
      colnames(eqtlfile)[p]         <- "p_eqtl"
      colnames(eqtlfile)[snpchr]    <- "snpchr"
      colnames(eqtlfile)[probechr]  <- "probechr"
      colnames(eqtlfile)[probe]     <- "probe"
      colnames(eqtlfile)[gene]      <- "gene"

      list_out <- list()

      #Match eSNPs to query_snps
      if (query_snps == "all") {
          sig_overall.1 <- eqtlfile
        } else {
          sig_overall.1 <- eqtlfile[which(eqtlfile[,rsid] %in% query_snps),]
        }

      #Significant eQTL associations
      sig_overall.2 <- sig_overall.1[which(sig_overall.1[,p] < as.numeric(sig)),]

      #Classify trans vs. cis eQTLs
      sig_overall.3 <- sig_overall.2
      sig_overall.3$Diff_chromo <- ifelse(sig_overall.2[,snpchr] != sig_overall.2[,probechr], "YES", "NO")
      sig_overall.3$Trans_eQTL <- sig_overall.3$Diff_chromo

      #Table of all significant eQTLs
      if (keep_allsig == TRUE) {
        list_out[["sig_eqtls"]] <- sig_overall.3      
      } else if (keep_allsig == FALSE) {
        list_out[["sig_eqtls"]] <- "NA as input was keep_allsig = TRUE"      
      }

      #Just trans-eQTLs
      trans_eQTLs <- sig_overall.3[which(sig_overall.3$Trans_eQTL == "YES"),]

      #Cister eQTLs
      cis <- sig_overall.3[which(sig_overall.3$Diff_chromo == "NO"),]
      cister_eqtls <- cis[which(cis[,rsid] %in% trans_eQTLs[,rsid]),]

    #If function is not being used for resampling
    if (resamp == FALSE) {

      #Trans-eQTLs and cister-eqtls to list_out
      list_out[["trans_eqtls"]] <- trans_eQTLs
      list_out[["cister_eqtls"]] <- cister_eqtls

      #Make trans-eqtl identifiers for later use
      trans_eQTLs_aut <- trans_eQTLs[which(as.integer(trans_eQTLs[,probechr]) < 23),]
      trans_eQTLs.id <- trans_eQTLs_aut[,c(probe,gene)]
      trans_eQTLs.id$qqid <- paste(paste(paste(trans_eQTLs_aut[,rsid],trans_eQTLs_aut[,snpchr], sep = " "),
                                               paste(trans_eQTLs_aut[,gene],trans_eQTLs_aut[,probechr], sep = " "), sep = " / "), 
                                               paste(" p=", trans_eQTLs_aut[,p], sep = ""), sep = ",") #legend [,5], p-value should distinguish
      trans_eQTLs.id$qqid <- as.character(as.matrix(trans_eQTLs.id$qqid)) #just to eliminate factor problems
      colnames(trans_eQTLs.id) <- c("probe", "gene", "qqid")
      list_out[["trans_eqtls.id"]] <- trans_eQTLs.id

      #Create tables if resamp == FALSE and write == TRUE
      if (write == TRUE) {
        setwd(out_dir)
        write.table(sig_overall.3, file = paste(name, "sig_cis_trans_eQTLs.txt", sep = "_"), 
          row.names = F, col.names = T, quote = F)
        write.table(trans_eQTLs, file = paste(name, "trans_eQTLs.txt", sep = "_"), 
          row.names = F, col.names = T, quote = F)        
        write.table(cister_eqtls, file = paste(name, "cister_eqtls.txt", sep = "_"), 
          row.names = F, col.names = T, quote = F)        
        write.table(trans_eQTLs.id, file = paste(name, "trans_eQTLs_id.txt", sep = "_"), 
          row.names = F, col.names = T, quote = F)                
      }

    #If function is being used for resampling
    } else if (resamp == TRUE) {

      if (nrow(trans_eQTLs) != 0) {

        list_out[["trans_eqtls"]] <- trans_eQTLs
        list_out[["trans_eqtls_num"]] <- nrow(trans_eQTLs)        
        list_out[["trans_esnps_num"]] <- length(unique(trans_eQTLs[,rsid]))
        list_out[["cister_eqtls_num"]] <- nrow(cister_eqtls)
        list_out[["cister_esnps_num"]] <- length(unique(cister_eqtls[,rsid]))
        if (is.null(length(unique(cister_eqtls[,rsid])))) {
        list_out[["cister_esnps_num"]] <- NA
        }

      } else if (nrow(trans_eQTLs) == 0) {

        list_out[["trans_eqtls"]] <- matrix(data = NA, nrow = 1, ncol = ncol(eqtlfile))
        list_out[["trans_eqtls_num"]] <- 0 
        list_out[["trans_esnps_num"]] <- 0
        list_out[["cister_eqtls_num"]] <- NA        
        list_out[["cister_esnps_num"]] <- NA

      }

    }

  return(list_out)

}


