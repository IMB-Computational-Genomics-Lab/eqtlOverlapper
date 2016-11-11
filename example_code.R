#==============================================================================#
#==============================================================================#
#           			       						                                           #
#           			       EXAMPLE "DISCOVERY" ANALYSIS                          #
#                                                                              #
# This is an example script showing a workflow for a discovery trans-eQTLs     #
# analysis, as discussed in:                                                   #
#                                                                              #
#       "Independent trans-eQTL identified in whole blood have limited         #
#                   influence on complex disease biology"                      #
#                                                                              #
# This workflow requires functions that are available at <http://_______>.     #
#                                                                              #
# It follows 4 general steps:                                                  #
#                                                                              #
#   1) Obtaining trans-eQTLs and cis-eQTLs                                     #
#       Function: trans_cister_eqtls(..., resamp = FALSE)                      #
#                                                                              #
#   2) Extracting SNPs within the trans-probe region                           #
#       Function: transprobe(..., resamp = FALSE)                              #
#                                                                              #
#      ... (Optional LD clumping) ...                                          #
#                                                                              #
#   3) Producing a QQ plot, showing the GWAS p-values for trans-probe SNPs     #
#       Function: qqplot_comp() - all transprobe SNPs on one graph             #
#       Function: sep_qq(...) - each transprobe region = one graph             #
#       N.B. These functions BOTH rely on functions adapted from Hamlin::pQQ   #
#          a) pQQ2.fun.lite() - show aggregated discovery transprobe SNP p     #
#          b) points_pQQ2.fun() - add points per transprobe region             #
#                                                                              #
#   4) Producing a circle plot, illustrating the trans-eQTL associations       #
#       Function: circle_plot()                                                #
#                                                                              #
#                           Author: Chloe X Yap                                #
#                                                                              #
#==============================================================================#
#==============================================================================#

#-------------------------------------------------------------------------------
# 1.0 Obtain trans-eQTLs and cis-eQTLs
#-------------------------------------------------------------------------------

	# 1.1: Data - eQTLs

    # 1.1.1: Read-in eQTL file
    eqtlfile <- read.csv("eqtl.csv", header = T, as.is = T)

	# 1.2: Function

    # Notes on eqtl_cols arg: order is rsid_p_snpc_probec_probe_gene
    # If using BSGS database, try c(2,15,1,10,7,9)

    # May choose to set working directory
    # setwd("../outdir")

    # 1.2.1: Call function
	trans_eqtls <- trans_cister_eqtls(eqtlfile, eqtl_cols = c(2,15,1,10,7,9),
		query_snps = "all", sig = 1e-6, write = TRUE, out_dir = "./", 
		name = "disease", resamp = FALSE, keep_allsig = FALSE)

#-------------------------------------------------------------------------------
# 2.0 Find SNPs within the trans-gene region
#-------------------------------------------------------------------------------

	# 2.1: Data

    # 2.1.1: Data - GWAS SNPs and positional coordinates.
      # Only need rsid and p-value up to Step 4.0.
      # Step 4.0 Circle plot also requies chr, bp.
      # Therefore, retaining chr and bp is advised.
      # N.B. The column names of the dataframe are critical!

    # 2.1.1.1: Read in file
    gwas_pos <- read.table("BMI_SNP_posfile.txt", header = T, as.is = T)

    # 2.1.1.2: Clean-up, taking note of gwas_pos colnames!
	gwas_col <- c(2, 10, 1, 4) #order rsid, p-value, chr, bp
    gwas_pos <- gwas_pos[,gwas_col]
    gwas_pos <- gwas_pos[order(gwas_pos[,3], gwas_pos[,4]),]
    colnames(gwas_pos) <- c("rsid", "p_snp", "snpchr", "snpbp")

    # 2.1.2: Data - SNP positions, as a Genomic Ranges file

    # 2.1.2.1: Load FDb.UCSC.snp137common.hg19 package
    # Use provided .rds version of the SNPs for speed
    snp137common <- readRDS("snp137common.rds")
    # OR load from bioconductor. This takes some time.
	# source("http://bioconductor.org/biocLite.R")
    # biocLite("FDb.UCSC.snp137common.hg19")
    # print("... Calling features of UCSC common SNP package. Takes time!")
    # snp137common <- features(FDb.UCSC.snp137common.hg19)

    # 2.1.2.2: To reduce file, filter to snp137common entries that match GWAS
    	#As enrichment tests are limited by GWAS SNPs
    snp137common <- snp137common[which(names(snp137common) %in% gwas_pos[,1])] 

    # 2.1.3: Data - Probe coordinates, as a GRangesList

    # 2.1.3.1 Load GRangesList of probe coordinates
    probe_coord <- readRDS("IlluminaHT12v4_GRlistProbesFromBSGS_CHR6HAP.rds")

  # 2.2: Function - Extract trans-probe region SNPs

    # 2.2.1: Libraries
    require(plyr)
    require(GenomicRanges)
  
    # 2.2.2: Call transprobe()
    qqsnps <- transprobe(trans_eqtls[["trans_eqtls.id"]], gwas = gwas_pos, 
    	snp_gr = snp137common, probe_gr = probe_coord, trans_gene = 5e4, 
    	resamp = FALSE, write = TRUE, out_dir = "", name = "out")

    # 2.2.3: Write table of qqsnps

    # May choose to set working directory
    # setwd("../outdir")

    write.table(qqsnps, file = "discovery_qqplotSNPs.txt", 
      sep = "\t", row.names = F, col.names = T, quote = F)

#-------------------------------------------------------------------------------
# (Optional: LD clumping)
#-------------------------------------------------------------------------------

  # 2.2.2: To LD clump, write to a table before sending to plink
  snps2clump <- qqsnps[,c("snpchr", "p_snp", "rsid")]

  # An example clumping function is provided at <http://____>

  # Below is a conceptual outline, overwriting qqsnps for continuity:
  # Eg. qqsnps <- ld_clump_fun(snps2clump, p1p2r2kb = c(0.1, 0.1, 0.1, 1000),
  #                            clump_dir = "../clump", plink_dir = "/plink", 
  #                            bfile_dir = "../bfile", pre_exclude = 0.5)
  #     setwd("outdir") # if you moved directory for clumping ouputs
  #                     # move back to outdir to ensure continuity of pipeline

#-------------------------------------------------------------------------------
# 3.0 Plot SNPs within trans-gene region on QQ-plot
#-------------------------------------------------------------------------------

  # 3.1: Function - 1 QQ-plot, compartmentalised by transprobe region

    # 3.1.1: Libraries
    require(colorspace)

    # 3.1.2: Prepare for plot

    # May choose to set working directory
    # setwd("../outdir")

    pdf(file = paste(name, "_QQplot_colour.pdf", sep = ""), 
    	height = 8, width = 9)
    par(xpd = TRUE, mar = par()$mar + c(0,0,0,10)) #space for legend
  
    # 3.1.2: Plot QQ plot
    qq_out <- qqplot_comp(qqsnps, resamp = FALSE, out_dir = ".", name = "out")
    
    # 3.1.3: Create legend for QQ plot
    legend("right", c("All probes", qq_out[["id_ordered"]]$qqid), 
    	col = c("gray20", qq_out[["colouryay"]]), 
    	title = expression(bold("trans-eQTL ID and assoc p-value")), 
    	title.adj = 0.2, xpd = TRUE, horiz = FALSE, xjust = 0, 
    	inset = c(-0.4,0), bty = "n", pch = 16, cex = 0.75)
    
    # 3.1.4: Turn off graphics to save
    dev.off()

    # 3.1.5: May choose to keep the table of sentinel trans-gene SNPs
    # N.B. if there were multiple sentinels, only first is taken

    # May choose to set working directory
    # setwd("../outdir")

    write.table(qq_out["id_ordered"], row.names = F, col.names = T, 
    	file = paste(out_dir, name, "qqplothits.txt", sep = "_"),
    	quote = F, sep = "\t")

  # 3.2: Function - separate QQ-plots per transprobe region

    # May choose to set working directory
    # setwd("../outdir")

    # 3.2.1: Call
    sep_qq(snps = qqsnps, id_ordered = qq_out[["id_ordered"]], 
      colouryay = qq_out[["colouryay"]])

#-------------------------------------------------------------------------------
# 4.0 Circle plot
#-------------------------------------------------------------------------------
# N.B. If the number of trans-eQTLs represented on the QQ-plot and the 
# circle-plot do not match, this is likely because there were too few
# SNPs in the trans-gene region corresponding to that probe.

  # 4.1: Libraries

    # 4.1.1: Load bioconductor packages
    #source("http://bioconductor.org/biocLite.R")
    #biocLite("ggbio") #includes biocLite("IRanges and GenomicRanges")
    require(ggbio)
    require(GenomicRanges)
  
    # 4.1.2: Load biocondcutor genome data for chr names and lengths to plot.
    # Used as GRanges seqnames and seqlengths.
    # Need to match chromosome names used for SNP and probe coordinates.
    # Genome builds also need to match.
    require(GenomeInfoDb)
    hg19chr <- fetchExtendedChromInfoFromUCSC("hg19")
    seq_name <- hg19chr[1:22,"NCBI_seqlevel"]
    seq_length <- hg19chr[1:22, "UCSC_seqlength"]

    # 4.1.3: Load other libraries
    require(plyr)

  # 4.2: Cleaning and prepping data, get trans-eQTL SNP, probe pos and p-values

    # N.B. This example code only plots trans-eQTLs that do have SNPs
    # in the trans-gene region. It means that the colours used in the 
    # circle and QQ plot will match per transprobe region

    # 4.2.1: Obtain trans-eQTLs corresponding to transprobe regions (autosomal)
    circle_ind <- which(trans_eqtls[["trans_eqtls.id"]]$qqid 
      %in% qq_out[["id_ordered"]]$qqid)
    trans_eqtls_circle <- trans_eqtls[["trans_eqtls"]][circle_ind,]
  
    # 4.2.1: Obtain SNP and probe coordinates
      # SNP: Can use GWAS positional information 
      # Alternatively, could draw upon snp137common (see earlier)
    snps <- gwas_pos[which(gwas_pos[,"rsid"] %in% trans_eqtls_circle$rsid),]
  
    # 4.2.2: Probes: Can use probe GenomicRanges file (see earlier)
    probe_gr <- as.data.frame(probe_coord[trans_eqtls_circle$probe])[,2:4]
    probe_gr$seqnames <- gsub("chr", "", probe_gr$seqnames)
    probe_gr$p_eqtl <- trans_eqtls_circle$p_eqtl
    probe_gr$qqid <- trans_eqtls[["trans_eqtls.id"]]$qqid
    colnames(probe_gr) <- c("probe", "probechr", "probebp", "p_eqtl", "qqid")
  
    # 4.2.3: Combine the SNP and probe coordinates together, per trans-eQTL
    to_circleplot <- cbind(snps, probe_gr)
  
    # 4.2.4: Optional: Get positions and p-values of other GWS SNPs
    gws <- read.table("gws_snps.txt") # a single column of GWS snps
    gws_snps_df <- gwas_pos[which(gwas_pos[,"rsid"] %in% gws[,1]),]

  # 4.3: Function - create circle plot

    # May choose to set working directory
    # setwd("../outdir")

    # 4.3.1: If a dataframe of gws_snps are not input:
    p <- circle_plot(to_circleplot, gws_snps = NA, gen_name = "hg19",
                id_order = qq_out[["id_ordered.row"]], 
                hg_chr = seq_name,  hg_len = seq_length)

    # 4.3.2: If a dataframe of gws_snps are input (gws_snps_df):
    p <- circle_plot(to_circleplot, gws_snps = gws_snps_df, gen_name = "hg19", 
                id_order = qq_out[["id_ordered.row"]], legend = "both",
                hg_chr = seq_name,  hg_len = seq_length)

  # 4.4: Save plots + formatting
  
    # 4.4.1: If you do want a legend:
    p + ggtitle("trans-eQTL Circle Plot") +
      theme(plot.title = element_text(lineheight = 10, face = "bold"))
    ggsave(filename = "circleplot_legend.pdf", 
      plot = last_plot(), width = 8.2, height = 8.2)

    # 4.4.2: If you don't want a legend:
    # p <- p + theme(legend.position="none")  
    # p + ggtitle("trans-eQTL Circle Plot") +
    #   theme(plot.title = element_text(lineheight = 10, face = "bold"))
    # ggsave(filename = "circleplot_nolegend.pdf", 
    #   plot = last_plot(), width = 7.2, height = 7.2)

    
#==============================================================================#
#==============================================================================#
#                                                                              #
#                       EXAMPLE SNP RESAMPLING ANALYSIS                        #
#                                                                              #
# This is an example script showing a workflow for a SNP trans-eQTLs analysis  #
# as discussed in:                                                             #
#                                                                              #
#    "Independent trans-eQTL identified in whole blood have                    # 
#         limited influence on complex disease biology "                       #
#                                                                              #
# This workflow requires functions that are available at <http://_______>.     #
#                                                                              #
# Difference between discovery and resampling analyses:                        #
#   This analysis requires a list ("eqtlfile.ls") of extracted eQTLs,          #
#   each entry containing a dataframe representing one resampling              #
#   (dataframe specifications are the same in discovery analysis "eqtlfile" )  #
#                                                                              #
# It follows 5 general steps:                                                  #
#                                                                              #
#   1) Obtaining trans-eQTLs and cis-eQTLs                                     #
#       Function: trans_cister_eqtls(..., resamp = TRUE)                       #
#                                                                              #
#   2) Extracting SNPs within the trans-probe region                           #
#       Function: transprobe(..., resamp = TRUE)                               #
#                                                                              #
#      ... (Optional LD clumping) ...                                          #
#                                                                              #
#   3) Producing a QQ plot, showing the GWAS p-values for transprobe SNPs      #
#       Function: pQQ2.fun.lite() - "base" with discovery transprobe SNPs      #
#       Function: points_pQQ2.fun() - add points per resampling                #
#                                                                              #
#   4) Find the discovery transprobe SNPs' lambda compare to resampling dist   #
#       Function: lambda_snps() - calculates lamda, per resampling             #
#       Function: arrow_hist() - calculates discovery lambda + histograms      #
#                                                                              #
#   5) Create summary tables                                                   #
#                                                                              #
#                                                                              #
#                           Author: Chloe X Yap                                #
#                                                                              #
#==============================================================================#
#==============================================================================#

#-------------------------------------------------------------------------------
# 0.0 Load things
#-------------------------------------------------------------------------------

  # 0.1 Data: eQTL file
  # N.B. The columns that will be used within these functions can't be factors
  # (ie. probe, gene, chr)
  eqtlfile.ls <- readRDS("BMI_eqtl_perm.RData")

  # 0.2 Data: GWAS SNP summary statistics file
  gwas_pos <- read.table("BMI_SNP_posfile.txt", header = T, as.is = T)    
  gwas_col <- c(2, 10, 1, 4) #order rsid, p-value, chr, bp
  gwas_pos <- gwas_pos[,gwas_col]
  gwas_pos <- gwas_pos[order(gwas_pos[,3], gwas_pos[,4]),]
  names(gwas_pos) <- c("rsid", "p_snp", "snpchr", "snpbp")

  # 0.3 Data: Load FDb.UCSC.snp137common.hg19 package
  # source("http://bioconductor.org/biocLite.R")
  # biocLite("FDb.UCSC.snp137common.hg19")
  print("... Calling features of UCSC common SNP package. Takes time!")
  snp137common <- features(FDb.UCSC.snp137common.hg19)
  # For speed, filter to snp137common entries that match GWAS
  # as enrichment tests are limited by SNPs in GWAS dataset
  snp137common <- snp137common[which(names(snp137common) %in% gwas_pos[,1])] 

  # 0.4 Data: Load GRangesList of probe coordinates
  probe_coord <- readRDS("IlluminaHT12v4_GRlistProbesFromBSGS_CHR6HAP.rds")

  # 0.5 Data: Transprobe SNPs from the discovery analysis
  disc_snps <- read.delim("discovery_qqplotSNPs.txt", header = T, as.is = T)

#-------------------------------------------------------------------------------
# 1.0 Obtain trans-eQTLs and cis-eQTLs
#-------------------------------------------------------------------------------

  # 1.1: Function to convert all factors to characters
  out1.ls <- lapply(eqtlfile.ls, trans_cister_eqtls, 
		resamp = TRUE, keep_allsig = FALSE)

  # 1.2: To keep all trans-eQTLs, extract them from out1.ls
  trans_eqtls.ls <- lapply(1:length(out1.ls), 
		function(x) out1.ls[[x]][["trans_eqtls"]])

#-------------------------------------------------------------------------------
# 2.0 Find SNPs within the trans-gene region
#-------------------------------------------------------------------------------

  # 2.1: Load libraries
  require(GenomicRanges)
  require(plyr)

  # 2.2: Extract trans-probe region SNPs. Note: this takes a while!
  # If you want to perform LD clumping on the extracted SNPs, and the .bam files
  # are by chromosome, we recommend having a column for chromosome in the "gwas"
  # file to accelerate this process. 
  # An example ld clumping script is provided (in R, calling PLINK)
  qqsnp.ls <- lapply(trans_eqtls.ls, transprobe, gwas = gwas_pos, 
    snp_gr = snp137common, probe_gr = probe_coord, resamp = TRUE,  
    trans_gene = 5e4, out_dir = "../Output/", name = "bmi")

  # Alternatively for speed if there are a lot of probes that are repeated, 
  # you could aggregate all of the unique probes and input them instead of 
  # trans_eqtls.ls (here, resamp = FALSE, as not in list form). 
  # Then, you could reconstitute probes back to their respective resampling.

#-------------------------------------------------------------------------------
# (Optional: LD clumping)
#-------------------------------------------------------------------------------

  # First, get list elements where there *are* SNPs 
  # (ie. don't have "NA", 0.5)
  notempty <- which(unlist(lapply(1:length(qqsnp.ls), 
      function(x) isTRUE(qqsnp.ls[[x]][1,1] == NA))) == FALSE)

  # An example clumping function is provided at <http://____>

  # You would then need to re-read in the clumped table
  # Below is a conceptual outline:
  # Eg. clump[notempty] <- lapply(qqsnps.ls[notempty], ld_clump_fun)
  #     qqsnps.ls[notempty] <- clump[notempty]
  #     setwd("outdir") # if you moved directory for clumping ouputs
  #                     # move back to outdir to ensure continuity of pipeline

  # N.B. If you decide to LD clumping SNP resampling transprobe regions,
  # then you should also LD clump the discovery transprobe regions so that the
  # distributions are comaprable on the QQ-plot.

#-------------------------------------------------------------------------------
# 3.0 QQ-plot generation
#-------------------------------------------------------------------------------

  # 3.1: Create base plot with discovery SNPs

  # 3.1.1: Plot the transprobe SNPs from the discovery analysis
  png(paste(disease_name, "_snpresamp.png", sep = ""), res = 300, 
    height = 29.7, width = 21, units = "cm")

  # 3.1.2: Call plot
  pQQ2.fun.lite(disc_snps$p_snp, main = paste(disease_name, "QQ Plot", sep = "_"), 
    limx = c(0, log10(max(
                          max(unlist(lapply(qqsnp.ls,nrow))),
                          length(disc_snps$p_snp)))), 
    lim = c(0,-log10(min(
                          min(unlist(lapply(lapply(qqsnp.ls,`[`,i="p_snp"), min))),
                          min((disc_snps$p_snp))))))

  # 3.2: lapply() points from resampling on
  lapply(lapply(qqsnp.ls, function(x) x[,2]), points_pQQ2.fun)

#-------------------------------------------------------------------------------
# 4.0: Lambda calculation and histograms
#-------------------------------------------------------------------------------

  # 4.1: Calculate lambda (mean and median) per SNP resampling
  lambda_list <- lapply(qqsnp.ls, lambda_snps)
  lambda_out_med <- unlist(lapply(lambda_list, function(x) x[1]))
  lambda_out_mea <- unlist(lapply(lambda_list, function(x) x[2]))
  lambda_out <- as.data.frame(cbind(lambda_out_med, lambda_out_mea))
  colnames(lambda_out) <- c("median_lambda", "mean_lambda")

  # 4.2: Plot lambda histograms, one each for mean and median lambda
  actual_lambda <- arrow_hist(disc_snps$p_snp, lambda_out, 
    resamp_n = 1000, out_dir = "")
  colnames(actual_lambda) <- c("median_lambda", "mean_lambda")

#-------------------------------------------------------------------------------
# 5.0: Output summary tables
#-------------------------------------------------------------------------------

  # 5.1: Create summary statistics table
  # From trans_cister_eqtls() output
  trans_eqtls_num <- unlist(lapply(1:length(out1.ls), 
  		function(x) out1.ls[[x]][["trans_eqtls_num"]]))
  trans_esnps_num <- unlist(lapply(1:length(out1.ls), 
  		function(x) out1.ls[[x]][["trans_esnps_num"]]))
  cister_eqtls_num <- unlist(lapply(1:length(out1.ls), 
  		function(x) out1.ls[[x]][["cister_eqtls_num"]]))	
  cister_esnps_num <- unlist(lapply(1:length(out1.ls), 
  		function(x) out1.ls[[x]][["cister_esnps_num"]]))

  # Make summary table	
  summ_tab <- cbind(1:nrow(lambda_out), trans_eqtls_num, trans_esnps_num, 
  	cister_eqtls_num, cister_esnps_num, lambda_out)
  colnames(summ_tab) <- c("resamp", "trans_eqtls_num", "trans_esnps_num",
    "cister_eqtls_num", "cister_esnps_num", "median_lambda", "mean_lambda")

  # Print the lambda rank of discovery vs distribution
  lambda_tab <- data.frame(length(which(summ_tab$median_lambda > actual_lambda[1,1])),
    length(which(summ_tab$mean_lambda > actual_lambda[1,2])))
  colnames(lambda_tab) <- c("median_lambda_rank", "mean_lambda_rank")
  write.table(lambda_tab, paste(out, "lambda_rank.txt", sep = ""),
    col.names = T, row.names = F, sep = "\t", quote = F)

  # Add discovery analysis lambda to the summary table
  disc <- unlist(c(0, rep("NA", 4), actual_lambda))
  summ_tab <- rbind(disc, summ_tab)

  # Write table
  write.table(summ_tab, "snp_resamp_summ_table.txt",
    col.names = T, row.names = F, sep = "\t", quote = F)

  # 5.2: Create table containing all trans-eQTLs
  saveRDS(trans_eqtls.ls, "resamp_trans_eqtls.rds")

#==============================================================================#
#==============================================================================#
#                                                                              #
#                      EXAMPLE PROBE RESAMPLING ANALYSIS                       #
#                                                                              #
# This is an example script showing a workflow for a probe trans-eQTLs         #
# analysis as discussed in:                                                    #
#                                                                              #
#    "Independent trans-eQTL identified in whole blood have limited            # 
#             influence on complex disease biology "                           #
#                                                                              #
# This workflow requires functions that are available at <http://_______>.     #
#                                                                              #
# Difference between discovery and resampling analyses:                        #
#   This analysis requires the probe coordinates file that is a GRangesList    #
#   object, one list entry per probe, and the entry named by the probe ID      #
#                                                                              #
# It follows 5 general steps:                                                  #
#                                                                              #
#   1) Generate the probe resamplings                                          #
#       (code is in this file)                                                 #
#                                                                              #
#   2) Extracting SNPs within the trans-probe region                           #
#       Function: transprobe(..., resamp = TRUE)                               #
#                                                                              #
#      ... (Optional LD clumping) ...                                          #
#                                                                              #
#   3) Producing a QQ plot, showing the GWAS p-values for transprobe SNPs      #
#       Function: pQQ2.fun.lite() - "base" with discovery transprobe SNPs      #
#       Function: points_pQQ2.fun() - add points per resampling                #
#                                                                              #
#   4) Find the discovery transprobe SNPs' lambda compare to resampling dist   #
#       Function: lambda_snps() - calculates lamda, per resampling             #
#       Function: arrow_hist() - calculates discovery lambda + histograms      #
#                                                                              #
#   5) Create summary tables                                                   #
#                                                                              #
#                                                                              #
#                           Author: Chloe X Yap                                #
#                                                                              #
#==============================================================================#
#==============================================================================#
#-------------------------------------------------------------------------------
# 0.0 Load things
#-------------------------------------------------------------------------------

  # 0.1 Data: GWAS SNP summary statistics file
  gwas_pos <- read.table("BMI_SNP_posfile.txt", header = T, as.is = T)    
  gwas_col <- c(2, 10, 1, 4) #order rsid, p-value, chr, bp
  gwas_pos <- gwas_pos[,gwas_col]
  gwas_pos <- gwas_pos[order(gwas_pos[,3], gwas_pos[,4]),]
  names(gwas_pos) <- c("rsid", "p_snp", "snpchr", "snpbp")

  # 0.2 Data: Load FDb.UCSC.snp137common.hg19 package
  # source("http://bioconductor.org/biocLite.R")
  # biocLite("FDb.UCSC.snp137common.hg19")
  print("... Calling features of UCSC common SNP package. Takes time!")
  # snp137common <- features(FDb.UCSC.snp137common.hg19)
  # For speed, filter to snp137common entries that match GWAS
  # as enrichment tests are limited by SNPs in GWAS dataset
  # snp137common <- snp137common[which(names(snp137common) %in% gwas_pos[,1])] 
  # Or just load in the provided .rds file for speed
  snp137common <- readRDS("snp137common.rds")

  # 0.3 Data: Load GRangesList of probe coordinates
  probe_coord <- readRDS("IlluminaHT12v4_GRlistProbesFromBSGS_CHR6HAP.rds")

  # 0.4 Data: Transprobe SNPs from the discovery analysis
  disc_snps <- read.delim("discovery_qqplotSNPs.txt", header = T, as.is = T)

#-------------------------------------------------------------------------------
# 1.0 Probe resampling
#-------------------------------------------------------------------------------

  # 1.0 Choose settings

  # 1.0.1 Choose number of resamplings
  resamplings <- 1000

  # 1.0.2 Choose number of probes to sample = number of sig trans_eqtls
  # For BMI, this was 2
  n_hits <- 2

  # 1.2 Sample probes, without replacement
  # Do resampling, put probes in list, one entry per resampling
  probe_resamp.ls <- replicate(resamplings, 
    list(as.data.frame(sample(names(probe_coord), n_hits, replace = F))))
  # Change col names of each entry to "probe", so transprobe_qq can be used
  probe_resamp.ls <- lapply(probe_resamp.ls, setNames, "probe")

#-------------------------------------------------------------------------------
# 2.0 Find SNPs within the trans-gene region
#-------------------------------------------------------------------------------

  # 2.1: Load libraries
  require(GenomicRanges)
  require(plyr)

  # 2.2: Extract trans-probe region SNPs. Note: this takes a while!
  # If you want to perform LD clumping on the extracted SNPs, and the .bam files
  # are by chromosome, we recommend having a column for chromosome in the "gwas"
  # file to accelerate this process. 
  # An example ld clumping script is provided (in R, calling PLINK)
  qqsnp.ls <- lapply(probe_resamp.ls, transprobe, gwas = gwas_pos, 
    snp_gr = snp137common, probe_gr = probe_coord, resamp = TRUE,  
    trans_gene = 5e4, out_dir = "Output/", name = "bmi")

  # Alternatively for speed if there are a lot of probes that are repeated, 
  # you could aggregate all of the unique probes and input them instead of 
  # trans_eqtls.ls (here, resamp = FALSE, as not in list form). 
  # Then, you could reconstitute probes back to their respective resampling.

  # 2.3 Replace resamplings with no SNPs in the transgene region/s
  # Look for data.frame("NA", 0.5); proxy is the number of rows
  rows4missing <- unlist(lapply(qqsnp.ls, nrow))
  missing <- which(rows4missing == 1)

  # while statement: continue until no data.frame("NA", 0.5) in the list
  while (length(missing) > 0) {

  probe_resamp.ls[missing] <- replicate(length(missing), 
    list(as.data.frame(sample(names(probe_coord), n_hits, replace = F))))
  # Change col names of each entry to "probe", so transprobe_qq can be used
  probe_resamp.ls <- lapply(probe_resamp.ls, setNames, "probe")

  # Re-extract trans-probe region SNPs
  qqsnp.ls[missing] <- lapply(probe_resamp.ls[missing], transprobe, gwas = gwas_pos, 
    snp_gr = snp137common, probe_gr = probe_coord, resamp = TRUE,  
    trans_gene = 5e4, out_dir = "Output/", name = "bmi")

  rows4missing <- unlist(lapply(qqsnp.ls, nrow))
  missing <- which(rows4missing == 1)

  }

#-------------------------------------------------------------------------------
# 3.0 QQ-plot generation
#-------------------------------------------------------------------------------

  # 3.1: Create base plot with discovery SNPs

  # 3.1.1: Plot the transprobe SNPs from the discovery analysis
  png(paste(disease_name, "_proberesamp.png", sep = ""), res = 300, 
    height = 29.7, width = 21, units = "cm")

  # 3.1.2: Call plot
  pQQ2.fun.lite(disc_snps$p_snp, main = paste(disease_name, "QQ Plot", sep = "_"), 
    limx = c(0, log10(max(
                          max(unlist(lapply(qqsnp.ls,nrow))),
                          length(disc_snps$p_snp)))), 
    lim = c(0,-log10(min(
                          min(unlist(lapply(lapply(qqsnp.ls,`[`,i="p_snp"), min))),
                          min((disc_snps$p_snp))))))

  # 3.2: lapply() points from resampling on
  lapply(lapply(qqsnp.ls, function(x) x[,2]), points_pQQ2.fun)

  dev.off()

#-------------------------------------------------------------------------------
# 4.0: Lambda calculation and histograms
#-------------------------------------------------------------------------------

  # 4.1: Calculate lambda (mean and median) per SNP resampling
  lambda_list <- lapply(qqsnp.ls, lambda_snps)
  lambda_out_med <- unlist(lapply(lambda_list, function(x) x[1]))
  lambda_out_mea <- unlist(lapply(lambda_list, function(x) x[2]))
  lambda_out <- as.data.frame(cbind(lambda_out_med, lambda_out_mea))
  colnames(lambda_out) <- c("median_lambda", "mean_lambda")

  # 4.2: Plot lambda histograms, one each for mean and median lambda
  actual_lambda <- arrow_hist(disc_snps$p_snp, lambda_out, 
    resamp_n = 1000, out_dir = "")
  colnames(actual_lambda) <- c("median_lambda", "mean_lambda")

#-------------------------------------------------------------------------------
# 5.0: Output summary tables
#-------------------------------------------------------------------------------

  # 5.1: Create summary statistics table
  # Make summary table  
  summ_tab <- cbind(1:nrow(lambda_out), lambda_out)
  colnames(summ_tab) <- c("resamp", "median_lambda", "mean_lambda")

  # Print the lambda rank of discovery vs distribution
  lambda_tab <- data.frame(resamplings-(length(which(summ_tab$median_lambda > actual_lambda[1,1]))),
    resamplings-(length(which(summ_tab$mean_lambda > actual_lambda[1,2]))))
  colnames(lambda_tab) <- c("median_lambda_rank", "mean_lambda_rank")
  write.table(lambda_tab, paste(out, "lambda_rank.txt", sep = ""),
    col.names = T, row.names = F, sep = "\t", quote = F)

  # Add discovery analysis lambda to the summary table
  disc <- unlist(c(0, actual_lambda))
  summ_tab <- rbind(disc, summ_tab)

  # Write table
  write.table(summ_tab, "probe_resamp_summ_table.txt",
    col.names = T, row.names = F, sep = "\t", quote = F)