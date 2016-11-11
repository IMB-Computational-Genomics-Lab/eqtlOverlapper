#============================================================================
# Example LD clumping function (in R) to lapply() over
#
# Variables:
#   qqplotSNPs:   dataframe of snpchr, rsid, p_snp (columns in that order)
#   clump_dir:    directory to output the clumping results
#   out_dir:      directory for overall analysis pipeline output 
#   plink_dir:    directory of PLINK software
#   bfile_dir:    directory of bfiles
#   p1p2r2kb:     vector of ld clumping settings in order of:
#                 clump-p1, clump-p2, clump-r2, clump-kb
#                 Default: c(0.1, 0.1, 0.1, 1000)
#   pre_exclude:  to reduce number of SNPs to clump, can set a p-value
#                 pre-exclusion threshold
#============================================================================

ld_clump_fun <- function(qqplotSNPs, 
                        clump_dir = clump_dir, 
                        plink_dir = plink_dir,
                        bfile_dir = bfile_dir,
                        p1p2r2kb = c(0.1, 0.1, 0.1, 1000),
                        pre_exclude = 0.5) {

  # Extract settings
  clump_p1 <- p1p2r2kb[1]
  clump_p2 <- p1p2r2kb[2]
  clump_r2 <- p1p2r2kb[3]
  clump_kb <- p1p2r2kb[4]

  #1.0 LD clumping
    setwd(clump_dir) #write qqplotSNPs2plink to ld_clumping direcotry

    #Reduce size of qqplotSNPs2plink file to speed up clumping
    qqplotSNPs2plink <- qqplotSNPs[which(qqplotSNPs$p_snp < pre_exclude),] 
        #need chr for ld clump bychr.txt
        #pre_exclude ok as the clump-p1 and p2 values are less than this anyway
    colnames(qqplotSNPs2plink) <- c("Chr", "P", "SNP")
    write.table(qqplotSNPs2plink, file = "qqplotSNPs2plink.txt", sep = "\t", 
      row.names = F, col.names = T, quote = F)

    #Call LD clumping
    #Make foo - column of chr in resampling
      system("touch clump_mergebychr.txt; sed '1d' qqplotSNPs2plink.txt | cut -f 1 | sort | uniq > foo")

    #Use foo to split up to bychr.txt --> name bychr$n.txt
      system("
        while read line
        do 
          echo \"$line\" 
          n=\"$line\"

          awk -v var=\"$n\" '$1 == var' qqplotSNPs2plink.txt > bychr$n.txt
          awk 'NR==1' qqplotSNPs2plink.txt | cat - bychr$n.txt > temp && mv temp bychr$n.txt
        done <foo
      ")

    #Then read over foo to read correct bychr$n.txt file --> use PLINK on that
      system(paste("
      while read line
      do 
        echo \"$line\" 
        n=\"$line\"
      
        ", plink_dir, 
        " --allow-no-sex --bfile ", bfile_dir, 
        " --clump ", clump_dir, "bychr$n.txt", 
        " --clump-p1 ", clump_p1, " --clump-p2 ", clump_p2, " --clump-r2 ", clump_r2, " --clump-kb ", clump_kb,
        " --out ", clump_dir, "probes_clumped$n

      done <foo
      ", sep = ""))
      
    #Then merge the files
      system("
      while read line
      do 
        echo \"$line\" 
        n=\"$line\"
      
        awk -v var1=\"$n\" '$1 == var1' probes_clumped$n.clumped > probes_clumpednoNA$n.clumped
        cat probes_clumpednoNA$n.clumped >> clump_mergebychr.txt
      done <foo 
      ")
      
    #Finish off
      system("awk '{print $3 \" \" $5}' clump_mergebychr.txt > clumped.txt; times")

    #Desired output = clump_FINAL.txt
    ld_clump <- read.table(paste(clump_dir, "clumped.txt", sep = ""), header = F)

  #2.0 Reconstitute SNPs
    #Gives all rows (even if repeated) where the ld_clump SNPids are found in qqplotSNPs
    qqplotSNPs <- qqplotSNPs[which(qqplotSNPs$rsid %in% ld_clump[,1]),]
      #Overwrite qqplotSNPs so no pipeline hiccups
      #Output: dataframe like qqplotSNPs, only with sentinels

return(qqplotSNPs)
}