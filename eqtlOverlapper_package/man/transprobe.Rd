% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transprobe_qq.R
\name{transprobe}
\alias{transprobe}
\title{Extract trans-probe region SNPs}
\usage{
transprobe(trans_eqtl_id, gwas = gwas_pos, snp_gr = snp137common,
  probe_gr = probe_coord, trans_gene = 50000, write = TRUE,
  out_dir = ".", name = "out", resamp = FALSE)
}
\arguments{
\item{trans_eqtl_id}{Dataframe with eqtl probe and qqid;
Eg. list_out[["trans_eQTLs.id"]] from trans_cister_eqtls() output;
Else, colnames as:
1) "Probe"   (colClass: character)
  And if resamp = FALSE: ...
2) "qqid"    (colClass: character) - a unique trans-eQTL identifier.}

\item{gwas}{Dataframe of gwas summary statistics.
gwas column *indices* in order: 
1) rsid        (colClass: character)
2) gwas p-val  (colClass: numeric)
Recommend also having chr and bp columns for output completeness
as columns after rsid and gwas p-val, though this is not essential. 
However, if planning on LD clumping, chr columns are required.}

\item{snp_gr}{GRanges object with SNP position coordinates;
names(snp_gr) corresponds with SNP rsIDs;
Ensure snp_gr genome build matches with probe_gr file.}

\item{probe_gr}{GRangesList object with probe position coordinates; 
names(probe_gr) corresponds with Illumina probe IDs;
Where there are multiple probe coordinates, the first will be used;
Ensure probe_gr genome build matches with snp_gr file.}

\item{trans_gene}{Numeric: extract SNPs around probe start +/- this region; 
Default = 5e4.}

\item{out_dir}{Character: output folder, followed by /; Default is "."}

\item{name}{Character: phenotype name.}

\item{resamp}{TRUE or FALSE}
}
\value{
Dataframe of trans-gene SNPs with columns titled:
 If resamp = FALSE:
   1) "rsid", 2) "p_snp", 3) "probe", 4) "qqid";
   If gwas input had other columns, these will also be incorporated.
 If resamp = TRUE:
   1) "snpchr", 2) "p_snp", 3) "rsid"
   Only actually need "p_snp" and "rsid"; "snpchr" is to 
   make clumping faster, if desired
}
\description{
Extract trans-probe region SNPs
}
\author{
Chloe X Yap
}

