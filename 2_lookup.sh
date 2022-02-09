#!/usr/bin/bash

export HbF=${HOME}/COVID-19/HbF
export pgwas=~/rds/results/public/proteomics

# ARIC, https://sites.cscc.unc.edu/aric/

cat <(awk -v OFS="\t" '{print "Gene",$0}' ${pgwas}/ARIC/glm.linear.hdr) \
    <(
       Rscript -e '
       options(width=200)
       suppressMessages(library(dplyr))
       HbF <- Sys.getenv("HbF")
       pgwas <- Sys.getenv("pgwas")
       snps <- read.delim(file.path(HbF,"work","hbf_GWAS_top_snps_long.txt")) %>%
               left_join(read.delim(file.path(pgwas,"ARIC","seqid.txt")),by=c("acc"="uniprot_id")) %>%
               select(seqid_in_sample,acc,entrezgenesymbol,snpid,snpid38,rs.ID,gene) %>%
               subset(!is.na(seqid_in_sample))
       write.table(select(snps,seqid_in_sample,rs.ID,gene),col.names=FALSE,row.names=FALSE,quote=FALSE)
       ' | \
       parallel -C' ' --env pgwas '
         zgrep -w {2} ${pgwas}/ARIC/EA/{1}.PHENO1.glm.linear.gz | \
         awk -v gene={3} -v OFS="\t" "{print gene,\$0}"
       '
     )

# deCODE

# Fenland

# GTEx, eQTL Catalog

# SCALLOP-CVD1

# ~/COVID-19/HbF/work/hbf_GWAS_top_snps_long.txt
snpid   Locus   snpid38 rs.ID   gene    p       acc     geneSynonyms    hgncSym ensGene prot
chr2:53993037_T_TTCTTTGGTGCTCCATTGATAAGAGCACCC  1       chr2:53765900_T_TTCTTTGGTGCTCCATTGATAAGAGCACCC  rs148978228     CHAC2   3.93e-15        Q8WUX2     >
chr2:60718043_G_T       2       chr2:60490908_G_T       rs1427407       BCL11A  3.01e-161       Q9H165  CTIP1; EVI9; KIAA1809; ZNF856   BCL11A  ENSG0000011>

# ${pgwas}/ARIC/seqid.txt
seqid_in_sample uniprot_id      entrezgenesymbol        chromosome_name transcription_start_site
SeqId_10000_28  P43320  CRYBB2  22      25212564
