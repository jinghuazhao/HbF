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

#Gene    #CHROM  POS     ID      REF     ALT     A1      A1_FREQ TEST    OBS_CT  BETA    SE      T_STAT  P       ERRCODE
#BACH2   6       90236760        rs4707609       T       C       C       0.351865        ADD     7213    0.0144076       0.0174087       0.827607        0.40792     .
#MYB     6       134938596       rs7772031       G       A       G       0.36961 ADD     7213    0.0103618       0.0172359       0.601175        0.547742   .
#MYB     6       135097778       rs7776054       A       G       G       0.269721        ADD     7213    0.0292828       0.0187838       1.55894 0.119054   .
#MYB     6       135129617       rs9376095       T       C       C       0.23118 ADD     7213    0.0152785       0.0197551       0.773394        0.439315   .

# AGES

R --no-save <<END
  HbF <- Sys.getenv("HbF")
  stables <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-27850-z/MediaObjects/41467_2021_27850_MOESM18_ESM.xlsx"
  st16 <- openxlsx::read.xlsx(stables, sheet=1, colNames=TRUE, skipEmptyRows=TRUE, startRow=3)
  write.table(st16,file=file.path(HbF,"work","AGES.txt"),col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
END

cd ${HbF}/work
sed '1d' hbf_GWAS_top_snps_long.txt | cut -f5 | sort | uniq | grep -f - -w AGES.txt | cut -f2
lftp -c mirror https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90086001-GCST90087000/GCST90086395
lftp -c mirror https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90086001-GCST90087000/GCST90086825
lfpt -c mirror https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90087001-GCST90088000/GCST90087207
bgzip GC*/*tsv
cat ${pgwas}/AGES/AGES.hdr
(
cat <<END
Accession rs.ID gene
GCST90086395 rs10128556 HBG1
GCST90086395 rs398114988 HBG1
GCST90087207 rs4707609 BACH2
GCST90086825 rs7772031 MYB
GCST90086825 rs7776054 MYB
GCST90086825 rs9376095 MYB
END
) | \
sed '1d' | \
parallel -C' ' '
  echo {1}-{2}-{3}
  zgrep -w {2} {1}/{1}_buildGRCh37.tsv.gz
'
cd -
#variant_id      p_value chromosome      base_pair_location      effect_allele   other_allele    effect_allele_frequency beta    standard_error
#GCST90086395-rs10128556-HBG1
#rs10128556      0.919811849322134       11      5263683 T       C       0.3493  -0.00174019     0.0172851
#GCST90086395-rs398114988-HBG1
#GCST90086825-rs7772031-MYB
#rs7772031       0.713143617980774       6       135259734       G       A       0.3357  0.00722091      0.0196404
#GCST90086825-rs7776054-MYB
#rs7776054       0.212966846431685       6       135418916       G       A       0.279   0.0257332       0.0206594
#GCST90086825-rs9376095-MYB
#rs9376095       0.536243449393508       6       135450755       C       T       0.2279  -0.0136041      0.0219938
#GCST90087207-rs4707609-BACH2
#rs4707609       0.0717480331291094      6       90946479        C       T       0.352   0.0340157       0.0188864

# deCODE

export deCODE=${pgwas}/deCODE
ls ${deCODE} | tr '_' '\t' | cut -f3 | sort | uniq | grep -f <(sed '1d' hbf_GWAS_top_snps_long.txt | cut -f5 | sort | uniq) -w -
ls ${deCODE} | grep -e MYB_MYB -e BACH2_BACH2 | grep gz
zgrep -w rs4707609 ${deCODE}/11618_83_MYB_MYB.txt.gz
zgrep -w -e rs7772031 -e rs7776054 -e rs9376095 ${deCODE}/12756_3_BACH2_BACH2.txt.gz
# Chrom   Pos     Name    rsids   effectAllele    otherAllele     Beta    Pval    minus_log10_pval        SE      N       ImpMAF
#chr6    90236760        chr6:90236760:C:T       rs4707609       C       T       0.0172  0.037333        1.42791 0.008261        35335   0.34794
#chr6    134938596       chr6:134938596:A:G      rs7772031       A       G       -0.0087 0.296918        0.52736 0.008341        35363   0.32928
#chr6    135097778       chr6:135097778:G:A      rs7776054       G       A       0.0241  0.006332        2.19846 0.008828        35363   0.27214
#chr6    135129617       chr6:135129617:C:T      rs9376095       C       T       0.0020  0.833272        0.07921 0.009501        35363   0.22820

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

# ${pgwas}/AGES/AGES.hdr
# variant_id	p_value	chromosome	base_pair_location	effect_allele	other_allele	effect_allele_frequency	beta	standard_error

# deCODE
# Chrom   Pos     Name    rsids   effectAllele    otherAllele     Beta    Pval    minus_log10_pval        SE      N       ImpMAF
