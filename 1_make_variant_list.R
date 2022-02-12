options(width=2000)

suppressMessages(require(pQTLtools))
suppressMessages(require(GenomicRanges))
suppressMessages(require(BiocGenerics))
require(openxlsx)
suppressMessages(require(dplyr))
suppressMessages(require(gap))
require(tidyr)
hbf_dir <- "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HbF"
hbf_xlsx <- "hbf_GWAS_top_snps.xlsx"
if(!dir.exists(file.path(hbf_dir,"work"))) dir.create(file.path(hbf_dir,"work"))
hbf_hits <- read.xlsx(file.path(hbf_dir,hbf_xlsx),sheet=1,colNames=TRUE,startRow=1,skipEmptyRows=TRUE) %>%
            rename(gene="Candidate.gene(s)",b=X9,p.value="p-value") %>%
            mutate(p=pvalue(b/SE),seqnames=paste0("chr",CHR),start=BP,end=BP)
# only as a reference, chr2:60718043_G_T rs1427407 indicates that it is in GRCh37/hg19 format
hbf_hits

liftRegion <- function(x,chain)
{
  gr <- as(x,"GRanges")
  seqlevelsStyle(gr) <- "UCSC"
  gr38 <- rtracklayer::liftOver(gr, chain)
  chr38 <- as.character(seqnames(gr38))
  start38 <- as.integer(start(gr38))
  end38 <- as.integer(end(gr38))
  invisible(data.frame(x,chr38,start38,end38))
}

f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
chain <- rtracklayer::import.chain(f)

hbf_hits_lifted <- liftRegion(select(hbf_hits,seqnames,start,end,REF,ALT),chain) %>%
                   mutate(snpid=chr_pos_a1_a2(seqnames,start,REF,ALT,prefix=""),
                          snpid38=if_else(chr38=="NA","NA",chr_pos_a1_a2(chr38,start38,REF,ALT,prefix=""))) %>%
                   left_join(hbf_hits) %>%
                   select(Locus,snpid,snpid38,rs.ID,gene,p) %>%
                   separate_rows(gene, sep=",", convert = TRUE) %>%
                   mutate(gene=gsub(" ","",gene))
hbf_hits_na <- left_join(hbf_hits_lifted,select(pQTLtools::SomaLogic160410,-chr,-start,-end),by=c('gene'='entGene')) %>%
               select(snpid,rs.ID,SOMAMER_ID,UniProt,Target,ensGene,extGene)
hbf_hits_long <- left_join(hbf_hits_lifted,data.frame(pQTLtools::hg19Tables),by=c('gene'='geneName')) %>%
                 select(snpid,names(hbf_hits_lifted),acc,uniprotName,geneSynonyms,hgncSym,ensGene) %>%
                 mutate(prot=gsub("_HUMAN","",uniprotName)) %>%
                 select(-uniprotName)
# a more useable form
data.frame(hbf_hits_long)

write.table(select(hbf_hits,-seqnames,-start,-end),file=file.path(hbf_dir,"work","hbf_GWAS_top_snps.txt"),row.names=FALSE,quote=FALSE,sep="\t")
write.table(hbf_hits_long,file=file.path(hbf_dir,"work","hbf_GWAS_top_snps_long.txt"),row.names=FALSE,quote=FALSE,sep="\t")

hbf_hits_region <- liftRegion(select(hbf_hits,seqnames,start,end,REF,ALT),chain) %>%
                   mutate(snpid=chr_pos_a1_a2(seqnames,start,REF,ALT,prefix=""),
                          snpid38=if_else(chr38=="NA","NA",chr_pos_a1_a2(chr38,start38,REF,ALT,prefix=""))) %>%
                   left_join(hbf_hits) %>%
                   mutate(pos38=start38) %>%
                   arrange(CHR,BP) %>%
                   select(CHR,BP,pos38,rs.ID,,REF,ALT,b,SE,p,Locus,snpid,snpid38,gene)
write.table(hbf_hits_region,file=file.path(hbf_dir,"work","hbf_hits.txt"),row.names=FALSE,quote=FALSE,sep="\t")
