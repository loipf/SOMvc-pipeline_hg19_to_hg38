
############################################
### check maftools

setwd(dirname(rstudioapi::getSourceEditorContext()$path))   ### only in RStudio

source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/small_functions.R")

pacman::p_load(data.table, GenomicRanges, dplyr)



############################################
### download 1-based gff3 and get 0-based bed file

# ### on ensembl - only exons
# wget -qO- ftp://ftp.ensembl.org/pub/release-101/gff3/homo_sapiens/Homo_sapiens.GRCh38.101.gff3.gz \
# | gunzip --stdout - \
# | awk '$3 == "exon"' - \
# | convert2bed -i gff --attribute-key="gene_name" - \
# > Homo_sapiens.GRCh38.exons.101.bed
s
# ### on ensembl - all cds
# https://github.com/loipf/SOMvc-pipeline/blob/main/scripts/data_acquisition.sh


############################################
### combine bed segments

bed <- '/home/stefanloipfinger/Desktop/other/bed_file_test/Homo_sapiens.GRCh38.exons.101.bed'
# bed <- '/home/stefanloipfinger/Desktop/other/bed_file_test/Homo_sapiens.GRCh38.exons.101.squeezed.bed'

### filter for transcript_support_level==1 like in martens?
bed.df <- fread(bed) %>% as.data.frame()
bed.df = bed.df[,c(1:3,6)]

colnames(bed.df) = c("chromosome","start","end","strand")

# bed.df = subset(bed.df, chromosome %in% c("X","Y","MT",1:22) )

bed.gr <- makeGRangesFromDataFrame(
  df = bed.df, keep.extra.columns = F,
  ignore.strand = TRUE,
  starts.in.df.are.0based = TRUE
)

# Total bases before combining:
width(bed.gr) %>% sum()

# Combinine overlapping regions:
bed.gr.reduced <- reduce(bed.gr)
exome_size = width(bed.gr.reduced) %>% sum()

# 1266076039  ### ensembl cds
#  133277685  ### gencode
#  147690727  ### ensembl 101
#   28711682  ### TSL==1


bed.gr.reduced_tidy = keepStandardChromosomes(bed.gr.reduced, pruning.mode = 'tidy')
width(bed.gr.reduced_tidy) %>% sum()  ## 1265712403

bed_reduced_tidy_bed = data.frame(bed.gr.reduced_tidy)

### remove pieces with length 1
bed_reduced_tidy_bed = bed_reduced_tidy_bed[bed_reduced_tidy_bed$width>1,]

fwrite(bed_reduced_tidy_bed, paste0(basename_core(bed, with_path = T),'.squeezed.bed'), sep = '\t', col.names = F)



