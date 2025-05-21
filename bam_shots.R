#!/usr/bin/env

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

.libPaths('/lustre/scratch126/tol/teams/lawniczak/users/jr35/sware/personal_R_libs')

##load libraries
library(Biostrings)
library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(Gviz)
library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(swissknife)


##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

chr <- args[1]
pos <- args[2]
gn <- args[3]
pf_id <- args[4]
windw<-args[5]
i_dir<-args[6]
o_dir<-args[7]

##read  in table with plotting details position and gene names
#pos_table <- read_csv(strain_top_snps)

##read  in table with plotting details position and gene names
#pos_table <- read_csv(strain_top_snps)

##Plotting script
gn_name <- paste(gn, pf_id, sep = ': ')
start_pos <- as.numeric(pos)-as.numeric(windw)
end_pos <- as.numeric(pos)+as.numeric(windw)

##get fasta
fasta_ref=readDNAStringSet('/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7_Genome.fasta')
fasta_ref@ranges@NAMES<-str_extract(fasta_ref@ranges@NAMES, 'Pf.*v3')

##get gtf granges
pf_gr <- prepareGTF('/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7.gtf')

options(ucscChromosomeNames=FALSE)

##construct general tracks
grtrack <- GeneRegionTrack(pf_gr, chromosome = chr, transcriptAnnotation = "gene", rotation.title = 0, transcriptAnnotation = 'gene', cex.title =  1, name = 'Gene', fontcolor.title = 'black')
gtrack <- GenomeAxisTrack(labelPos="beside", scale = 10, cex=1.1, col='black')
sTrack<-SequenceTrack(fasta_ref,cex=1.2, genome = 'Pf52', chromosome = chr)


altracks <- list.files(i_dir,
   pattern = "_sorted.bam$", ##more specific regex option to limit to stages and strains of interest "^(asex|fh|fl|m|Strain|Doublet).*_sorted.bam$"
   full.names = T,
   recursive = T,
   include.dirs = T)%>%
set_names(nm = (str_remove_all(.,paste0(i_dir,'/|_filtered_sorted.bam|_strain|')) %>%
#        str_replace(., '_Strain_', '\nStrain\n') %>%
#        str_replace(., '_', ' ') %>%
	str_replace(., '_G', '\nG') %>%
	str_replace(., '_', '') %>%
	tools::file_path_sans_ext()))

altracks_order <- names(altracks)[order(as.numeric(gsub("[^0-9]", "", names(altracks))))]
altracks <- altracks[altracks_order]

at <- list(altracks, altracks_order)%>%
pmap(~ ..1 %>% AlignmentsTrack(.,
        isPaired = TRUE,
        stacking = "squish",
        name = ..2,
        cex.title=1.22,
        cex.axis = 1.2,
        fontcolor.title = 'black',
        background.title="white",
        col.axis = 'black',
        col.border.title="lightgray"))

##Comment out highlighting with red bar
# htrack <- HighlightTrack(trackList = c(at, grtrack,gtrack,sTrack), start = pos, end = pos,chromosome = chr)

pdf(paste0(o_dir, str_replace(gn_name, ': ', '_'),'.pdf'))

##Comment out highlighting with red bar
#plotTracks(c(htrack),
plotTracks(c(at, grtrack,gtrack,sTrack),
       sizes=c(rep(5, length(altracks)), c(1,1.5,1)),
       # sizes=c(rep(5, length(altracks)), c(1,3)),
       from = start_pos,
       to = end_pos,
       chromosome = chr,
       type = "coverage",
       main = gn_name
      )
dev.off()



