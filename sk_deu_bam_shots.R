#!/usr/bin/env

##load personal library folder
.libPaths('/lustre/scratch126/tol/teams/lawniczak/users/jr35/sware/personal_R_libs')

##load libraries
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(Gviz)
library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(swissknife)
library(MetBrewer)


##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

stg_comp_tbl_path <- args[1]
stgs_comparison <- args[2]
i_dir<-args[3]
o_dir<-args[4]

stgs_comparison <- 'top_gns'
i_dir<-'/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/V3_4TP/sc_velo/data/processed/PfDB52e/all_days/stage_filtered_bam/lr_iso/bed_bwig/'
o_dir<-'/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/V3_4TP/sc_velo/data/processed/PfDB52e/all_days/stage_filtered_bam/lr_iso/bed_bwig/bw_plots'


##getting granges object from gtf
granges_gtf <- prepareGTF('/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7.gtf')

##creating gene id, gene name combined with exon variables from csv table with all the information
stg_comp_tbl <- read_csv(stg_comp_tbl_path)

#gn <- stg_comp_tbl %>%
#	filter(stg_comp == stgs_comparison & !str_detect(gene_id, '\\+')) %>%
#	pull(gene_id)
#gn_title <- stg_comp_tbl%>%
#	filter(stg_comp == stgs_comparison & !str_detect(gene_id, '\\+')) %>%
#        mutate( gn_title= paste0(Gene, ': Exons ', str_remove_all(DE_exons, 'E00|E0'))) %>%
#        pull(gn_title)

gn <- ('PF3D7_1438800' 'PF3D7_1403400' 'PF3D7_1122200' 'PF3D7_0926300' 'PF3D7_0609100' 'PF3D7_1006600' 'PF3D7_0306200' 'PF3D7_0916200' 'PF3D7_1247800' 'PF3D7_0518800' 'PF3D7_0216700' 'PF3D7_1145300' 'PF3D7_1469900')
gn_title <- ('PF3D7_1438800' 'PF3D7_1403400' 'PF3D7_1122200' 'PF3D7_0926300' 'PF3D7_0609100' 'PF3D7_1006600' 'PF3D7_0306200' 'PF3D7_0916200' 'PF3D7_1247800' 'PF3D7_0518800' 'PF3D7_0216700' 'PF3D7_1145300' 'PF3D7_1469900')

##Listing bigwig files for the different stages
bigwig_files = list.files(i_dir,
			  pattern = "*.bw$",           
                          full.names = T,           
                          recursive = T,           
                          include.dirs = T) %>%
set_names(nm = (str_remove_all(.,paste0(i_dir,'/|v3|.bw')) %>%
		str_replace(., '_plus', ' +') %>%                
                str_replace(., '_minus', ' -') %>%  
                tools::file_path_sans_ext())) 


branches1 = c("asexuals", "commited", "developing", "branching","ef", "lf","em","lm")

bigwig_files <- bigwig_files[paste(branches1, c('-','+'))]


##Assign variables for naming and coloring tracks
monet_col <- met.brewer('Monet',9,override.order = T)

col_branches1 = c("grey", monet_col[4], monet_col[5], monet_col[6], monet_col[2], monet_col[1], monet_col[8], monet_col[7])

bw_condtn_lab <- names(bigwig_files) %>% set_names(names(bigwig_files))
bw_condtn_cols <- rep(col_branches1, each = 2) %>% set_names(bw_condtn_lab)


##Loop to plot plots per stage pairwise comparison
pdf(paste0(o_dir ,stgs_comparison, '_deu_cov.pdf'))

for (i in 1:length(gn)){
	plotGeneRegion(
		       granges = granges_gtf,            
                       bigwigFiles = bigwig_files,            
                       bigwigCond = bw_condtn_lab,            
                       showgene = gn[i],            
                       geneTrackTitle = '',            
                       colorByStrand = TRUE,            
                       condColors = bw_condtn_cols,            
                       scaleDataTracks = F,            
                       plotTitle = gn_title[i],            
                       cex.title=1.2,            
                       cex.axis = 2,            
                       fontcolor.title = 'black',            
                       background.title="white",            
                       col.axis = 'black',            
                       col.border.title="lightgray",            
                       cex.main = 1.6,
                       frame=T,
                        labelPos="beside",
                        scale = 500, 
                        cex=1.3, 
                        #col='black',
                        showId = F
                                      
                       )
}
dev.off()

##Loop to plot plots per stage pairwise comparison for both the sense and antisense genes

for (strnd in c('minus', 'plus')) {
        
        ##Assign variables for naming and coloring tracks
	bw_hf_files <- bigwig_files[str_detect(bigwig_files, strnd)]
	#bw_hf_files <- bw_hf_files[branches1]
	bw_hf_condtn_lab <- names(bw_hf_files) %>% set_names(names(bw_hf_files))
	bw_hf_condtn_cols <- c(col_branches1) %>% set_names(bw_hf_condtn_lab)
        
	pdf(paste0(o_dir ,stgs_comparison,'_', strnd, '_deu_cov.pdf'))

        for (i in 1:length(gn)){
                plotGeneRegion(granges = granges_gtf,
                                   bigwigFiles = bw_hf_files,
                                   bigwigCond = bw_hf_condtn_lab,
                                   showgene = gn[i],
                                   geneTrackTitle = '',
                                   colorByStrand = TRUE,
                                   condColors = bw_hf_condtn_cols,
                                   scaleDataTracks = F,
                                   plotTitle = gn_title[i],
                                   cex.title=1.5,
                                   cex.axis = 2,
                                   fontcolor.title = 'black',
                                   background.title="white",
                                   col.axis = 'black',
                                   col.border.title="lightgray",
                                   cex.main = 1.6,
                                   frame=T,
                                   labelPos="beside",
                                   scale = 500, 
                                   cex=1.3, 
                                   #col='black',
                                   showId = F
                                   )
        }
        dev.off()
}


