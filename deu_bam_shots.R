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


##getting granges object from gtf
granges_gtf <- prepareGTF('/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7.gtf')

##creating gene id, gene name combined with exon variables from csv table with all the information
stg_comp_tbl <- read_csv(stg_comp_tbl_path)

gn <- stg_comp_tbl %>%
	filter(stg_comp == stgs_comparison & !str_detect(gene_id, '\\+')) %>%
	pull(gene_id)
gn_title <- stg_comp_tbl%>%
	filter(stg_comp == stgs_comparison & !str_detect(gene_id, '\\+')) %>%
        mutate( gn_title= paste0(Gene, ': Exons ', str_remove_all(DE_exons, 'E00|E0'))) %>%
        pull(gn_title)

##Listing bigwig files for the different stages
bigwig_files = list.files(i_dir,
			  pattern = "*.bw$",           
                          full.names = T,           
                          recursive = T,           
                          include.dirs = T) %>%
set_names(nm = (str_remove_all(.,paste0(i_dir,'/|v3|.bw')) %>%
		str_replace(., '_plus', ' +') %>%                
                str_replace(., '_minus', ' -') %>%  
                str_replace(., 'Strain_0 ', 'Strain 0\n') %>%
                str_replace(., 'Strain_1 ', 'Strain 1\n') %>%
                str_replace(., 'Strain_2 ', 'Strain 2\n') %>%
                str_replace(., 'Doublet ', 'Doublet\n') %>%
                tools::file_path_sans_ext()))


##Assign variables for naming and coloring tracks
if (unique(str_detect(bigwig_files, "stage_filtered_bam"))) {
	bw_condtn_lab <- names(bigwig_files) %>% set_names(names(bigwig_files))
    	bw_condtn_cols <- rep(c(met.brewer('Java')[4], met.brewer('Ingres')[7], met.brewer('Ingres')[5],  met.brewer('Nizami')[7]), each = 2) %>% set_names(bw_condtn_lab)
	
} else if (unique(str_detect(bigwig_files, "strain_b4qc_filtered_bam"))) {
	bw_condtn_lab <- names(bigwig_files) %>% set_names(names(bigwig_files))
    	bw_condtn_cols <- rep(c(met.brewer('Austria')[6], met.brewer('Austria')[3], met.brewer('Austria')[5], met.brewer('Austria')[1]), each = 2) %>% set_names(bw_condtn_lab)

} else {
	bw_condtn_lab <- names(bigwig_files) %>% set_names(names(bigwig_files))
    	bw_condtn_cols <- rep(c(met.brewer('Austria')[1], met.brewer('Austria')[3], met.brewer('Austria')[5]), each = 2) %>% set_names(bw_condtn_lab)

}


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
                       cex.axis = 1,            
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
        if (unique(str_detect(bigwig_files, "stage_filtered_bam"))) {
    		bw_hf_files <- bigwig_files[str_detect(bigwig_files, strnd)]
    		bw_hf_condtn_lab <- names(bw_hf_files) %>% set_names(names(bw_hf_files))
    		bw_hf_condtn_cols <- c(met.brewer('Java')[4], met.brewer('Ingres')[7], met.brewer('Ingres')[5],  met.brewer('Nizami')[7]) %>% set_names(bw_hf_condtn_lab)
   
	} else if (unique(str_detect(bigwig_files, "strain_b4qc_filtered_bam"))) {
		
    		bw_hf_files <- bigwig_files[str_detect(bigwig_files, strnd)]
    		bw_hf_condtn_lab <- names(bw_hf_files) %>% set_names(names(bw_hf_files))
    		bw_hf_condtn_cols <- c(met.brewer('Austria')[6], met.brewer('Austria')[3], met.brewer('Austria')[5], met.brewer('Austria')[1]) %>% set_names(bw_hf_condtn_lab)

        } else {
    
		bw_hf_files <- bigwig_files[str_detect(bigwig_files, strnd)]
    		bw_hf_condtn_lab <- names(bw_hf_files) %>% set_names(names(bw_hf_files))
    		bw_hf_condtn_cols <- c(met.brewer('Austria')[1], met.brewer('Austria')[3], met.brewer('Austria')[5]) %>% set_names(bw_hf_condtn_lab)

        }
    

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
                                   cex.axis = 1.3,
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


