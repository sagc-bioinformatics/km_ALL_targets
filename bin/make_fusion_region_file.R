
####################################
#                                  #
#   Create custom fusion targets   #
#                                  #
####################################

## variation on generate_fusion_target_sequence_regions using rtracklayer to read gtf file

# load required libraries
suppressMessages(library(tidyverse, warn.conflicts = F, quietly = T))
suppressMessages(library(rtracklayer, warn.conflicts = F, quietly = T))

# collect arguments
args <- commandArgs(trailingOnly = TRUE)

# args[1] text file containing gene and exon information for targets to be generated
# args[2] location of gtf file

# load target information file
target_file <- read_delim(file = args[1], delim = "\t", col_names = FALSE, comment = "#", 
                          col_types = cols(
                            X1 = col_character(),
                            X2 = col_character(),
                            X3 = col_character(),
                            X4 = col_character(),
                            X5 = col_double()
                          ))

# extract information listing exon ranges to include in targets
exon_ranges <- suppressWarnings(tidyr::separate(target_file, X2, c("gene1_exon_start","gene1_exon_end"), sep = ":") %>% 
  separate(X4, c("gene2_exon_start","gene2_exon_end"), sep = ":"))
# replace_na should some genes have only single exon listed
exon_ranges <- exon_ranges %>% 
  mutate(gene1_exon_end = ifelse(!is.na(gene1_exon_end), gene1_exon_end, paste(gene1_exon_start))) %>% 
  mutate(gene2_exon_end = ifelse(!is.na(gene2_exon_end), gene2_exon_end, paste(gene2_exon_start))) %>% 
  mutate_at(c('gene1_exon_start','gene1_exon_end','gene2_exon_start','gene2_exon_end'), as.numeric)
           


# load gtf file
gtf_file <- args[2]
gtf_obj <- import(gtf_file)


#############################
#                           #
#    Gene Fusion targets    #
#                           #
#############################

## Run the following function for each line in the target_file

for(i in 1:nrow(target_file)){

    ## extract transcript information for given gene pair, limiting to exonic regions
    # gene1
    gene1_exon_data <- data.frame(gtf_obj[gtf_obj$gene_name == target_file$X1[i] & gtf_obj$type == "exon", ])
    # gene2
    gene2_exon_data <- data.frame(gtf_obj[gtf_obj$gene_name == target_file$X3[i] & gtf_obj$type == "exon", ])
    
    ## limit regions to protein_coding transcripts
    ## exclude any duplicate exonic regions from separate transcripts
    
    ### Need an Error Check to confirm there is gene data & print message if not ###
    
    # gene1
    gene1_exon_data <- gene1_exon_data %>% 
      # filter(source == "protein_coding") %>% 
      distinct(start, end, .keep_all = TRUE) %>% 
      mutate(exon_number = as.numeric(exon_number))
      
    # gene2
    gene2_exon_data <- gene2_exon_data %>% 
      # filter(source == "protein_coding") %>% 
      distinct(start, end, .keep_all = TRUE) %>% 
      mutate(exon_number = as.numeric(exon_number))
      
    
    ########## calculating target regions ##########
    
    ## gene1 (+) strand OR gene2 (-) strand
    # mutate(gene1_region = paste(chromosome_name,":",(exon_chrom_end - length_gene_region),"-", exon_chrom_end, sep = ""))
    ## gene2 (+) strand OR gene1 (-) strand
    # mutate(gene2_region = paste(chromosome_name,":",exon_chrom_start,"-",(exon_chrom_start + length_gene_region), sep = ""))
    
    ################################################
    
    
    ## Gene 1
    gene1_exon_data <- gene1_exon_data %>% 
      # generate field containing target region
      mutate(gene1_region = ifelse(strand == "+", 
                                   paste(seqnames, ":", (end - (target_file$X5[i]/2 -1)), "-", (end), sep = ""), 
                                   paste(seqnames, ":", start, "-", (start + (target_file$X5[i]/2 -1)), sep = ""))) %>% 
      # limit target regions to requested exons
      filter(exon_number >= exon_ranges$gene1_exon_start[i] & exon_number <= exon_ranges$gene1_exon_end[i]) %>% 
      # generate field with exon_number
      mutate(exon_number = paste("exon", exon_number, sep = "")) %>% 
      # select required fields
      dplyr::select(gene1_region, gene_name, strand, transcript_name, exon_number) %>% 
      dplyr::rename(gene1_name = gene_name, gene1_strand = strand, gene1_transcript = transcript_name, gene1_exon = exon_number) %>% 
      # eliminate any duplicated regions
      dplyr::distinct(gene1_region, .keep_all = TRUE)
    
    ## Gene 2
    gene2_exon_data <- gene2_exon_data %>% 
      # generate field containing target region
      mutate(gene2_region = ifelse(strand == "+", 
                                   paste(seqnames, ":", start, "-", (start + (target_file$X5[i]/2 -1)), sep = ""), 
                                   paste(seqnames, ":", (end - (target_file$X5[i]/2 -1)) ,"-", end, sep = ""))) %>% 
      # limit target regions to requested exons
      filter(exon_number >= exon_ranges$gene2_exon_start[i] & exon_number <= exon_ranges$gene2_exon_end[i]) %>% 
      # generate field with exon_number
      mutate(exon_number = paste("exon", exon_number, sep = "")) %>% 
      # select required fields
      dplyr::select(gene2_region, gene_name, strand, transcript_name, exon_number) %>% 
      dplyr::rename(gene2_name = gene_name, gene2_strand = strand, gene2_transcript = transcript_name, gene2_exon = exon_number) %>% 
      # eliminate any duplicated regions
      dplyr::distinct(gene2_region, .keep_all = TRUE)
    
    
    # specify output fileName
    file = paste("FusionRegions", "_", target_file$X1[i], "-", target_file$X2[i], "_", target_file$X3[i], "-", target_file$X4[i], ".txt", sep = "")
    # outLocation = paste("~/km_ALL_targets/custom_targets/", sep = "")        
    # outFile  = paste(outLocation, file, sep = "")
    
    
    for(x in 1:nrow(gene1_exon_data)){
      
      for(y in 1:nrow(gene2_exon_data)){
        
        bind_cols(gene1_exon_data[x,], gene2_exon_data[y,]) %>% 
          write_tsv(file = file, 
                    append = TRUE)
        
      }
    }
    
}    
    
