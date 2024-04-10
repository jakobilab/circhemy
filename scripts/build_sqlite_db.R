#!/usr/bin/env Rscript

# Copyright (C) 2024 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# let's now get the command line arguments

options(echo=FALSE)
suppressMessages(suppressWarnings(library(R.utils)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(biomaRt)))

process_input_tables = function(data_file, column_names, species, bed_file, ribocirc_species){
  
  message(paste0("Processing species: ", species))
  
  message("Loading CircAtlas2 table")
  # read orginal CircAtlas2 data, remove first line with broken headers
  input <- read.csv(data_file, sep="\t", header=F, skip=1)
  
  # fix headers
  colnames(input) <- column_names

  message("Loading CircAtlas2 BED file")
  
  # read BED file for strand information
  bed <- read.csv(bed_file, sep="\t", header=T)
  
  # extract columns with name (for line matching) and strand
  strand_info <- bed[,c(4,5)]
  
  # fix column names
  colnames(strand_info) <- c("Strand","CircAtlas2")

  message("Merging CircAtlas2 main tables with BED file")
  # merge BED data with other data to recover strand information where possible
  # keep input column order -> https://stackoverflow.com/questions/17574882/how-to-prevent-merge-from-reordering-columns 
  input_with_strand <- merge(input,strand_info, all = T)[, union(names(input), names(strand_info))]
  
  # replace "-" in data frame with actual NAs for smoother SQLite import
  # however, make sure we only replace single "-", nothing that is part of
  # circRNA IDs
  input_with_strand <- as.data.frame(lapply(input_with_strand, function(y) gsub("^-$", NA, y)))

  # assign species name to column
  input_with_strand$Species <- species

  # only available for human and mouse (hg19/mm9)
  if (species == "homo_sapiens" || species == "mus_musculus"){
    message("Splicing in alternate circBase IDs")
    circbase <- read.csv(paste0("../circhemy/data/circbase/",species,"_circBase_alt.csv"), sep="\t", header=T)[,c("circBase","circBase_alt")]
    input_with_strand <- merge(circbase, input_with_strand, by.x=c("circBase"),  by.y=c("circBase"), all.y=T)[, union(names(circbase), names(input_with_strand))]
  } else {
    circbase <- as.data.frame(rep(NA,nrow(input_with_strand)))
    colnames(circbase) <- "circBase_alt"
    input_with_strand <-cbind(circbase,input_with_strand)
  }

  # fix column order
  input_with_strand <- input_with_strand[,c("Species", "CircAtlas2", column_names[3], column_names[4], "circBase", "circBase_alt", "circRNADb", "deepBase2", "Circpedia2", "Strand")]

  message("Pre-reading CSNv1 data to get gene names")

  csn <- read.csv(paste0("../circhemy/data/csnv1/",ribocirc_species,"_CSNv1.csv"), sep="\t", header=F)
  colnames(csn) <- c("CircAtlas2","CSNv1","Error","build2","Source")
  csn <- csn[,c("CSNv1","build2")]

  # We will get the gene names from the BLAST results
  csn$Gene <- csn$CSNv1
  csn$Gene <- gsub("circ", "", csn$Gene)
  csn$Gene <- gsub("\\(.*", "", csn$Gene)

  csn$CSNv1 <- NULL
  colnames(csn) <- c(colnames(input_with_strand)[3], "Gene")

  input_with_strand <- merge(csn,input_with_strand, by=colnames(input_with_strand)[3], all.y = T)[, union(names(input_with_strand), "Gene")]

  annotation <- paste0(substr(str_split(species,"_")[[1]][1],1,1),str_split(species,"_")[[1]][2],"_gene_ensembl")
  ensembl <- useMart("ensembl", dataset = annotation, host = "https://useast.ensembl.org")

  bm <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "wikigene_description"), values=input_with_strand$Gene,  mart = ensembl)

  idx_bm <- match(input_with_strand$Gene, bm$external_gene_name)
  input_with_strand$Description <- bm$wikigene_description[idx_bm]
  input_with_strand$ENSEMBL <- bm$ensembl_gene_id[idx_bm]
  input_with_strand$Entrez <- bm$entrezgene_id[idx_bm]

  message("Loading circ2disease table")
  
  # include circ2disease data
  circ2disease <- read.csv("../circhemy/data/circ2disease/circ2disease_association_curated.csv", sep=",", header=T)[,c("circBase","Pubmed")]
  
  # replace N/A string with actual NA value
  circ2disease <- as.data.frame(lapply(circ2disease, function(y) gsub("^N/A$", NA, y)))
  
  message("Merging circ2disease data")
  
  # merge in circ2disease
  # this results in a few circRNAs to be duplicated,
  # as they have more then one Pubmed result

  input_with_strand <- merge(circ2disease, input_with_strand, by.x=c("circBase"),  by.y=c("circBase"), all.y=T)[, union(names(circ2disease), names(input_with_strand))]

  message("Loading ribocirc table")
  
  ribocirc <- read.csv("../circhemy/data/ribocirc/All_condition_independent_riboCIRC_meta.txt", sep="\t", header=T)
  
  ribocirc$Genome_pos<-gsub("-","|",ribocirc$Genome_pos)
  
  # chromosome
  ribocirc$Chr <- gsub(':.*','',ribocirc[,5])
  
  # start & stop for genome build 1
  ribocirc$Start <- as.numeric(gsub('\\|.*','',gsub('.*:','',ribocirc[,5])))
  ribocirc$Stop <- gsub('.*\\|','',gsub('.*:','',ribocirc[,5]))
  
  ribocirc$Start <- ribocirc$Start+1
  
  message("Fixing ribocirc genome positions")
  
  ribocirc$Genome_pos <- paste0(ribocirc$Chr,":",ribocirc$Start,"|",ribocirc$Stop)
  
  message("Selecting ribocirc species")
  ribocirc <- ribocirc[ribocirc$Species==ribocirc_species,]
  
  ribocirc<- ribocirc[,c(3,5)]
  colnames(ribocirc) <- c("riboCIRC",colnames(input_with_strand[5]))
  
  message("Merging riboCIRC data")
  
  input_with_strand <- merge(ribocirc,input_with_strand, all.y = T)[, union(names(input_with_strand), names(ribocirc))]

  # integrate circbank, only human data available
  if (species == "homo_sapiens"){

    # we only need first two columns: circbank ID and circbase
    # as we can match 100% of the rows via circbase
    circbank <- read.csv(paste0("../circhemy/data/circbank/circBank_circrna_annotation.txt"), sep="\t", header=T)[,c(1,2)]
    colnames(circbank) <- c("circBank","circBase")

    input_with_strand <- merge(circbank,input_with_strand, by = c("circBase"), all.y = T)[, union(names(circbank), names(input_with_strand))]

  } else {
    circbank <- as.data.frame(rep(NA,nrow(input_with_strand)))
    colnames(circbank) <- "circBank"
    input_with_strand <-cbind(circbank,input_with_strand)
    # input_with_strand$circbank <- NA
  }

  message("Splitting by genome build")

  # two new data frames, each only with one set of coordinates
  build1 <- input_with_strand[,-c(6)]
  build2 <- input_with_strand[,-c(7)]

  names(build1)[c(6)] <- 'build1'
  names(build2)[c(6)] <- 'build2'
  build1$Genome <- column_names[4]
  build2$Genome <- column_names[3]

  message("Loading Arraystar data")
  
  arraystar <- read.csv(paste0("../circhemy/data/arraystar/",ribocirc_species,"_mapping_array_final.csv"), sep="\t", header=T)
  colnames(arraystar) <- c("build2","build1","Arraystar")

  message("Merging genome build 1")
  build1 <- merge(arraystar[,c(2,3)],build1, by=c("build1"), all.y = T)
  message("Merging genome build 2")
  build2 <- merge(arraystar[,c(1,3)],build2, by=c("build2"), all.y = T)

  message("Loading CircRNA Standard Nomenclature (CSN) v1 data")

  csn <- read.csv(paste0("../circhemy/data/csnv1/",ribocirc_species,"_CSNv1.csv"), sep="\t", header=F)
  colnames(csn) <- c("CircAtlas2","CSNv1","Error","build2","Source")
  csn <- csn[,c("CSNv1","build2")]

  message("Merging genome build 2 only")
  build2 <- merge(build2,csn, by=c("build2"), all.x = T)

  # add empty column for old genome build
  csn <- as.data.frame(rep(NA,nrow(build1)))
  colnames(csn) <- "CSNv1"
  build1 <-cbind(build1,csn)

  # integrate exorbase, only human data available
  if (species == "homo_sapiens"){
    
    message("Loading exorbase data")
    
    exorbase <- read.csv("../circhemy/data/exorbase/circRNAs_anno.txt", sep="\t", header=T)[,c("circID","Genomic.position")]

    # adjust position to match the CircAtlas2 format
    exorbase$Genomic.position<-gsub("-","|",exorbase$Genomic.position)

    # merge based on hg38 coordinates

    message("Merging genome build 1")
    colnames(exorbase) <- c("exoRBase2","build1")
    build1 <- merge(exorbase,build1, by = c("build1"), all.y = T)[, union(names(exorbase), names(build1))]

    message("Merging genome build 2")
    colnames(exorbase) <- c("exoRBase2","build2")
    build2 <- merge(exorbase,build2, by = c("build2"), all.y = T)[, union(names(exorbase), names(build2))]
    
  } else {
    exoRBase2 <- as.data.frame(rep(NA,nrow(build1)))
    colnames(exoRBase2) <- "exoRBase2"
    build1 <-cbind(exoRBase2,build1)

    exoRBase2 <- as.data.frame(rep(NA,nrow(build2)))
    colnames(exoRBase2) <- "exoRBase2"
    build2 <-cbind(exoRBase2,build2)
  }

  message("Creating multi-field genomic coordinates")
  
  # break up coordinate column into columns for chr/start/stop
  
  column <- 2
  
  message("Removing non-existent genomic coordinates for genome builds")

  # remove columns w/o coordinates for the corresponding genome build
  build1 <- build1[!is.na(build1[,c(column)]),]
  build2 <- build2[!is.na(build2[,c(column)]),]

  # chromosome
  build1$Chr <- gsub(':.*','',build1[,column])
  build2$Chr <- gsub(':.*','',build2[,column])

  # start & stop for genome build 1
  build1$Start <- gsub('\\|.*','',gsub('.*:','',build1[,column]))
  build1$Stop <- gsub('.*\\|','',gsub('.*:','',build1[,column]))
  
  # start & stop for genome build 2
  build2$Start <- gsub('\\|.*','',gsub('.*:','',build2[,column]))
  build2$Stop <- gsub('.*\\|','',gsub('.*:','',build2[,column]))
  
  build1[column] <- NULL
  build2[column] <- NULL

  return_data <- rbind(build1,build2)
  
  num_entries <- format(round(nrow(return_data), 1), nsmall=0, big.mark=",")
  
  message(paste0("Species tables for ",species," has ",num_entries, " entries."))
  
  return(return_data)
}

human <- c("Species",
           "CircAtlas2",
           "hg38",
           "hg19",
           "circBase",
           "circRNADb",
           "deepBase2",
           "Circpedia2"
)

mouse <- c("Species",
           "CircAtlas2",
           "mm10",
           "mm9",
           "circBase",
           "circRNADb",
           "deepBase2",
           "Circpedia2"
)

rat <- c("Species",
         "CircAtlas2",
         "rn6",
         "rn5",
         "circBase",
         "circRNADb",
         "deepBase2",
         "Circpedia2"
         )

human_data <- process_input_tables("../circhemy/data/circatlas/hg38_hg19_v2.0.txt", human, "homo_sapiens", "../circhemy/data/circatlas/human_bed_v2.0.txt","Human")
mouse_data <- process_input_tables("../circhemy/data/circatlas/mm10_mm9_v2.0.txt", mouse, "mus_musculus", "../circhemy/data/circatlas/mouse_bed_v2.0.txt","Mouse")
rat_data <- process_input_tables("../circhemy/data/circatlas/rn6_rn5_v2.0.txt", rat, "rattus_norvegicus", "../circhemy/data/circatlas/rat_bed_v2.0.txt","Rat")

message("Starting final data merge for CSV export")
# binning together species data
final_data <- list(human_data,mouse_data,rat_data)

# generate main dataframe
bind_dataframe <- do.call(rbind, final_data )

# reorder table fields
bind_dataframe <- bind_dataframe[,c(
  "Species",
  "Gene",
  "Description",
  "ENSEMBL",
  "Entrez",
  "circBase",
  "circBase_alt",
  "CircAtlas2",
  "circRNADb",
  "deepBase2",
  "Circpedia2",
  "circBank",
  "riboCIRC",
  "exoRBase2",
  "Arraystar",
  "CSNv1",
  "Chr",
  "Start",
  "Stop",
  "Strand",
  "Genome",
  "Pubmed"
)]

message("Writing CSV export")
# write to CSV file for easy SQLite3 export
# no column names since we set those in the SQL schema
write.table(bind_dataframe, file= "../circhemy/data/circhemy_data.csv", row.names=F, quote = F, col.names = F, sep="\t")

message("Bzip2 CSV export")
# create a nice, small package, it's all text data after all
bzip2("../circhemy/data/circhemy_data.csv", overwrite=T)
