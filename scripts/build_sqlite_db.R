#!/usr/bin/env Rscript

# Copyright (C) 2023 Tobias Jakobi
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

options(echo=TRUE)

process_input_tables = function(data_file, column_names, species, bed_file, ribocirc_species){
  
  message(paste0("Processing species: ", species))
  
  message("Loading circAtlas table")
  # read orginal circatlas data, remove first line with broken headers
  input <- read.csv(data_file, sep="\t", header=F, skip=1)

  
  # fix headers
  colnames(input) <- column_names
  
 
  
  message("Loading circAtlas BED file")
  
  # read BED file for strand information
  bed <- read.csv(bed_file, sep="\t", header=T)
  
  # extract columns with name (for line matching) and strand
  strand_info <- bed[,c(4,5)]
  
  # fix column names
  colnames(strand_info) <- c("Strand","CircAtlas")

  message("Merging circAtlas main tables with BED file")
  # merge BED data with other data to recover strand information where possible
  # keep input column order -> https://stackoverflow.com/questions/17574882/how-to-prevent-merge-from-reordering-columns 
  input_with_strand <- merge(input,strand_info, all = T)[, union(names(input), names(strand_info))]
  
  # replace "-" in data frame with actual NAs for smoother SQLite import
  # however, make sure we only replace single "-", nothing that is part of
  # circRNA IDs
  input_with_strand <- as.data.frame(lapply(input_with_strand, function(y) gsub("^-$", NA, y)))

  # assign species name to column
  input_with_strand$Species <- species

  
  message("Loading circ2disease table")
  
  # include circ2disease data
  circ2disase <- circ2disase <- read.csv("../circ2disease/circ2disease_association_curated.csv", sep=",", header=T)[,c("circBase","Pubmed")]
  
  
  # replace N/A string with actual NA value
  circ2disase <- as.data.frame(lapply(circ2disase, function(y) gsub("^N/A$", NA, y)))
  
  message("Merging circ2disease data")
  
  # merge in circ2disease
  # this results in a few circRNAs to be duplicated,
  # as they have more then one Pubmed result

  input_with_strand <- merge(circ2disase, input_with_strand, by.x=c("circBase"),  by.y=c("circBase"), all.y=T)[, union(names(circ2disase), names(input_with_strand))]
  

  message("Loading ribocirc table")
  
  ribocirc <- read.csv("../ribocirc/All_condition_independent_riboCIRC_meta.csv", sep="\t", header=T)
  
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


  message("Splitting by genome build")
  
  # two new data frames, each only with one set of coordinates
  build1 <- input_with_strand[,-c(5)]
  build2 <- input_with_strand[,-c(6)]
  
  names(build1)[c(5)] <- 'build1'
  names(build2)[c(5)] <- 'build2'
  build1$Genome <- column_names[3]
  build2$Genome <- column_names[4]
  
  message("Loading Arraystar data")
  
  arraystar <- read.csv(paste0("../arraystar/",ribocirc_species,"_mapping_array_final.csv"), sep="\t", header=T)
  colnames(arraystar) <- c("build2","build1","Arraystar")

  message("Merging genome build 1")
  build1 <- merge(arraystar[,c(2,3)],build1, by=c("build1"), all.y = T)
  message("Merging genome build 2")
  build2 <- merge(arraystar[,c(1,3)],build2, by=c("build2"), all.y = T)
  
  # integrate exorbase, only human data available
  if (species == "homo_sapiens"){
    
    message("Loading exorbase data")
    
    exorbase <- read.csv("../exorbase/circRNAs_anno.txt", sep="\t", header=T)[,c("circID","Genomic.position")]
    colnames(exorbase) <- c("exorBase2","Genomic.position")

    # adjust position to match the circatlas format
    exorbase$Genomic.position<-gsub("-","|",exorbase$Genomic.position)

    # merge based on hg38 coordinates

    message("Merging genome build 1")
    build1 <- merge(exorbase,build1, by.x = c("Genomic.position"), by.y = c("build1"), all.y = T)
    
    message("Merging genome build 2")
    build2 <- merge(exorbase,build2, by.x = c("Genomic.position"), by.y = c("build2"), all.y = T)
    
  } else {
    build1$exorBase2 <- NA
    build2$exorBase2 <- NA
  }

  
  message("Creating multi-field genomic coordinates")
  
  # break up coordinate column into columns for chr/start/stop
  
  column <- 1
  
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
           "CircAtlas",
           "hg19",
           "hg38",
           "circBase",
           "circRNADB",
           "Deepbase2",
           "Circpedia2"
)

mouse <- c("Species",
           "CircAtlas",
           "mm9",
           "mm10",
           "circBase",
           "circRNADB",
           "Deepbase2",
           "Circpedia2"
)

rat <- c("Species",
         "CircAtlas",
         "rn5",
         "rn6",
         "circBase",
         "circRNADB",
         "Deepbase2",
         "Circpedia2"
         )


human_data <- process_input_tables("../circatlas/hg38_hg19_v2.0.txt", human, "homo_sapiens", "../circatlas/human_bed_v2.0.txt","Human")
mouse_data <- process_input_tables("../circatlas/mm10_mm9_v2.0.txt", mouse, "mus_musculus", "../circatlas/mouse_bed_v2.0.txt","Mouse")
rat_data <- process_input_tables("../circatlas/rn6_rn5_v2.0.txt", rat, "rattus_norvegicus", "../circatlas/rat_bed_v2.0.txt","Rat")

# binning together species data
final_data <- list(human_data,mouse_data,rat_data)

# generate main dataframe
bind_dataframe <- do.call(rbind, final_data )

# reorder table fields
bind_dataframe <- bind_dataframe[,c(
  "Species",
  "circBase",
  "CircAtlas",
  "circRNADB",
  "Deepbase2",
  "Circpedia2",
  "riboCIRC",
  "exorBase2",
  "Arraystar",
  "Chr",
  "Start",
  "Stop",
  "Strand",
  "Genome",
  "Pubmed"
)]

# write to CSV file for easy SQLite3 export
# no column names since we set those in the SQL schema
write.table(file="../sqlite_export.csv",bind_dataframe, sep=",", row.names=F, quote = F, col.names = F)

