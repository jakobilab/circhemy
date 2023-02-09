#!/usr/bin/env python3
# Copyright (C) 2023 Tobias Jakobi
#
# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu
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

import argparse
import re
import gzip
import os
import urllib.request
from operator import itemgetter
import io

import pysam
from pybedtools import BedTool

gtf_dict = {
    "homo_sapiens": "https://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz",
    "mus_musculus": "https://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm39.90.gtf.gz",
    "rattus_norvegicus": "https://ftp.ensembl.org/pub/release-90/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf.gz"
    }

genome_dict = {
    "homo_sapiens": "https://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    "mus_musculus": "https://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz",
    "rattus_norvegicus": "https://ftp.ensembl.org/pub/release-90/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"
    }

line_dict = {}

parser = argparse.ArgumentParser(
    prog="naming_conversion_data_preparation.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    fromfile_prefix_chars="@",
    description="Prepares external ENSEMBL data sources for conversion of "
                "CircAtlas circRNA IDs to the new naming scheme proposed in "
                "Chen et al. 2023.\n"
                "\n"
                "Version 0.0.2\n"
                "\n"
                "https://github.com/jakobilab/circhemy\n"
                "https://jakobilab.org\n"
                "tjakobi@arizona.edu",

    usage=""" naming_conversion_data_preparation [<args>]"""
)

group = parser.add_argument_group("output parameters")

group.add_argument("-o",
                   "--output",
                   dest="output_directory",
                   default="./",
                   help="The output folder for files created the script",
                   required=True
                   )

group.add_argument("-s",
                   "--species",
                   dest="species",
                   help="Species to process",
                   choices=("mus_musculus", "homo_sapiens", "rattus_norvegicus"),
                   default=["mus_musculus", "homo_sapiens", "rattus_norvegicus"]
)

args = parser.parse_args()


def download_data(file_url):

    # get base file name
    file_name = file_url.split("/")[-1]

    file_full_path = args.output_directory + "/" + file_name

    # check if already downloaded (check unzipped file)
    if os.path.isfile(file_full_path.replace(".gz", "")):
        print("Found required file " + file_full_path.replace(".gz", ""))
    else:
        print("Downloading " + file_url)
        urllib.request.urlretrieve(file_url, file_full_path)
        print("Unpacking")
        os.system("gzip -d " + file_full_path)
        print("Done")

    return file_full_path.replace(".gz", "")


def process_fasta(bed_file, genome):

    # check if we have a fasta index
    # create index if not
    if not os.path.isfile(genome+""):
        pysam.index(genome)

    os.system("rm -f " + bed_file + ".fasta")

    o = open(bed_file + ".fasta", "a")

    pysam_obj = pysam.FastaFile(genome)

    with open(bed_file, 'r') as bed:
        for line in bed:
            line = line.strip().split("\t")
            sam_pos = line[0]+":"+str(line[1])+"-"+str(line[2])

            # no NOT write out gene bodies into FASTA
            if line[3].split("!")[5] == "G":
                continue

            seq = pysam_obj.fetch(region=sam_pos)

            o.write(">" + line[3] + "\n")
            o.write(seq + "\n")

    o.close()

    return bed_file + ".fasta"


def process_gzipped_gtf(gzipped_gtf):

    prev_gene = ""
    exon_num = 1
    five_num = 1
    three_num = 1
    line_num = 1

    os.system("rm -f " + gzipped_gtf + ".bed")
    os.system("rm -f " + gzipped_gtf + "_exon.bed")
    os.system("rm -f " + gzipped_gtf + "_gene.bed")
    os.system("rm -f " + gzipped_gtf + "_final.bed")

    exon_output_file = open(gzipped_gtf + "_exon.bed", "a")
    gene_output_file = open(gzipped_gtf + "_gene.bed", "a")

    # TODO: remove temp files

    print("Cleaning up")

    # remove old file, we have to append to we have to make sure we start fresh

    print("Sort GTF by coordinates to ensure correct exon numbering")

    # poor man's sorting by chromosome and start coordinate
    # also remove comment lines from gtf
    # will only work with unix/linux, not windows

    os.system("cat " +
              gzipped_gtf +
              " | grep -v \"#!\" " +
              "| awk '{if ($7==\"+\"){print}}'" +
              "|  sort -k1,1 -k4,4n > " + gzipped_gtf +
              ".plus")

    os.system("cat " +
              gzipped_gtf +
              " | grep -v \"#!\" " +
              "| awk '{if ($7==\"-\"){print}}'" +
              "|  sort -k1,1 -k4,4rn > " + gzipped_gtf +
              ".minus")

    todo_list = [gzipped_gtf + ".plus", gzipped_gtf + ".minus"]

    for gtf_file in todo_list:

        print("Processing GTF "+gtf_file)

        with open(gtf_file, 'r') as file:
            for line in file:
                print("processed " + str(line_num) + " lines",  end="\r")
                line_num = line_num + 1

                # split by whitespace, not only tab to directly get
                # access to the description fields
                index = re.split('\s', line)

                # skip non-standard chromosomes for now
                # too long chromosome names will interfere with
                # BLAST DB creation as only 50 characters are allowed for
                # the name tag
                if len(index[0]) > 6:
                    continue

                # only exon lines
                # make sure we're skipping the header
                if index[2] == "exon":

                    # get desired fields
                    chr, start, stop, strand, gene_id, gene = \
                        itemgetter(0, 3, 4, 6, 9, 19)(index)

                    # remove quotation marks and ;
                    gene_id = gene_id.replace("\"", "")
                    gene_id = gene_id.replace(";", "")
                    gene_id = gene_id[-6:]
                    gene = gene.replace("\"", "")
                    gene = gene.replace(";", "")
                    gene = gene[0:12]

                    idx = "!".join([chr, start, stop])

                    # still name gene, increase exon number
                    # but only if we did not see this exon yet!
                    if idx not in line_dict and prev_gene == gene:
                        exon_num = exon_num + 1
                        line_dict[idx] = 1
                        name = "!".join([gene, gene_id, idx, "E", str(exon_num)])

                        output = "\t".join([chr,
                                            start,
                                            stop,
                                            name,
                                            str(0),
                                            strand])

                        exon_output_file.write(output+"\n")

                    # new gene, reset exon numbers and gene name
                    elif idx not in line_dict and prev_gene != gene:
                        prev_gene = gene
                        exon_num = 1
                        five_num = 1
                        three_num = 1
                        line_dict[idx] = 1
                        name = "!".join([gene, gene_id, idx, "E", str(exon_num)])

                        output = "\t".join([chr,
                                            start,
                                            stop,
                                            name,
                                            str(0),
                                            strand])

                        exon_output_file.write(output+"\n")

                elif index[2] == "three_prime_utr":

                    # get desired fields
                    chr, start, stop, strand, gene_id, gene = \
                        itemgetter(0, 3, 4, 6, 9, 17)(index)

                    # remove quotation marks and ;
                    gene_id = gene_id.replace("\"", "")
                    gene_id = gene_id.replace(";", "")
                    gene_id = gene_id[-6:]

                    gene = gene.replace("\"", "")
                    gene = gene.replace(";", "")
                    gene = gene[0:12]

                    idx = "!".join([chr, start, stop])

                    # still name gene, increase exon number
                    # but only if we did not see this exon yet!
                    if idx not in line_dict and prev_gene == gene:

                        name = "!".join([gene, gene_id, idx,  "U3", str(three_num)])

                        three_num = three_num + 1
                        line_dict[idx] = 1

                        output = "\t".join([chr,
                                            start,
                                            stop,
                                            name,
                                            str(0),
                                            strand])

                        exon_output_file.write(output+"\n")

                    # new gene, reset exon numbers and gene name
                    elif idx not in line_dict and prev_gene != gene:
                        prev_gene = gene
                        exon_num = 1
                        five_num = 1
                        three_num = 1
                        line_dict[idx] = 1
                        name = "!".join([gene, gene_id, idx, "U3", str(three_num)])

                        output = "\t".join([chr,
                                            start,
                                            stop,
                                            name,
                                            str(0),
                                            strand])

                        exon_output_file.write(output+"\n")

                elif index[2] == "five_prime_utr":

                    # get desired fields
                    chr, start, stop, strand, gene_id, gene = \
                        itemgetter(0, 3, 4, 6, 9, 17)(index)

                    # remove quotation marks and ;
                    gene_id = gene_id.replace("\"", "")
                    gene_id = gene_id.replace(";", "")
                    gene = gene.replace("\"", "")
                    gene = gene.replace(";", "")
                    gene = gene[0:12]
                    gene_id = gene_id[-6:]

                    idx = "!".join([chr, start, stop])

                    # still name gene, increase exon number
                    # but only if we did not see this exon yet!
                    if idx not in line_dict and prev_gene == gene:

                        name = "!".join([gene, gene_id, idx,  "U5", str(five_num)])

                        # print("okay")
                        five_num = five_num + 1
                        line_dict[idx] = 1

                        output = "\t".join([chr,
                                            start,
                                            stop,
                                            name,
                                            str(0),
                                            strand])

                        exon_output_file.write(output + "\n")

                    # new gene, reset exon numbers and gene name
                    elif idx not in line_dict and prev_gene != gene:
                        prev_gene = gene
                        exon_num = 1
                        five_num = 1
                        three_num = 1
                        line_dict[idx] = 1
                        name = "!".join([gene, gene_id, idx,  "U5", str(five_num)])

                        output = "\t".join([chr,
                                            start,
                                            stop,
                                            name,
                                            str(0),
                                            strand])

                        exon_output_file.write(output + "\n")

                elif index[2] == "gene" and len(index) >= 18:

                    # get desired fields
                    chr, start, stop, strand, gene_id, gene, type = \
                        itemgetter(0, 3, 4, 6, 9, 13, 17)(index)

                    # remove quotation marks and ;
                    gene_id = gene_id.replace("\"", "")
                    gene_id = gene_id.replace(";", "")
                    gene_id = gene_id[-6:]

                    gene = gene.replace("\"", "")
                    gene = gene.replace(";", "")
                    gene = gene[0:10]

                    idx = "!".join([chr, start, stop])

                    name = "!".join([gene, gene_id, idx, "G"])

                    output = "\t".join([chr,
                                        start,
                                        stop,
                                        name,
                                        str(0),
                                        strand])

                    gene_output_file.write(output + "\n")

        line_num = 1

    print("")
    print("Done")

    # close files, done with first step
    exon_output_file.close()
    gene_output_file.close()

    # now, lets build an intron list by removing exons from the genes
    # using bedtools subtract

    exons = BedTool(gzipped_gtf + "_exon.bed")
    genes = BedTool(gzipped_gtf + "_gene.bed")

    print("Extracting introns")

    introns = str(genes.subtract(exons, s=True))

    final_output_file = open(gzipped_gtf + "_final.bed", "a")

    final_output_file.write(str(exons))
    final_output_file.write(str(genes))

    for line in io.StringIO(introns):

        index = re.split('\s', line)

        chr, start, stop, name, score, strand = \
            itemgetter(0, 1, 2, 3, 4, 5)(index)

        gene, id = itemgetter(0, 1)(name.split("!"))

        newname = "!".join([gene, id, chr, start, stop, "RI"])

        output = "\t".join([chr,
                            start,
                            stop,
                            newname,
                            str(0),
                            strand])

        final_output_file.write(output + "\n")

    return gzipped_gtf + "_final.bed"


def initialize_blast_database(fasta_file):

    print("Initializing BLAST database")

    # create BLAST database
    # use hash_index , DB version 5 and parse-seqid to allow
    # targeted search using the -seqidlist parameter
    os.system("makeblastdb -in "
              + fasta_file +
              " -dbtype nucl -hash_index -blastdb_version 5 -parse_seqids")


def process_species(species):

    print("Working on species "+species)

    # download data
    gtf_file = download_data(gtf_dict[species])
    genome_file = download_data(genome_dict[species])

    # preprocess GTF file
    bed_file = process_gzipped_gtf(gtf_file)

    # process FASTA file, generate input for BLAST database
    fasta_file = process_fasta(bed_file, genome_file)

    # generate BLAST database
    initialize_blast_database(fasta_file)


# main program loop

if isinstance(args.species, list):
    for item in args.species:
        process_species(item)
else:
    process_species(args.species)
