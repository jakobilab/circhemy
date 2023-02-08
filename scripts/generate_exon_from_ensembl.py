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

import re
import gzip
import os
import urllib.request
from operator import itemgetter
import io

import pysam

import pybedtools
from pybedtools import BedTool

gtf_dict = {
    "homo_sapiens": "https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz",
    "mus_musculus": "https://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz",
    # "rattus_norvegicus": "https://ftp.ensembl.org/pub/release-90/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf.gz"
    "rattus_norvegicus": "/mnt/big_data/genomes/Rnor_6_0_96/Rnor_6.0.96.gtf.gz"
    }

genome_dict = {
    "homo_sapiens": "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz",
    "mus_musculus": "https://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz",
    # "rattus_norvegicus": "https://ftp.ensembl.org/pub/release-90/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf.gz"
    "rattus_norvegicus": "/mnt/big_data/genomes/Rnor_6_0_96/Rnor_6_0_96.fa"
    }

line_dict = {}

bad_gene = {"DTX2P1-UPK3BP1-PMS2P11"}


def process_fasta(bed_file, genome):

    # check if we have a fasta index
    # create index if not
    if not os.path.isfile(genome+""):
        pysam.index(genome)

    os.system("rm " + bed_file + ".fasta")

    o = open(bed_file + ".fasta", "a")

    pysam_obj = pysam.FastaFile(genome)

    with open(bed_file, 'r') as bed:
        for line in bed:
            line = line.strip().split("\t")
            sam_pos = line[0]+":"+str(line[1])+"-"+str(line[2])

            # no NOT write out gene bodies into FASTA
            if line[3].split("!")[5] == "g":
                continue

            seq = pysam_obj.fetch(region=sam_pos)

            o.write(">" + line[3] + "\n")
            o.write(seq + "\n")

    o.close()

    return


def process_gzipped_gtf(uri):

    prev_gene = ""
    exon_num = 1
    five_num = 1
    three_num = 1
    line_num = 1

    # TODO: set up download

    # print("Downloading "+uri)
    # file_name, headers = urllib.request.urlretrieve(uri)
    # print("Done")

    os.system("rm "+uri+".bed")
    os.system("rm "+uri + "_exon.bed")
    os.system("rm "+uri + "_gene.bed")
    os.system("rm "+uri + "_final.bed")

    exon_output_file = open(uri + "_exon.bed", "a")
    gene_output_file = open(uri + "_gene.bed", "a")

    # TODO: remove temp files

    print("Cleaning up")

    # remove old file, we have to append to we have to make sure we start fresh

    print("Sort GTF by coordinates to ensure correct exon numbering")

    # poor man's sorting by chromosome and start coordinate
    # also remove comment lines from gtf
    # will only work with unix/linux, not windows

    # os.system("zcat " +
    #           uri +
    #           " | grep -v \"#!\" " +
    #           "| awk '{if ($7==\"+\"){print}}'" +
    #           "|  sort -k1,1 -k4,4n | gzip -c > " + uri +
    #           ".plus")
    #
    # os.system("zcat " +
    #           uri +
    #           " | grep -v \"#!\" " +
    #           "| awk '{if ($7==\"-\"){print}}'" +
    #           "|  sort -k1,1 -k4,4rn | gzip -c > " + uri +
    #           ".minus")

    # os.system("rm "+uri)

    todo_list = [uri + ".plus", uri + ".minus"]

    for gtf_file in todo_list:

        print("Processing GTF "+gtf_file)

        with gzip.open(gtf_file, 'rb') as file:
            for line in file:
                print("processed " + str(line_num) + " lines",  end="\r")
                line_num = line_num + 1

                # split by whitespace, not only tab to directly get
                # access to the description fields
                index = re.split('\s', line.decode('utf-8'))

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
                        # print("okay")
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

                        # print("okay")
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

                elif index[2] == "gene":

                    # get desired fields
                    chr, start, stop, strand, gene_id, gene, type = \
                        itemgetter(0, 3, 4, 6, 9, 13, 17)(index)

                    # remove quotation marks and ;
                    gene_id = gene_id.replace("\"", "")
                    gene_id = gene_id.replace(";", "")
                    gene_id = gene_id[-6:]

                    gene = gene.replace("\"", "")
                    gene = gene.replace(";", "")
                    gene = gene[0:12]
                    #
                    # type = type.replace("\"", "")
                    # type = type.replace(";", "")
                    # type = type.replace("protein_coding", "coding")
                    # type = type[0:6]


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

    # close files, done with first step
    exon_output_file.close()
    gene_output_file.close()

    # now, lets build an intron list by removing exons from the genes
    # using bedtools subtract

    exons = BedTool(uri + "_exon.bed")
    genes = BedTool(uri + "_gene.bed")

    print("Extracting introns")

    introns = str(genes.subtract(exons, s=True))

    final_output_file = open(uri + "_final.bed", "a")

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

    #os.system("rm "+uri+".fixed")

    return uri + "_final.bed"

#process_gzipped_gtf(gtf_dict['rattus_norvegicus'])

# process_fasta("/mnt/big_data/genomes/", "/mnt/big_data/genomes/Rnor_6_0_96/Rnor_6_0_96.fa")
process_gzipped_gtf("/home/tjakobi/repos/jakobilab/circhemy/data/blast/GRCh38.90.gtf.gz")
process_fasta("/home/tjakobi/repos/jakobilab/circhemy/data/blast/GRCh38.90.gtf.gz_final.bed", "/home/tjakobi/repos/jakobilab/circhemy/data/blast/GRCh38_90.fa")






