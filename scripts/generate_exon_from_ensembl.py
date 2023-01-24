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

import pysam

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


def process_fasta(gtf, genome):

    # check if we have a fasta index
    # create index if not
    if not os.path.isfile(genome+""):
        pysam.index(genome)

    o = open(gtf + ".fasta", "a")

    pysam_obj = pysam.FastaFile(genome)

    with open(gtf, 'r') as gtf_file:
        for line in gtf_file:
            line = line.strip().split("\t")
            sam_pos = line[0]+":"+str(line[1])+"-"+str(line[2])

            # print()

            seq = pysam_obj.fetch(region=sam_pos)

            # print(seq)

            o.write(">" + line[3] + "\n")
            o.write(seq + "\n")

    o.close()

    return


def process_gzipped_gtf(uri):

    prev_gene = ""
    exon_num = 1

    # print("Downloading "+uri)
    # file_name, headers = urllib.request.urlretrieve(uri)
    # print("Done")

    o = open(uri + ".bed", "a")

    with gzip.open(uri, 'rb') as file:
        for line in file:

            # split by whitespace, not only tab to directly get
            # access to the description fields
            index = re.split('\s', line.decode('utf-8'))

            # only exon lines
            # make sure we're skipping the header
            if not index[0].startswith("#") and index[2] == "exon":

                # get desired fields
                chr, start, stop, strand, gene_id, gene = \
                    itemgetter(0, 3, 4, 6, 9, 19)(index)

                # remove quotation marks and ;
                gene_id = gene_id.replace("\"", "")
                gene_id = gene_id.replace(";", "")
                gene = gene.replace("\"", "")
                gene = gene.replace(";", "")

                idx = "_".join([chr, start, stop])

                # still name gene, increase exon number
                # but only if we did not see this exon yet!
                if idx not in line_dict and prev_gene == gene:
                    # print("okay")
                    exon_num = exon_num + 1
                    line_dict[idx] = 1
                    name = "__".join([gene, gene_id, str(exon_num)])

                    output = "\t".join([chr,
                                        start,
                                        stop,
                                        name,
                                        strand])
                    o.write(output+"\n")


                # new gene, reset exon numbers and gene name
                elif idx not in line_dict and prev_gene != gene:
                    prev_gene = gene
                    exon_num = 1
                    line_dict[idx] = 1
                    name = "__".join([gene, gene_id, str(exon_num)])

                    output = "\t".join([chr,
                                        start,
                                        stop,
                                        name,
                                        strand])

                    o.write(output+"\n")

    o.close()

                # print(output)

    return uri+".bed"

#process_gzipped_gtf(gtf_dict['rattus_norvegicus'])

# process_fasta("/mnt/big_data/genomes/", "/mnt/big_data/genomes/Rnor_6_0_96/Rnor_6_0_96.fa")
process_gzipped_gtf("/mnt/big_data/genomes/GRCh38_90/GRCh38.90.gtf.gz")
process_fasta("/mnt/big_data/genomes/GRCh38_90/GRCh38.90.gtf.gz.bed", "/mnt/big_data/genomes/GRCh38_90/GRCh38_90.fa")






