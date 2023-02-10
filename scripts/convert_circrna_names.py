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
# The code for the multiprocessing text part is taken from
# https://stackoverflow.com/a/12293094

# Other code authored by Tobias Jakobi

import multiprocessing as mp
import sys
import tempfile

import io
import os
import re
import argparse

from pybedtools import BedTool
import subprocess


def run_blast_query(database,
                    fasta_query_file,
                    limit_bsl_file,
                    blast_output_file):

    # run customized BLAST query
    subprocess.run(["blastn", "-query", fasta_query_file,
                    "-db", database,
                    "-outfmt", "6 qseqid sseqid qstart qend "
                               "sstart send sstrand gaps mismatch pident qseq",
                    "-seqidlist", limit_bsl_file,
                    "-out", blast_output_file])
    return


def beautify_bed_names(bed_name_input, circrna_id):

    # check if intron or exon
    bed_name_input = bed_name_input.split("!")

    ensembl_prefix = ""

    # check which species we are in to get the correct ENSEMBL ID
    if circrna_id.split("-")[0] == "hsa":
        ensembl_prefix = "ENSG"
    elif circrna_id.split("-")[0] == "mmu":
        ensembl_prefix = "ENSMUSG"
    elif circrna_id.split("-")[0] == "rno":
        ensembl_prefix = "ENSRNOG"

    gene = bed_name_input[0]
    ensembl_id = bed_name_input[1]
    chr = bed_name_input[2]
    start = bed_name_input[3]
    end = bed_name_input[4]
    type = bed_name_input[5]

    output_name = gene + " | " + ensembl_prefix + "00000" + ensembl_id \
                  + " | " + type + " | " + chr + ":"+start+"-"+end

    return output_name



def get_circrna_boundaries_from_bedfile(bedfile):

    lines = bedfile.split("\n")

    # Coordinates from BED file are coming from BLAST
    # BLAST takes coordinates directly from the input GTF file
    # I.e. we are working with 1-BASED COORDINATES
    # circhemy uses circatlas as primary end point, e.g. coordinates
    # are also 1-based

    # first line
    # circ_start = int(lines[0].split("\t")[1])
    circ_start = int(lines[0].split("\t")[3].split("!")[3])

    # last line... very last item does not work, take -2 then
    # subtract 1 to get circRNA lined up with circatlas
    # again: circatlas: 1-based coordinates
    # circ_stop = int(lines[-2].split("\t")[2])-1
    circ_stop = int(lines[-2].split("\t")[3].split("!")[4])

    return circ_start, circ_stop


def process_blast_results(file_name, debug_bed_file, intersect_bed_file, bed_strand):

    # init empty string to hold bedtools string
    bedtools_string = ""

    exon_dict = {}
    intron_dict = {}

    distance_watcher_start = 0
    distance_watcher_stop = 0

    with open(file_name, 'r') as blast_result:

        item_dict = {}

        for line in blast_result:

            # get line in BLAST file
            line = line.strip().split("\t")

            # this is the SUBJECT name, e.g. the exon or introns
            # we parse location and other infos from there
            name = line[1].split("!")

            item_type = name[5]
            if item_type == "U5" or item_type == "U3" or item_type == "G":
                continue

            if line[1] in item_dict and item_dict[line[1]] >= float(line[9]):
                continue
            else:
                item_dict[line[1]] = float(line[9])

            # get start from the QUERY part of the BLAST results, i.e. the
            # circRNA name and location in the database
            circ_start, circ_stop = line[0].split("@")[1].split(":")[1].split("|")
            circ_chr = line[0].split("@")[1].split(":")[0].replace("chr", "")

            # get start and stop ane make INTs
            circ_start = int(circ_start)
            circ_stop = int(circ_stop)

            # we need the strand for orientation
            circ_strand = bed_strand

            # start/stop from QUERY (circRNA)
            item_chr = name[2]
            item_start = int(name[3])
            item_stop = int(name[4])

            # start/stop from SUBJECT (exons, intron)
            blast_start = int(line[4])
            blast_stop = int(line[5])

            if line[6] == "minus":
                # get the amount of the exon/intron we actually aligned
                # allows to see if we have a new internal exon
                # which is only part of an intron
                aligned_length = blast_start-blast_stop

                circ_item_start = item_start + blast_stop - 1
                circ_item_stop = item_start + blast_start

            else:
                aligned_length = blast_stop-blast_start

                circ_item_start = item_start + blast_start - 1
                circ_item_stop = item_start + blast_stop

            # will intron/exon dicts
            # used to later on detect too short exons that are continuing
            # to align into an intron

            # this is an exon
            if item_type == "E":
                tmp = {'start': circ_item_start,
                       'stop': circ_item_stop,
                       'id': line[1]
                       }
                # we use BED chr|start|stop as key because an intron may
                # match multiple times with different coordinates
                # i.e. an internal new exon AND an extension of an exon
                exon_dict[item_chr + "|" +
                          str(circ_item_start) + "|" +
                          str(circ_item_stop)] = tmp

            # this is an intron
            if item_type == "RI":
                tmp = {'start': circ_item_start,
                       'stop': circ_item_stop,
                       'id': line[1]
                       }

                intron_dict[item_chr + "|" +
                          str(circ_item_start) + "|" +
                          str(circ_item_stop)] = tmp

            # this is the exon/intron annotated overall length
            annotated_length = item_stop-item_start

            perc_aligned = aligned_length/annotated_length

            # if we have less than 50% coverage of the intron
            # AND we are in an intron we assume it's a new exon

            if item_type == "RI" and perc_aligned < 0.5:
                line[1] = line[1].replace("!RI",  "!NE")
                # print("aligned: " + str(perc_aligned))
                # print("not fully aligned intron, making it novel exon")

            name_tag = beautify_bed_names(bed_name_input=line[1],
                                          circrna_id=line[0].split("@")[0])

            bedtools_string = bedtools_string + "\t".join([item_chr,
                                                           str(circ_item_start),
                                                           str(circ_item_stop),
                                                           line[1]]
                                                          ) + "\n"

    # we are done with BLAST processing

    # now assign the exon numbers and magically discover too short and
    # too long exons

    for exon in exon_dict:

        exon_start = exon_dict[exon]['start']
        exon_stop = exon_dict[exon]['stop']
        exon_id = exon_dict[exon]['id']

        for intron in intron_dict:

            intron_start = intron_dict[intron]['start']
            intron_stop = intron_dict[intron]['stop']

            # now we can check potential overlaps

            # in this case the exon is extended by the intron on the right side
            if intron_start < exon_stop < intron_stop:
                # depending on the strand we have to put the L for long exon
                # left or right
                # if strand == plus: this would e.g. be 15L
                # if strand == minus: this would be L15, since at the end
                # for minus strand circRNAs the exons are brought into
                # numerical order and (L15, INTRON, 14) becomes
                # (14, INTRON, L15). This is why the L goes on the other side

                # search and replace exon ID in the bedtools string to
                # put the L on the correct side (eee above)

                if circ_strand == "+":

                    # build the replacement string (for plus strand)

                    replacement = exon_id+"L"
                    bedtools_string = bedtools_string.replace(exon_id,
                                                              replacement)

                else:
                    # build the replacement string (for minus strand)

                    replacement = exon_id.replace("!E!", "!E!L")
                    bedtools_string = bedtools_string.replace(exon_id,
                                                              replacement)

            # in this case the intron leaks into the exon on the left side
            # see above for details, here the switcheroo for the L is the
            # opposite compared to before

            elif intron_stop > exon_start > intron_start:

                if circ_strand == "+":
                    # build the replacement string (for plus strand)

                    replacement = exon_id.replace("!E!", "!E!L")
                    bedtools_string = bedtools_string.replace(exon_id,
                                                              replacement)

                else:
                    # build the replacement string (for minus strand)

                    replacement = exon_id+"L"
                    bedtools_string = bedtools_string.replace(exon_id,
                                                              replacement)

    # make new bedtools object
    bed_obj = BedTool(bedtools_string, from_string=True)

    # sort by coordinates
    bed_obj = bed_obj.sort()

    single_intron_check = False

    # merge contained exons, i.e. if exon 5 is 500 bp and exon 3 if 100 bp
    # completely within exon 5 only keep exon 5
    merged_blast_results = bed_obj.merge(c="4", o="first")

    debug_bed_file.write(str(merged_blast_results))

    return_list = []

    for line in io.StringIO(str(merged_blast_results)):
        line = line.strip().split("\t")

        tmp_list = line[3].split("!")[5:8]

        # remove E from exon annotation
        while "E" in tmp_list:
            tmp_list.remove("E")

        return_list.append("".join(tmp_list))

    # here we look at the actual start and stop position of the flanking
    # circRNA exons from the annotated exons from BLAST results

    canonical_exon_start, canonical_exon_stop = get_circrna_boundaries_from_bedfile(str(bed_obj))

    if len(return_list) == 1 and return_list[0] == "RI":
        single_intron_check = True

    # we now compare the exons boundaries with the circRNA annotation
    # this is strand independent

    # this just detects the outer border junctions of the circRNA
    # and compares with the references annotation
    # Exons that are overlapping with introns IN the circRNA are covered above

    if circ_start > canonical_exon_start and return_list[0] != "NE":
        # circRNA starts BEHIND annotated exon
        # i.e. first exon is shorter

        dist = circ_start - canonical_exon_start

        if dist > distance_watcher_start:
            distance_watcher_start = dist

        dist = str(dist)

        if circ_strand == "+":
            return_list[0] = "[" + dist + "]" + "L" + return_list[0]
        else:
            return_list[0] = return_list[0] + "L" + "[" + dist + "]"

    elif circ_start < canonical_exon_start and return_list[0] != "NE":
        # circRNA starts BEFORE annotated exon
        # i.e. first exon is longer

        dist = canonical_exon_start - circ_start

        if dist > distance_watcher_start:
            distance_watcher_start = dist

        dist = str(dist)

        if circ_strand == "+":
            return_list[0] = "[" + dist + "]" + "S" + return_list[0]
        else:
            return_list[0] = return_list[0] + "S" + "[" + dist + "]"

    if circ_stop > canonical_exon_stop and return_list[-1] != "NE":
        # circRNA stop BEHIND annotated exon
        # i.e. last exon is longer

        dist = circ_stop - canonical_exon_stop

        if dist > distance_watcher_stop:
            distance_watcher_stop = dist

        dist = str(dist)

        if circ_strand == "+":
            return_list[-1] = return_list[-1] + "S" + "[" + dist + "]"
        else:
            return_list[-1] = "[" + dist + "]" + "S" + return_list[-1]

    elif circ_stop < canonical_exon_stop and return_list[-1] != "NE":
        # circRNA starts BEFORE annotated exon
        # i.e. last exon is shorter

        dist = canonical_exon_stop - circ_stop

        if dist > distance_watcher_stop:
            distance_watcher_stop = dist

        dist = str(dist)

        if circ_strand == "+":
            return_list[-1] = return_list[-1] + "L" + "[" + dist + "]"
        else:
            return_list[-1] = "[" + dist + "]" + "L" + return_list[-1]

    # we should create a back-up for incomplete sequences
    # i.e. what if we can only match exon junction with BLAST?
    # in that case we should anyway insert a new exon NA into the name
    # if we have a somewhat fitting exon from the previous intersect

    # last check to get exons right

    # do we only have one exon and does that exon have a very high error
    # 5' or 3'? If so, we might want to fix that by adding a fitting exon
    # on 3'/5' as that is more probable than one exon that is > 20kb long.

    # we create another intersect with the exact start/stop circRNA position
    # and see which exons fits best

    start_fit = 0
    stop_fit = 0

    changed_start = False
    changed_stop = False

    if max(distance_watcher_start, distance_watcher_stop) > 100:

        bedtools_string = "\t".join([circ_chr,
                                     str(circ_start - 1000),
                                     str(circ_stop + 1000),
                                     "exon_fix"]) + "\n"

        circrna_junction = BedTool(bedtools_string, from_string=True)
        fix = intersect_bed_file.intersect(circrna_junction)

        start_exon_dist_dict = {}
        stop_exon_dist_dict = {}

        start_exon = 0
        stop_exon = 0

        for bed_line in io.StringIO(str(fix)):

            bed_line = bed_line.strip().split("\t")

            fix_start = int(bed_line[3].split("!")[3])
            fix_stop = int(bed_line[3].split("!")[4])
            fix_type = bed_line[3].split("!")[5]

            # we don't want introns here, just exons
            if fix_type == "U5" or fix_type == "U3" or fix_type == "RI" or fix_type=="G":
                continue

            exon_num = bed_line[3].split("!")[6]

            # look for start
            if distance_watcher_start > 100:
                # direct hit for start
                if fix_start == circ_start:
                    start_exon = exon_num
                else:
                    start_exon_dist_dict[circ_start-fix_start] = exon_num

            # look for stop
            if distance_watcher_stop > 100:
                # direct hit for stop
                if fix_stop == circ_stop:
                    stop_exon = exon_num
                else:
                    stop_exon_dist_dict[circ_stop-fix_stop] = exon_num

        tmp_exon = -1
        if distance_watcher_stop > 100 and len(stop_exon_dist_dict.keys()) > 0 and stop_exon == 0:
            fit = int(min(list(stop_exon_dist_dict.keys()), key=abs))

            tmp_exon = stop_exon_dist_dict[fit]

            if fit < 0 and bed_line[5] == "+":
                stop_exon = stop_exon_dist_dict[fit] + "L[" + str(abs(fit)) + "]"
            elif fit < 0 and bed_line[5] == "-":
                stop_exon = "[" + str(abs(fit)) + "]L" + stop_exon_dist_dict[fit]
            elif fit > 0 and bed_line[5] == "+":
                stop_exon = stop_exon_dist_dict[fit] + "S[" + str(abs(fit)) + "]"
            elif fit > 0 and bed_line[5] == "-":
                stop_exon = "[" + str(abs(fit)) + "]S" + stop_exon_dist_dict[fit]
            stop_fit = fit

        if distance_watcher_start > 100 and len(start_exon_dist_dict.keys()) > 0 and start_exon == 0:

            fit = int(min(list(start_exon_dist_dict.keys()), key=abs))

            # we have the same exon, have
            # to add the modifier to the same string
            if start_exon_dist_dict[fit] != tmp_exon:

                if fit > 0 and bed_line[5] == "-":
                    start_exon = start_exon_dist_dict[fit] + "L[" + str(abs(fit)) + "]"
                elif fit > 0 and bed_line[5] == "+":
                    start_exon = "[" + str(abs(fit)) + "]L" + start_exon_dist_dict[fit]
                elif fit < 0 and bed_line[5] == "+":
                    start_exon = "[" + str(abs(fit)) + "]S" + start_exon_dist_dict[fit]
                elif fit < 0 and bed_line[5] == "-":
                    start_exon = start_exon_dist_dict[fit] + "S[" + str(abs(fit)) + "]"

            else:
                # we directly modify the stop exon and add our modifier
                if fit > 0 and bed_line[5] == "-":
                    stop_exon = stop_exon + "L[" + str(abs(fit)) + "]"
                elif fit > 0 and bed_line[5] == "+":
                    stop_exon = "[" + str(abs(fit)) + "]L"  + stop_exon
                elif fit < 0 and bed_line[5] == "+":
                    stop_exon = "[" + str(abs(fit)) + "]S" + stop_exon
                elif fit < 0 and bed_line[5] == "-":
                    stop_exon = stop_exon + "S[" + str(abs(fit)) + "]"
                # will be unified in the next step
                start_exon = 0

            start_fit = fit

        # only remove start/stop part (i.e. the > 100bp error)
        # the other side might be okay ith a small error, e.g. 10bp,
        # so we do not want to remove that
        # also keep the modification if we did not find a better exon
        work_exon = ""

        if distance_watcher_start > 100 and start_exon != 0:
            work_exon = return_list[0]
            # first exon

            if bed_strand == "+":
                work_exon = re.sub("\[.*\][LS]{1}", "", work_exon)
            else:
                work_exon = re.sub("[LS]{1}\[.*\]", "", work_exon)

        if distance_watcher_stop > 100 and stop_exon != 0:

            if not work_exon:
                work_exon = return_list[-1]

            if bed_strand == "+":
                work_exon = re.sub("[LS]{1}\[.*\]", "", work_exon)
            else:
                work_exon = re.sub("\[.*\][LS]{1}", "", work_exon)

        # we looked for new 5' and 3'
        if distance_watcher_start > 100 and distance_watcher_stop > 100:

            # let's start with the start exon
            # we found both new exons, so the middle old exon MUST have no
            # mods anymore

            if start_exon != 0 and stop_exon != 0:
                changed_start = True
                changed_stop = True
                # only if we have at least 2 elements
                if len(return_list) >= 2:
                    del return_list[0]
                    del return_list[-1]
                else:
                    return_list.clear()
                    return_list.append(work_exon)

                return_list.insert(0, str(start_exon))
                return_list.append(str(stop_exon))

            # we got a new start but not stop
            elif start_exon != 0 and stop_exon == 0:
                changed_start = True
                del return_list[0]
                return_list.insert(0, str(start_exon))

            # we got a new stop but not start
            elif start_exon == 0 and stop_exon != 0:
                changed_stop = True
                del return_list[-1]
                return_list.append(str(stop_exon))

        # next case: only looked found a new start exon
        if distance_watcher_start > 100 and start_exon != 0 and not changed_start:
                changed_start = True
                # only if we have at least 2 elements
                if len(return_list) >= 2:
                    del return_list[0]
                    return_list.insert(0, start_exon)
                    # del return_list[-1]
                    # return_list.append(str(work_exon))
                else:
                    return_list.clear()
                    return_list.append(str(start_exon))
                    return_list.append(str(work_exon))

        # next case: only looked found a new stop exon
        if distance_watcher_stop > 100 and stop_exon != 0 and not changed_stop:

                changed_stop = True
                # only if we have at least 2 elements
                if len(return_list) >= 2:
                    del return_list[-1]
                    return_list.append(str(stop_exon))
                else:
                    return_list.clear()
                    return_list.append(str(work_exon))
                    return_list.append(str(stop_exon))

    # either still nothing, or, distance still > 1000 insert new exons
    if abs(start_fit) > 1000 and abs(stop_fit) > 1000:
        # still no good fit, let's assume we have new exons

        if bed_strand == "+":
            return_list[0] = re.sub("\[.*\][LS]{1}", "", return_list[0])
        else:
            return_list[0] = re.sub("[LS]{1}\[.*\]", "", return_list[0])

        if bed_strand == "+":
            return_list[-1] = re.sub("[LS]{1}\[.*\]", "", return_list[-1])
        else:
            return_list[-1] = re.sub("\[.*\][LS]{1}", "", return_list[-1])

        return_list.insert(0, "NE")
        return_list.append("NE")


    elif abs(start_fit) > 1000 and abs(stop_fit) < 1000:
        # create new exon for the start, as we did not find any good match
        # keep old okay exon and insert new start exon
        if len(return_list) >= 2:
            del return_list[0]

        if bed_strand == "+":
            return_list[0] = re.sub("\[.*\][LS]{1}", "", return_list[0])
        else:
            return_list[0] = re.sub("[LS]{1}\[.*\]", "", return_list[0])

        return_list.insert(0, "NE")

    elif abs(start_fit) < 1000 and abs(stop_fit) > 1000:
        # create new exon for the stop, as we did not find any good match
        # keep old okay exon and insert new end exon
        if len(return_list) >= 2:
            del return_list[-1]

        if bed_strand == "+":
            return_list[-1] = re.sub("[LS]{1}\[.*\]", "", return_list[-1])
        else:
            return_list[-1] = re.sub("\[.*\][LS]{1}", "", return_list[-1])
        return_list.append("NE")

    if bed_strand == "-":
        return_list.reverse()

    if single_intron_check:
        return "circ" + name[0] + "(NE)", \
            max(distance_watcher_start, distance_watcher_stop)
    else:
        return "circ"+name[0]+"(" + ",".join(return_list) + ")",\
            max(start_fit, stop_fit)


# process file function
def process_circrna_file(filename,
                         start=0,
                         stop=0,
                         coordinate_dict=None,
                         result_dict = None,
                         blast_db=None,
                         bedfile=None,
                         stats=None):

    with open(filename, 'r') as fp:

        lines = fp.readlines()[start:stop+1]

        for line in lines:

            line = line.strip().split("\t")

            # we got a circRNA with sequence
            circrna_id = line[1]
            circrna_sequence = line[2]

            # get coordinates from lookup table
            circrna_coordinates = coordinate_dict[circrna_id]

            # only partial sequence, nothing to do here
            if line[2] == "partial":
                stats[circrna_id] = {"coordinates": str(circrna_coordinates),
                                     "source": "PARTIAL"}
                continue

            # we do not look at intergeni circRNAs (yet)
            # TODO: add intergenic circRNA support
            if "intergenic" in line[1]:
                stats[circrna_id] = {"coordinates": str(circrna_coordinates),
                                     "source": "INTERGENIC"}
                continue

            # prepare TMP files for BLASTing
            # just plain write mode, not binary mode as this
            # pipeline is mainly intended for Linux systems
            fasta_query = tempfile.NamedTemporaryFile(mode="w")
            limit_text = tempfile.NamedTemporaryFile(mode="w")
            limit_bed_db = tempfile.NamedTemporaryFile(mode="w")
            limit_bed_query = tempfile.NamedTemporaryFile(mode="w")

            limit_bsl = tempfile.NamedTemporaryFile(mode="w")
            blast_output = tempfile.NamedTemporaryFile(mode="w")

            # FASTA header
            fasta_query.write(">" + circrna_id + "@" + circrna_coordinates + "\n")

            # sequence
            fasta_query.write(circrna_sequence)

            debug_file_name = "out/" +circrna_id+".bed"

            fasta_query.flush()

            # wrote query, now we have to prepare the limit file
            # for -seqidlist limiting

            # get gene name from circAtlas ID
            gene_name = circrna_id.split("-")[1].split("_")[0]

            # check we got any lines in our FASTA file for BLAST
            # if not, skip BLAST step

            # we also check that no sequences are in the seq list that are
            # massively out of circRNA range.
            # This might happen in case of repetitive intron sequences that
            # match dozens of times.

            circ_start = int(circrna_coordinates.split("|")[0].split(":")[1])
            circ_chr = circrna_coordinates.split("|")[0].split(":")[0].replace("chr","")
            circ_stop = int(circrna_coordinates.split("|")[1])

            # make a BED file from the potentially BLASTable regions
            # used for intersect to filter out far away regions
            os.system("egrep \"\s" +
                      gene_name +
                      "!\" " + bedfile +
                      " > "+limit_bed_db.name)

            # make sure data is written to file
            limit_bed_db.flush()

            # nothing in the BLAST prefilter, we can stop here already
            # this has to be fixed later with bedtools intersect if at all
            # possible
            if os.path.getsize(limit_bed_db.name) == 0:
                stats[circrna_id] = {"coordinates": str(circrna_coordinates),
                                     "source": "EMPTY_BLAST_DB"}
                continue

            debug_bed_file = open(debug_file_name, "w")

            # this is our second BED file just containing the circRNA
            # coordinates + 1000 flanking BP as buffer

            circrna_line = "\t".join([circ_chr,
                                             str(circ_start-200),
                                             str(circ_stop+200),
                                             circrna_id])

            limit_bed_query.write(circrna_line+"\n")

            limit_bed_query.flush()

            # create BEDtool instances
            bed_obj_query = BedTool(limit_bed_query.name)
            bed_obj_db = BedTool(limit_bed_db.name)

            # get exons & introns that somehow intersect with the circRNA
            result = bed_obj_db.intersect(bed_obj_query)

            row_count = 0

            bed_strand = ""

            for bed_line in io.StringIO(str(result)):
                content = bed_line.strip().split("\t")

                # write name of intersected sequences in limit file
                # print(content[3])
                limit_text.write(content[3]+"\n")
                # debug_bed_file.write("\t".join([content[0],
                #                                 content[1],
                #                                 content[2],
                #                                 content[3]+"_intersect",
                #                                 content[4],
                #                                 content[5]])+"\n")
                row_count += 1
                bed_strand = content[5]

            limit_text.flush()

            if row_count == 0:
                # TODO: this should be a branch out call for standard
                # bedtools interect to get exon numbers based on coordinates
                stats[circrna_id] = {"coordinates": str(circrna_coordinates),
                                     "source": "EMPTY_BED_INTERSECT"}
                continue

            # convert limit file so BLAST can read it
            os.system("blastdb_aliastool -seqid_file_in " +
                      limit_text.name +
                      " -seqid_file_out " +
                      limit_bsl.name)

            limit_bsl.flush()

            # start blast query for this circRNA
            run_blast_query(
                database=blast_db,
                fasta_query_file=fasta_query.name,
                limit_bsl_file=limit_bsl.name,
                blast_output_file=blast_output.name,
            )

            blast_output.flush()

            if os.path.getsize(blast_output.name) == 0:
                stats[circrna_id] = {"coordinates": str(circrna_coordinates),
                                     "source": "NO_BLAST_HIT"}
                continue

            # perform some alchemy to magically conjure the new IDs
            new_id = process_blast_results(file_name=blast_output.name,
                                           debug_bed_file=debug_bed_file,
                                           intersect_bed_file=result,
                                           bed_strand=bed_strand)

            circrna_line = "\t".join([circ_chr,
                                             str(circ_start),
                                             str(circ_stop),
                                             str(new_id[0])
                                      + " / " + circrna_id])

            debug_bed_file.write(circrna_line+"\n")

            # close TMP files in order get clean up
            # will be automatically deleted
            limit_text.close()
            limit_bsl.close()
            fasta_query.close()
            blast_output.close()
            debug_bed_file.close()

            # remove from circRNA dict and mark as done
            del coordinate_dict[circrna_id]

            result_dict[circrna_id] = {"new_id": new_id[0],
                                       "error": new_id[1],
                                       "coordinates": str(circrna_coordinates),
                                       "source": "BLAST_HIT"}

            stats[circrna_id] = {"coordinates": str(circrna_coordinates),
                                 "source": "BLAST_HIT"}

            sys.stdout.write("\r"+f'{len(result_dict):,}' +
                             " circRNAs processed, "+
                             f'{len(coordinate_dict):,}'+
                             " circRNAs remaining")
            sys.stdout.flush()


def process_remaining_circrnas(filename,
                               start=0,
                               stop=0,
                               coordinate_dict=None,
                               result_dict=None,
                               bedfile=None,
                               stats=None):

    with open(filename, 'r') as fp:

        lines = fp.readlines()[start:stop+1]

        for line in lines:

            line = line.strip().split("\t")

            circrna = line[1]

            # this circRNA has already been processed with the BLAST pipeline
            if circrna in result_dict:
                continue

            # fix for file header
            if circrna == "circpedia2":
                continue

            # these are hg19 only circRNAs that have no circatlas ID nor
            # hg38 coordinates.
            # those won't be used or converted (also for rn5 or mm9)
            if circrna == "-":
                stats[line[3]] = {"coordinates": str(line[3]),
                                     "source": "OLD_GENOME"}
                continue

            # check what kind of circRNA ID we have:
            # a) normal circatlas: xyz-Gene_AAAA
            # b) no id: chrX:start|stop

            # this is case a)
            if "-" in circrna:
                gene_name = circrna.split("-")[1].split("_")[0]
            # this is case b)
            else:
                gene_name = "NA"

            bed_items_to_intersect = tempfile.NamedTemporaryFile(mode="w")

            coordinates = coordinate_dict[circrna]

            circ_start = int(coordinate_dict[circrna].split("|")[0].split(":")[1])
            circ_chr = coordinate_dict[circrna].split("|")[0].split(":")[0].replace("chr", "")
            circ_stop = int(coordinate_dict[circrna].split("|")[1])

            circ_length = circ_stop - circ_start

            bedtools_string = "\t".join([circ_chr,
                                         str(circ_start - 1000),
                                         str(circ_stop + 1000),
                                         "exon_fix"]) + "\n"

            # get only parts of the BED file that correspond to the current gene
            if gene_name != "NA":
                os.system("egrep \"\s" +
                          gene_name +
                          "!\" " + bedfile +
                          " > " + bed_items_to_intersect.name)
                bed_items_to_intersect.flush()
                bed_obj_db = BedTool(bed_items_to_intersect.name)
            else:
                # well, no gene name, we have to run the full intersect...
                bed_obj_db = BedTool(bedfile)

            bed_obj_query = BedTool(bedtools_string, from_string=True)

            # sort by coordinates
            bed_obj_query = bed_obj_query.sort()

            # merge contained exons, i.e. if exon 5 is 500 bp and exon 3 if 100 bp
            # completely within exon 5 only keep exon 5
            bed_obj_query = bed_obj_query.merge(c="4", o="first")

            intersection = bed_obj_db.intersect(bed_obj_query)

            # debug_file_name = "out/" +circrna+".bed"
            # with open(debug_file_name, mode="w", ) as debug_file:
            #     debug_file.write(str(intersection))
            # debug_file.close()

            start_exon = 0
            stop_exon = 0

            start_fit = 0
            stop_fit = 0

            # in case on inexact matches, holds the best distances for each exon
            start_exon_dist_dict = {}
            stop_exon_dist_dict = {}
            exon_dict = {}

            bed_gene_name = ""

            bed_gene_name_list = []

            for bed_line in io.StringIO(str(intersection)):

                bed_line = bed_line.strip().split("\t")

                # get positions from the name string, NOT from the actual BED
                # positions as we know those already from above

                name = bed_line[3]

                bed_name_start_pos = int(name.split("!")[3])
                bed_name_stop_pos = int(name.split("!")[4])
                bed_name_type = name.split("!")[5]

                if bed_gene_name != name.split("!")[0]:
                    bed_gene_name_list.append(name.split("!")[0])
                    bed_gene_name = name.split("!")[0]

                # only if we have an exon
                if bed_name_type == "E":
                    bed_name_exon_num = name.split("!")[6]

                    # we need this to resolve gene names for unknown circRNAs
                    # based on the matched exons
                    exon_dict[bed_name_exon_num] = bed_gene_name
                elif bed_name_type == "G" and gene_name == "NA":
                    bed_gene_name_list.append(name.split("!")[0])

                # exon processing

                # direct hit for start
                if bed_name_type == "E" and bed_name_start_pos == circ_start:
                    start_exon = bed_name_exon_num
                    start_exon_dist_dict[0] = bed_name_exon_num

                elif bed_name_type == "E":
                    start_exon_dist_dict[
                        circ_start - bed_name_start_pos] = bed_name_exon_num

                # if not elif because we might have a single-exon circRNA

                # direct hit for stop
                if bed_name_type == "E" and bed_name_stop_pos == circ_stop:
                    stop_exon = bed_name_exon_num
                    # print("perfect end")
                    stop_exon_dist_dict[0] = bed_name_exon_num


                elif bed_name_type == "E":
                    stop_exon_dist_dict[
                        circ_stop - bed_name_stop_pos] = bed_name_exon_num

                # intron processing
                # because we did not hit any exons in the area

                if bed_name_type == "RI" and start_exon == 0:
                    # print("intron start match")
                    start_intron = "NE"

                if bed_name_type == "RI" and stop_exon == 0:
                    # print("intron stop match")
                    stop_intron = "NE"

            tmp_exon = ""

            min_val = 99999999

            # print(start_exon_dist_dict)
            # print(stop_exon_dist_dict)

            # check the abs distance to exon borders if start and stop are
            # on the same exon, i.e. make a single exon circRNA
            # we memorize the smallest combined absolute distance

            for start_key in start_exon_dist_dict:
                for stop_key in stop_exon_dist_dict:
                    if start_exon_dist_dict[start_key] == stop_exon_dist_dict[stop_key]:
                        if abs(start_key) + abs(stop_key) < min_val:
                            start_fit = start_key
                            stop_fit = stop_key
                            min_val = abs(start_key) + abs(stop_key)

            # we migh already found a perfect hit above
            # reset the start/stop fit to 0 to make sure we do not stick
            # to the probably worse fit if using only one exon

            if start_exon != 0:
                start_fit = 0

            if stop_exon != 0:
                stop_fit = 0

            # perfect hit yet, let's look at the suboptimal hits
            if stop_exon == 0 and len(stop_exon_dist_dict.keys()) > 0:

                if stop_fit != 0 and min(list(stop_exon_dist_dict.keys()), key=abs) < abs(stop_fit):
                    stop_fit = min(list(stop_exon_dist_dict.keys()), key=abs)


                # print("stop fit "+str(list(stop_exon_dist_dict.keys())))

                tmp_exon = stop_exon_dist_dict[stop_fit]

                if stop_fit < 0 and bed_line[5] == "+":
                    stop_exon = stop_exon_dist_dict[stop_fit] + "L[" + str(
                        abs(stop_fit)) + "]"
                elif stop_fit < 0 and bed_line[5] == "-":
                    stop_exon = "[" + str(abs(stop_fit)) + "]L" + \
                                stop_exon_dist_dict[stop_fit]
                elif stop_fit > 0 and bed_line[5] == "+":
                    stop_exon = stop_exon_dist_dict[stop_fit] + "S[" + str(
                        abs(stop_fit)) + "]"
                elif stop_fit > 0 and bed_line[5] == "-":
                    stop_exon = "[" + str(abs(stop_fit)) + "]S" + \
                                stop_exon_dist_dict[stop_fit]

            if start_exon == 0 and len(start_exon_dist_dict.keys()) > 0:

                # print("start fit " + str(start_fit))

                if start_fit != 0 and min(list(start_exon_dist_dict.keys()), key=abs) < abs(start_fit):
                    start_fit = min(list(start_exon_dist_dict.keys()), key=abs)

                    # we have the same exon, have
                # to add the modifier to the same string
                if start_exon_dist_dict[start_fit] != tmp_exon:

                    if start_fit > 0 and bed_line[5] == "-":
                        start_exon = start_exon_dist_dict[start_fit] + "L[" + str(
                            abs(start_fit)) + "]"
                    elif start_fit > 0 and bed_line[5] == "+":
                        start_exon = "[" + str(abs(start_fit)) + "]L" + \
                                     start_exon_dist_dict[start_fit]
                    elif start_fit < 0 and bed_line[5] == "+":
                        start_exon = "[" + str(abs(start_fit)) + "]S" + \
                                     start_exon_dist_dict[start_fit]
                    elif start_fit < 0 and bed_line[5] == "-":
                        start_exon = start_exon_dist_dict[start_fit] + "S[" + str(
                            abs(start_fit)) + "]"

                else:
                    # we directly modify the stop exon and add our modifier
                    if start_fit > 0 and bed_line[5] == "-":
                        stop_exon = stop_exon + "L[" + str(abs(start_fit)) + "]"
                    elif start_fit > 0 and bed_line[5] == "+":
                        stop_exon = "[" + str(abs(start_fit)) + "]L" + stop_exon
                    elif start_fit < 0 and bed_line[5] == "+":
                        stop_exon = "[" + str(abs(start_fit)) + "]S" + stop_exon
                    elif start_fit < 0 and bed_line[5] == "-":
                        stop_exon = stop_exon + "S[" + str(abs(start_fit)) + "]"

                        # will be unified in the next step
                        start_exon = 0

            # no exon matched, try introns
            if start_exon == 0 and stop_exon == 0:
                start_exon = "NE"

            if stop_exon == 0:
                stop_exon = "NE"

            if start_exon == stop_exon and start_exon != "NE":
                return_list = [start_exon]
            elif start_exon == stop_exon == "NE" and circ_length > 1000:
                return_list = [start_exon]
            elif start_exon != 0 and stop_exon != 0:
                return_list = [start_exon, stop_exon]
            else:
                return_list = [stop_exon]

            # either still nothing, or, distance still > 1000 insert new exons
            if abs(start_fit) > 1000 and abs(stop_fit) > 1000:
                # still no good fit, let's assume we have new exons
                return_list.clear()
                return_list.insert(0, "NE")
                return_list.append("NE")

            elif abs(start_fit) > 1000 and abs(stop_fit) < 1000:

                del return_list[0]
                return_list.insert(0, "NE")

            elif abs(start_fit) < 1000 and abs(stop_fit) > 1000:

                del return_list[-1]
                return_list.append("NE")

            if bed_line[5] == "-":
                return_list.reverse()

            # recover gene name if possible:
            if gene_name == "NA":

                if len(exon_dict) > 0:
                    bed_gene_name_list.clear()
                    for exon in exon_dict:
                        bed_gene_name_list.append(exon_dict[exon])

                gene_name = "|".join(list(set(bed_gene_name_list)))

            new_id = "circ" + gene_name + "(" + ",".join(return_list) + ")"

            # remove from circRNA dict and mark as done
            del coordinate_dict[circrna]

            result_dict[circrna] = {"new_id": new_id,
                                    "error": max(abs(start_fit), abs(stop_fit)),
                                    "coordinates": coordinates,
                                    "source": "BEDTOOLS"}

            stats[line[3]] = {"coordinates": str(coordinates),
                              "source": "BEDTOOLS"}

            sys.stdout.write("\r"+f'{len(result_dict):,}' +
                             " circRNAs processed, "+
                             f'{len(coordinate_dict):,}'+
                             " circRNAs remaining")
            sys.stdout.flush()


parser = argparse.ArgumentParser(
    prog="convert_circrna_names",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    fromfile_prefix_chars="@",
    description="Maps CircAtlas circRNA IDs to the Chen et al. 2023 naming "
                "scheme. "
                "Uses data prepared by naming_conversion_data_preparation\n"
                "\n"
                "Version 0.0.2\n"
                "\n"
                "https://github.com/jakobilab/circhemy\n"
                "https://jakobilab.org\n"
                "tjakobi@arizona.edu",

    usage=""" naming_conversion_data_preparation [<args>]"""
)

group = parser.add_argument_group("output parameters")

group.add_argument("-c",
                   "--circatlas",
                   dest="circatlas",
                   help="Path to the CircAtlas data table",
                   required=True
                   )

group.add_argument("-s",
                   "--sequence",
                   dest="circatlas_sequences",
                   help="Path to the CircAtlas sequences file",
                   required=True
                   )


group.add_argument("-o",
                   "--output",
                   dest="output_file",
                   default="./circrna_mapping.csv",
                   help="The output file for the circRNA name mapping"
                   )

group.add_argument("-b",
                   "--blast",
                   dest="blast_db",
                   help="Path to BLAST database",
                   required=True
                   )

group.add_argument("-B",
                   "--bedfile",
                   dest="bedfile",
                   help="Path to the BED file containing exons, "
                        "introns, and genes.",
                   required=True
                   )

group.add_argument("-t",
                   "--threads",
                   dest="cpu_threads",
                   help="Number of CPUs to use for multiprocessing",
                   default=4,
                   type=int
                   )

group.add_argument("-C",
                   "--chunksize",
                   dest="chunk_size",
                   help="Number of lines to process in batch per thread",
                   default=10000,
                   type=int
                   )


args = parser.parse_args()


if __name__ == "__main__":

    manager = mp.Manager()

    # generate dictionary to resolve circatlas IDs to coordinates
    coordinate_dict = manager.dict()

    # hold circRNA ID (circAtlas) as key and novel naming scheme ID as value
    result_dict = manager.dict()

    # holds stats for each circRNA, i.e. found via BLAST, bedtools, or why no
    # it was excluded
    stats_dict = manager.dict()

    total_circrna_count = 0
    print("Preparing circAtlas ID -> coordinate lookup table")
    with open(args.circatlas) as f:
        next(f)
        # skip header
        for line in f:
            line = line.split("\t")
            if line[1] == "-":
                # this is an hg19/rn5/mm9-only circRNA
                # we skip these
                continue
            coordinate_dict[line[1]] = line[2]
            total_circrna_count += 1

    cpu_count = args.cpu_threads

    # get file size and set chuck size
    line_chunk_size = args.chunk_size

    print(f'{total_circrna_count:,}' + " circRNAs to process")
    print("Starting phase 1: assignment of exons for"
          " circRNAs with full sequence information")

    # determine if it needs to be split
    if total_circrna_count > line_chunk_size:

        # create pool, initialize chunk start location (cursor)
        pool = mp.Pool(cpu_count)
        cursor = 1
        results = []

        # for every chunk in the file...
        for chunk in range((total_circrna_count // line_chunk_size) + 1):

            # determine where the chunk ends, is it the last one?
            if cursor + line_chunk_size > total_circrna_count:
                end = total_circrna_count
            else:
                end = cursor + line_chunk_size

            # add chunk to process pool, save reference to get results
            proc = pool.apply_async(process_circrna_file,
                                    args=[args.circatlas_sequences,
                                          cursor,
                                          end,
                                          coordinate_dict,
                                          result_dict,
                                          args.blast_db,
                                          args.bedfile,
                                          stats_dict])
            # setup next chunk
            cursor = end + 1

        # close and wait for pool to finish
        for proc in results:
            processfile_result = proc.get()

        pool.close()
        pool.join()

    else:
        process_circrna_file(filename=args.circatlas_sequences,
                             start=1,
                             stop=total_circrna_count,
                             coordinate_dict=coordinate_dict,
                             result_dict=result_dict,
                             blast_db=args.blast_db,
                             bedfile=args.bedfile,
                             stats=stats_dict)

    # here starts fallback code for
    #  - circRNAs without BLAST results
    #  - intergenic circRNAs
    #  - circRNAs without sequence in the sequence table (or partial)

    print("")
    print("Done")
    print("Starting phase 2: assignment of exons for"
          " circRNAs without sequence information")

    # determine if it needs to be split
    if total_circrna_count > line_chunk_size:

        # create pool, initialize chunk start location (cursor)
        pool = mp.Pool(cpu_count)
        cursor = 1
        results = []

        # for every chunk in the file...
        for chunk in range((total_circrna_count // line_chunk_size) + 1):

            # determine where the chunk ends, is it the last one?
            if cursor + line_chunk_size > total_circrna_count:
                end = total_circrna_count
            else:
                end = cursor + line_chunk_size

            # add chunk to process pool, save reference to get results
            proc = pool.apply_async(process_remaining_circrnas,
                                    args=[args.circatlas,
                                          cursor,
                                          end,
                                          coordinate_dict,
                                          result_dict,
                                          args.bedfile,
                                          stats_dict])
            # results.append(proc)

            # setup next chunk
            cursor = end + 1

        # close and wait for pool to finish
        pool.close()
        pool.join()

    else:
        process_remaining_circrnas(filename=args.circatlas,
                                   start=1,
                                   stop=total_circrna_count,
                                   coordinate_dict=coordinate_dict,
                                   result_dict=result_dict,
                                   bedfile=args.bedfile,
                                   stats=stats_dict)

    print("")
    print("Done")
    print("Writing output file")
    with open(args.output_file, mode="w", ) as out_file:
        for circrna in result_dict:
            out_file.write("\t".join([circrna,
                                      result_dict[circrna]['new_id'],
                                      str(result_dict[circrna]['error']),
                                      result_dict[circrna]['coordinates'],
                                      result_dict[circrna]['source']
                                      ]) + "\n")
    print("Done")

    tmp_dict = {'NO_BLAST_HIT': 0,
                'BLAST_HIT': 0,
                'OLD_GENOME': 0,
                'BEDTOOLS': 0,
                'PARTIAL': 0,
                'INTERGENIC': 0,
                'EMPTY_BLAST_DB': 0,
                'EMPTY_BED_INTERSECT': 0
                }

    # summarize stats data
    for circrna in stats_dict:
        tmp_dict[stats_dict[circrna]['source']] += 1

    print("Stats:")
    print("-----------------")
    # print out data
    for endpoint in tmp_dict:
        print(endpoint + ": " + str(tmp_dict[endpoint]))