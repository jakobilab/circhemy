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
import os
import time
from datetime import datetime

# own util functions
import circhemy.common.util as common

# create util instance for the web app
util = common.Util

keyword_db_columns = [
    "CircRNA_ID",
    "CSNv1",
    "Gene",
    "ENSEMBL",
    "circBase",
    "CircAtlas2",
    "circRNADb",
    "circBank",
    "deepBase2",
    "Circpedia2",
    "riboCIRC",
    "exoRBase2",
    "Arraystar"
]

max_rows = 50000

filename_prefix = "sitemap"

filename_suffix = ".xml"

sitemap_base = "https://circhemy.jakobilab.org/sitemap/"


def make_filename(directory: str, counter: int):

    return directory + filename_prefix + "_" + f"{counter:04d}" + "" + filename_suffix


def generate_sitemap(output_dir: str, db_connection, base_url: str):

    now = int(time.time())

    print("Output directory:\t" + output_dir)

    print("SQLite3 DB:\t"+str(db_connection))

    lastmod_time = datetime.utcfromtimestamp(now).strftime('%Y-%m-%dT%H:%M:%SZ')

    print("Circhemy DB time stamp:\t"+str(lastmod_time))

    # will be used to count the max number of entries
    per_file_counter = 0
    sitemap_counter = 1

    overall_counter = 0
    circrna_counter = 0

    # sql = "SELECT CircRNA_ID FROM circhemy LIMIT 0,100000"
    sql = "SELECT CircRNA_ID FROM circhemy"

    # yes, we could also iterate over all IDs with range,
    # but maybe at some point there might be deleted rows
    # thus we start directly with the correct list
    # db_connection.row_factory = lambda cursor, row: row[0]
    
    circrna_id_list = db_connection.execute(sql).fetchall()

    output_file = open(make_filename(output_dir, sitemap_counter), "w")

    # print(circrna_id_list)

    # write header for first file
    print('<?xml version="1.0" encoding="UTF-8"?>', file=output_file)
    print('<urlset xmln ="http://www.sitemaps.org/schemas/sitemap/0.9">', file=output_file)

    # this is the main loop
    for circrna_id in circrna_id_list:

        sql = "SELECT " + ",".join(keyword_db_columns) + " FROM circhemy WHERE CircRNA_ID = ?"
        circrna_data = db_connection.execute(sql, (str(circrna_id[0]),)).fetchall()

        if circrna_counter % 100000 == 0:
            print(f'{circrna_counter:,}' + " circRNA IDs processed.")
        circrna_counter += 1

        for current_id in range(0, len(keyword_db_columns)):

            # make sure entry is not "None"
            if circrna_data[0][current_id]:
                if per_file_counter == max_rows:

                    # "close" last file with /urlset
                    print("</urlset>", file=output_file)
                    output_file.flush()

                    # gzip file

                    os.system("gzip -f " + make_filename(output_dir, sitemap_counter))

                    # increase file counter
                    sitemap_counter += 1

                    # reset entry counter
                    overall_counter += per_file_counter
                    per_file_counter = 0

                    # open a new output file
                    print(f'{max_rows:,}' + " rows reached, switching to new file " + make_filename(output_dir, sitemap_counter))

                    output_file = open(make_filename(output_dir, sitemap_counter), "w")

                    # write new header
                    print('<?xml version="1.0" encoding="UTF-8"?>', file=output_file)
                    print('<urlset xmln ="http://www.sitemaps.org/schemas/sitemap/0.9">', file=output_file)

                print("\t<url>", file=output_file)
                print("\t\t<loc>" + base_url + str(circrna_data[0][current_id]) + "</loc>", file=output_file)
                print("\t\t<lastmod>" + lastmod_time + "</lastmod>", file=output_file)
                print("\t</url>", file=output_file)
                per_file_counter += 1

    # close final file with urlset
    print("</urlset>", file=output_file)
    output_file.flush()

    # gzip last file
    os.system("gzip -f " + make_filename(output_dir, sitemap_counter))

    overall_counter += per_file_counter
    print("Done.")
    print(f'{overall_counter:,}' + " locations written.")

    output_file = open(output_dir + "/" + "sitemap_index.xml", "w")

    # print(circrna_id_list)

    # write header for first file
    print('<?xml version="1.0" encoding="UTF-8"?>', file=output_file)
    print('<sitemapindex xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">', file=output_file)

    tmp_counter = 1

    # now we have to create the index file
    for sitemap in range(0, sitemap_counter):
        print("\t<sitemap>", file=output_file)
        print("\t\t<loc>" + make_filename(sitemap_base, tmp_counter) + ".gz</loc>", file=output_file)
        print("\t\t<lastmod>" + lastmod_time + "</lastmod>", file=output_file)
        print("\t</sitemap>", file=output_file)

        tmp_counter += 1
    print('</sitemapindex>', file=output_file)
    output_file.flush()

    os.system("gzip -f " + output_dir + "/" + "sitemap_index.xml")


parser = argparse.ArgumentParser(
    prog="generate_circhemy_sitemap.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    fromfile_prefix_chars="@",
    description="Generates a sitemap with all circRNA IDs for circhemy web interface SEO\n"
                "\n"
                "Version 0.0.1\n"
                "\n"
                "https://github.com/jakobilab/circhemy\n"
                "https://jakobilab.org\n"
                "tjakobi@arizona.edu",

    usage=""" generate_circhemy_sitemap.py [<args>]"""
)


group = parser.add_argument_group("input parameters")

group.add_argument("-d",
                   "--database",
                   dest="database",
                   default="../circhemy/data/circhemy.sqlite3",
                   help="The SQLite3 database file",
                   required=True
                   )

group.add_argument("-o",
                   "--output_dir",
                   dest="output",
                   help="Output directory, will hold XML sitemaps",
                   required=True
                   )

group.add_argument("-b",
                   "--base_url",
                   dest="url",
                   help="Base URL for the circhemy installation",
                   required=True
                   )

args = parser.parse_args()

util.setup_database(util, args.database)

sqlite_db = util.db_connection.execute("pragma query_only = ON;")

generate_sitemap(args.output, sqlite_db, args.url)
