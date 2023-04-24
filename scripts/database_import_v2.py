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
import bz2
import time
from datetime import datetime

# own util functions
import circhemy.common.util as common

# create util instance for the web app
util = common.Util


def process_input_data(input_data: str, db_connection):

    now = int(time.time())

    print("Input file:\t" + input_data)

    print("SQLite3 DB:\t"+str(db_connection))

    print_time = datetime.utcfromtimestamp(now).strftime('%Y-%m-%d')

    db_version_from_time = datetime.utcfromtimestamp(now).strftime('%Y.%m')

    print("Circhemy DB time stamp:\t"+str(print_time))

    print("Circhemy DB version:\t"+str(db_version_from_time))

    db_connection.execute("BEGIN TRANSACTION")

    sqlite_insert_query = """INSERT INTO circhemy_db_info
                      (Version, Date)
                      VALUES (?, ?);"""

    # insert corresponding log netry
    db_connection.execute(sqlite_insert_query, [db_version_from_time, now])

    db_connection.execute("COMMIT TRANSACTION")

    generated_db_id = db_connection.lastrowid

    print("Database version is now: " + str(generated_db_id))

    source_file = bz2.open(input_data, "rt")
    count = 0

    for line in source_file:

        # pseudocode:

        # go through each line of the input

        # check if chr:start:stop:csnv1 is already in db

        # if not:
        # 1) insert line in circhemy table
        # 2) insert line in circhemy_log with current db version as added

        # if yes:
        # todo: find out what field id different
        # insert line in circhemy_log with current db version as modified

        # add line in db info with current version

        # start transaction
        db_connection.execute("BEGIN TRANSACTION")

        # get all fields of the file
        line = line.rstrip()
        entry = line.split("\t")

        # replace NAs and empty fields with None
        # None is translated into NULL fot SQLite
        entry = [None if x == 'NA' or x == '' else x for x in entry]

        sqlite_insert_query = """INSERT INTO circhemy
                          (Species, Gene, Description, ENSEMBL, Entrez, circBase, CircAtlas2, circRNADb, deepBase2,
                           Circpedia2, circBank, riboCIRC, exoRBase2, Arraystar, CSNv1, Chr, Start, Stop, Strand,
                           Genome, Pubmed)
                          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"""

        db_connection.execute(sqlite_insert_query, entry)

        # get ID
        circrna_id = int(db_connection.lastrowid)

        if circrna_id % 10000 == 0:

            print(str(circrna_id) + " data rows processed.")

        sqlite_insert_query = """INSERT INTO circhemy_log
                          (CircRNA_ID, Action, DB_ID)
                          VALUES (?, ?, ?);"""

        # insert corresponding log entry
        # 1 for action == addition for now
        # more numbers to follow
        db_connection.execute(sqlite_insert_query, [circrna_id, 1, generated_db_id])

        # finish transaction
        db_connection.execute("COMMIT TRANSACTION")

    source_file.close()


parser = argparse.ArgumentParser(
    prog="database_import_v2.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    fromfile_prefix_chars="@",
    description="Inserts data prepared with build_sqlite_db.R into SQLite database\n"
                "\n"
                "Version 0.0.1\n"
                "\n"
                "https://github.com/jakobilab/circhemy\n"
                "https://jakobilab.org\n"
                "tjakobi@arizona.edu",

    usage=""" database_import_v2 [<args>]"""
)

group = parser.add_argument_group("output parameters")

group.add_argument("-d",
                   "--database",
                   dest="database",
                   default="../circhemy/data/circhemy.sqlite3",
                   help="The SQLite3 database file",
                   required=True
                   )

group.add_argument("-i",
                   "--input",
                   dest="input",
                   help="Input data, CSV-formatted",
                   required=True
                   )

args = parser.parse_args()

util.setup_database(util, args.database)

sqlite_db = util.db_connection.execute("pragma query_only = OFF;")

process_input_data(args.input, sqlite_db)
