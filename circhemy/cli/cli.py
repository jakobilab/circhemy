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
import sys
import circhemy.common.util as common

util = common.Util


def main():
    parser = argparse.ArgumentParser(
        prog=util.program_name,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        fromfile_prefix_chars="@",
        description="Convert circular RNA identifiers from different databases.\n"
                    "\n"
                    "Version " + util.software_version + "\n"
                    "\n"
                    "https://github.com/jakobilab/TOOLNAME\n"
                    "https://jakobilab.org\n"
                    "tjakobi@arizona.edu",

        usage=util.program_name+""" [-V] <command> [<args>]
    
        Available commands:
    
           convert: convert circRNA IDs
           query:   query local circRNA database
        """)
    parser.add_argument("command", help="Command to run")

    parser.add_argument("-V",
                        "--version",
                        action="version",
                        version=util.software_version
                        )

    # parse_args defaults to [1:] for args, but you need to
    # exclude the rest of the args too, or validation will fail

    args = parser.parse_args(sys.argv[1:2])

    if args.command == "convert":

        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            fromfile_prefix_chars="@",
        )
        group = parser.add_argument_group("input parameters")

        group.add_argument("-q",
                           dest="query_data",
                           nargs="+",
                           help="file with IDs to read. One ID per line; "
                                "Use -q STDIN for STDIN direct input",
                           required=True
                           )

        group.add_argument("-i",
                           dest="input_field",
                           help="type of input circular RNA ID, e.g. circBase",
                           choices=util.db_columns,
                           required=True
                           )

        group.add_argument("-o",
                           dest="output_fields",
                           help="desired output fields; "
                                "for multiple fields use space-separated "
                                "list of field names",
                           nargs="+",
                           required=True
                           )

        group = parser.add_argument_group("output parameters")

        group.add_argument("-O",
                           dest="output_file",
                           help="output file location; default: STDOUT",
                           default="STDOUT"
                           )

        group.add_argument("-S",
                           dest="separator_char",
                           help="specify the separator character for output; "
                                "default: tab (\\t)",
                           default="\t"
                           )

        group.add_argument("-E",
                           dest="empty_char",
                           help="specify the placeholder for "
                                "empty database fields; "
                                "default: NA",
                           default="NA"
                           )

        args = parser.parse_args(sys.argv[2:])

        # done with CLI parsing

        # make sure we only work on sanitized field names to minimize SQL errors
        util.check_output_field_names(util, args.output_fields)
        util.check_input_field_name(util, args.input_field)

        # running in STDIN mode, convert data for use
        if args.query_data == ["STDIN"]:

            # tmp list to move STDIN to list
            stdin = []

            for line in sys.stdin:
                if 'Exit' == line.rstrip():
                    break
                stdin.append(line.rstrip())
            args.query_data = stdin

            # done with STDIN preprocessing

        # setup db, get cursor
        util.setup_database(util, util.database_location)

        output = util.run_simple_select_query(util,
                                              args.output_fields,
                                              args.query_data,
                                              args.input_field)

        # process output

        sql_output = util.process_sql_output(output,
                                             seperator=args.separator_char,
                                             empty_char=args.empty_char)

        # default output to console via STDOUT
        if args.output_file == "STDOUT":
            print(sql_output, end='')
        # user specified file output, try to write file
        else:
            try:
                with open(args.output_file, 'w') as f:
                    f.write(sql_output)
            except FileNotFoundError:
                print("Output file" + args.output_file + " could not be created")
                exit(-1)

        # done with main program

        # close db connection
        util.db_connection.close()

    elif args.command == "query":

        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            fromfile_prefix_chars="@",
        )
        group = parser.add_argument_group("database queries")

        group.add_argument("-c",
                           dest="circbase_query",
                           help="specify circbase-related query. "
                                "Use -c keyword for exact search or "
                                "-c *keyword for partial matches."
                           )

        group.add_argument("-a",
                           dest="circatlas_query",
                           help="specify circatlas-related query. "
                                "Use -a keyword for exact search or "
                                "-a *keyword for partial matches."
                           )

        group.add_argument("-d",
                           dest="deepbase2_query",
                           help="specify deepbase2-related query. "
                                "Use -d keyword for exact search or "
                                "-d *keyword for partial matches."
                           )

        group.add_argument("-e",
                           dest="circpedia2_query",
                           help="specify circpedia2-related query. "
                                "Use -e keyword for exact search or "
                                "-e *keyword for partial matches."
                           )

        group.add_argument("-b",
                           dest="circbank_query",
                           help="specify circbank-related query. "
                                "Use -b keyword for exact search or "
                                "-b *keyword for partial matches."
                           )

        group.add_argument("-m",
                           dest="arraystar_query",
                           help="specify arraystar-related query. "
                                "Use -m keyword for exact search or "
                                "-m *keyword for partial matches."
                           )

        group.add_argument("-r",
                           dest="circrnadb_query",
                           help="specify circrnadb-related query. "
                                "Use -r keyword for exact search or "
                                "-r *keyword for partial matches."
                           )

        group = parser.add_argument_group("species & genome build queries")

        group.add_argument("-s",
                           dest="species_query",
                           help="specify species name",
                           choices=util.database_species_list
                           )

        group.add_argument("-g",
                           dest="genome_query",
                           help="specify genome build",
                           choices=util.database_genome_list
                           )

        group = parser.add_argument_group("genomic location & gene queries")

        group.add_argument("-C",
                           dest="chr_query",
                           help="specify chromosome-related query. "
                                "Use -h keyword for exact search or "
                                "-h *keyword for partial matches."
                           )

        group.add_argument("-t",
                           dest="start_query",
                           help="specify start-related query. "
                                "Use -t keyword for exact search or "
                                "-t *keyword for partial matches."
                           )

        group.add_argument("-T",
                           dest="stop_query",
                           help="specify stop-related query. "
                                "Use -T keyword for exact search or "
                                "-T *keyword for partial matches."
                           )

        group.add_argument("-G",
                           dest="gene_query",
                           help="specify gene-related query. "
                                "Use -G keyword for exact search or "
                                "-G *keyword for partial matches."
                           )

        group = parser.add_argument_group("output parameters")

        group.add_argument("-o",
                               dest="output_fields",
                               help="desired output fields; "
                                    "for multiple fields use space-separated "
                                    "list of field names",
                               nargs="+",
                               required=True
                               )

        group.add_argument("-O",
                           dest="output_file",
                           help="output file location; default: STDOUT",
                           default="STDOUT"
                           )

        group.add_argument("-S",
                           dest="separator_char",
                           help="specify the separator character for output; "
                                "default: tab (\\t)",
                           default="\t"
                           )

        group.add_argument("-E",
                           dest="empty_char",
                           help="specify the placeholder for "
                                "empty database fields; "
                                "default: NA",
                           default="NA"
                           )

        args = parser.parse_args(sys.argv[2:])

        # done with CLI parsing

        # make sure we only work on sanitized field names to minimize SQL errors
        util.check_output_field_names(util, args.output_fields)
        # util.check_input_field_name(util, args.input_field)

        # running in STDIN mode, convert data for use
        # if args.query_data == ["STDIN"]:
        #
        #     # tmp list to move STDIN to list
        #     stdin = []
        #
        #     for line in sys.stdin:
        #         if 'Exit' == line.rstrip():
        #             break
        #         stdin.append(line.rstrip())
        #     args.query_data = stdin
        #
        #     # done with STDIN preprocessing

        # setup db, get cursor
        util.setup_database(util, util.database_location)

        input_dict = dict(circbase=args.circbase_query,
                          circatlas=args.circatlas_query,
                          deepbase2=args.deepbase2_query,
                          circpedia2=args.circpedia2_query,
                          circbank=args.circbank_query,
                          arraystar=args.arraystar_query,
                          circrnadb=args.circrnadb_query,
                          species=args.species_query,
                          genome=args.genome_query,
                          chr=args.chr_query,
                          start=args.start_query,
                          stop=args.stop_query,
                          gene=args.gene_query
                          )

        # print(input_dict)

        sql_query = " "

        additional_output_fields = list()

        for query_item in input_dict:
            if input_dict[query_item]:

                additional_output_fields.append(query_item)

                if str(input_dict[query_item]).startswith("*"):

                    tmp = str(input_dict[query_item]).replace("*", "", 1)
                    sql_query += query_item + " LIKE \"%" \
                                 + tmp + "%\" AND "

                else:

                    sql_query += query_item + " == \"" \
                                 + input_dict[query_item] + "\" AND "

        # remove last AND from query
        sql_query = sql_query[:-4]

        output = util.run_keyword_select_query(util,
                                              args.output_fields +
                                               additional_output_fields,
                                              sql_query)

        # process output

        sql_output = util.process_sql_output(output,
                                             seperator=args.separator_char,
                                             empty_char=args.empty_char)

        # default output to console via STDOUT
        if args.output_file == "STDOUT":
            print(sql_output, end='')
        # user specified file output, try to write file
        else:
            try:
                with open(args.output_file, 'w') as f:
                    f.write(sql_output)
            except FileNotFoundError:
                print("Output file" + args.output_file + " could not be created")
                exit(-1)

        # done with main program

        # close db connection
        util.db_connection.close()

    else:
        print("Unknown command:", args.command)
        exit(-1)

