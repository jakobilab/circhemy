# Copyright (C) 2022 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either self.version 3 of the License, or
# (at your option) any later self.version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import os
import sqlite3
import circhemy


class Util(object):

    # global settings
    software_version = "0.0.3-dev"

    database_version = "2023.04"

    program_name = "circhemy"

    program_name_long = "The circhemy circRNA ID database"

    support_email = program_name+"@jakobilab.org"

    support_web = "https://github.com/jakobilab/" + program_name + "/issues/new"

    news_url = "https://redmine.jakobilab.org/projects/circhemy/" \
               "news.atom?key=c616b9fb231445ac4ca65d94db1d3207382798e9"

    database_location = circhemy.__path__[0]+"/data/circhemy.sqlite3"

    database_table_name = "circhemy"

    database_species_list = ["homo_sapiens",
                             "mus_musculus",
                             "rattus_norvegicus"]

    database_genome_list = ["hg19",
                            "hg38",
                            "mm9",
                            "mm10",
                            "rn5",
                            "rn6"]

    select_db_columns = [
        "CSNv1",
        "Gene",
        "ENSEMBL",
        "Entrez",
        "Description",
        "circBase",
        # "circBase_alt",
        "CircAtlas2",
        "circRNADb",
        "circBank",
        "deepBase2",
        "Circpedia2",
        "riboCIRC",
        "exoRBase2",
        "Arraystar",
        "Pubmed"
    ]

    active_db_columns = {
        "CSNv1": True,
        "Gene": True,
        "ENSEMBL": False,
        "Entrez": False,
        "Description": False,
        "circBase": True,
        # "circBase_alt": True,
        "CircAtlas2": True,
        "circRNADb": True,
        "circBank": True,
        "deepBase2": True,
        "Circpedia2": True,
        "riboCIRC": True,
        "exoRBase2": True,
        "Arraystar": True,
        "Pubmed": False
    }

    all_db_columns = [
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
    ]

    db_columns = ["Species"] + select_db_columns + ["Chr",
                                                    "Start",
                                                    "Stop",
                                                    "Strand",
                                                    "Genome"]

    external_db_urls = {"circBase": "http://www.circbase.org/cgi-bin/singlerecord.cgi?id=",
                        "circBase_alt": "http://www.circbase.org/cgi-bin/singlerecord.cgi?id=",
                        "CircAtlas2": "http://159.226.67.237:8080/new/circ_detail.php?ID=",
                        "Circpedia2": "/circrna/",
                        "circBank": "http://www.circbank.cn/infoCirc.html?id=",
                        "deepBase2": "https://rna.sysu.edu.cn/deepbase3/subpages/ViewDetail_circRNA.php?spe=hg19&name=",
                        "Arraystar": "/circrna/",
                        "circRNADb": "/circrna/",
                        "riboCIRC": "http://www.ribocirc.com/rna_detail.php?dependent=Condition-independent&circ_id=",
                        "Circ2Disease": "http://bioinformatics.zju.edu.cn/Circ2Disease/browse_circRNA_result.php?circRNA=",
                        "Pubmed": "https://pubmed.ncbi.nlm.nih.gov/",
                        "exoRBase2": "http://www.exorbase.org/exoRBaseV2/detail/detailInfo?kind=circRNA&id=",
                        "Chr": "",
                        "Start": "",
                        "Stop": "",
                        "Strand": "",
                        "Coordinates": "",
                        "Unspliced length": "",
                        "Genome": "",
                        "CSNv1": "https://genome.ucsc.edu/cgi-bin/hgTracks?",
                        "Genome-Browser": "https://genome.ucsc.edu/cgi-bin/hgTracks?",
                        "ENSEMBL": "http://www.ensembl.org/id/",
                        "Entrez": "https://www.ncbi.nlm.nih.gov/gene/",
                        "Description": "",
                        "Species": "",
                        "Gene": "",
                        "Stable circhemy database ID": "",
                        "CircRNA_ID":  "/circrna/"
                        }

    db_action_codes = {
        1: "Added",
        2: "Changed",
        3: "Deleted"
    }

    db_connection = ""
    db_cursor = ""

    def check_input_field_name(self, field):
        if field not in self.db_columns:
            print(field + " is not a valid input field name")
            exit(-1)
        return

    def check_output_field_names(self, field_list):
        for field in field_list:
            if field not in self.db_columns:
                print(field + " is not a valid output field name")
                exit(-1)
        return

    def prepare_coordinates(self, coord_list):

        # check line by line
        fixed_coord_list = []

        for line in coord_list:
            if ":" in line and "|" in line:
                fixed_coord_list.append(line)
            elif "\t" in line:
                tmp = line.split("\t")
                if len(tmp) == 3:
                    fixed_coord_list.append(tmp[0] +
                                            ":" + tmp[1] +
                                            "|" + tmp[2])

        return fixed_coord_list

    def check_input_return_found_circ_number(self, query_data, input_field):

        # this is a special case, we treat "Coordinates" as some kind of meta input
        # we break the input into chr, start and stop for the following query
        if input_field == "Coordinates":
            # build SQL string
            coords = self.prepare_coordinates(self, query_data)

            sql = "SELECT count(distinct(Chr || ':' || " \
                  "Start || '|' || Stop)) as Coordinates FROM " +\
                  self.database_table_name + \
                  " WHERE Chr || ':' || " \
                  "Start || '|' || Stop in ({seq})".format(
                    seq=','.join(['?'] * len(coords)))
            sql_output = self.db_cursor.execute(sql, coords).fetchall()

        else:
            # build SQL string
            sql = "SELECT count(distinct("+input_field+")) FROM " +\
                  self.database_table_name + \
                  " WHERE " + input_field + " in ({seq})".format(
                    seq=','.join(['?'] * len(query_data)))
            sql_output = self.db_cursor.execute(sql, query_data).fetchall()

        # return ratio (0->1)
        return sql_output[0][0]/len(query_data), sql_output[0][0]

    def database_stats(self):

        sql = "SELECT count() FROM "+self.database_table_name+";"

        dbsize = self.db_cursor.execute(sql).fetchall()[0][0]

        sql = "SELECT count(Genome), Genome, Species FROM "+\
              self.database_table_name+" group by Genome"

        sql_output = self.db_cursor.execute(sql).fetchall()

        chart_dict = {
#            'title': {'text': f"{int(dbsize):,}" + " CircRNAs in Database"},
            'title': {'text': ""},
            'chart': {'type': 'bar'},
            'yAxis': {'title': {'text': '# circRNAs'}},
            'xAxis': {'title': {'text': ''}, 'categories': [],
                      'labels': {'enabled': False}, 'tickLength': 0,
                      'tickInterval': 0},
            'series': [
                {'name': sql_output[0][1], 'data': [int(sql_output[0][0])]},
                {'name': sql_output[1][1], 'data': [int(sql_output[1][0])]},
                {'name': sql_output[2][1], 'data': [int(sql_output[2][0])]},
                {'name': sql_output[3][1], 'data': [int(sql_output[3][0])]},
                {'name': sql_output[4][1], 'data': [int(sql_output[4][0])]},
                {'name': sql_output[5][1], 'data': [int(sql_output[5][0])]}
            ],
        }

        chart2_dict = {
   #         'title': {'text': "Distribution of CircRNA Database IDs"},
            'title': {'text': ""},

            'chart': {'type': 'bar'},
            'yAxis': {'title': {'text': '# circRNA IDs'}, 'type': 'logarithmic'},
            'xAxis': {'title': {'text': ''}, 'categories': [],
                      'labels': {'enabled': False}, 'tickLength': 0,
                      'tickInterval': 0},
            'series': [
                {'name': 'circBase', 'data': [130248]},
                {'name': 'CircAtlas', 'data': [1996471]},
                {'name': 'Circpedia2', 'data': [399143]},
                {'name': 'CircBank', 'data': [127174]},
                {'name': 'Deepbase2', 'data': [200845]},
                {'name': 'Arraystar', 'data': [54394]},
                {'name': 'CircRNA DB', 'data': [59621]},
            ],
        }

        return chart_dict, dbsize, chart2_dict

    def setup_database(self, database):

        # check if there is a .bz2 version of the database
        # if yes, this is the first time circhemy runs
        # we have to unpack it one time

        if os.path.isfile(database+".bz2"):
            print("This is the first run of circhemy, local SQLite3 database "
                  "is being unpacked for the first use.")
            print("You should only see this message once.")
            os.system("bunzip2 " + database+".bz2")
            os.system("rm " + database+".bz2")
            print("Done.")
            # test and write note

        self.db_connection = sqlite3.connect(database)

        # SQLite optimizations from
        # https://phiresky.github.io/blog/2020/sqlite-performance-tuning/

        # conn.execute("pragma journal_mode = WAL;")
        self.db_connection.execute("pragma journal_mode = OFF;")
        self.db_connection.execute("pragma locking_mode = EXCLUSIVE;")
        self.db_connection.execute("pragma synchronous = OFF;")
        self.db_connection.execute("pragma temp_store = memory;")
        self.db_connection.execute("pragma mmap_size = 30000000000;")

        # setting db to read only
        self.db_connection.execute("pragma query_only = ON;")

        # getting db cursor
        self.db_cursor = self.db_connection.cursor()

    @staticmethod
    def process_sql_output(sql_output, seperator="\t", empty_char="NA"):

        processed_output = ""

        for entry in sql_output:

            to_print = [empty_char if i is None else str(i) for i in entry]
            processed_output += (seperator.join(to_print))+"\n"

        return processed_output

    def run_simple_select_query(self, output_field_list, query_data, input_field):

        # build SQL string from sanitized(!) field names
        sql_output_field_list = ",".join(output_field_list)

        if input_field == "Coordinates":
            # build SQL string
            coords = self.prepare_coordinates(self, query_data)

            sql = "SELECT " +sql_output_field_list+" FROM " +\
                  self.database_table_name + \
                  " WHERE Chr || ':' || " \
                  "Start || '|' || Stop in ({seq})".format(
                    seq=','.join(['?'] * len(coords)))
            sql_output = self.db_cursor.execute(sql, coords).fetchall()
        else:
            # build SQL string
            sql = "SELECT " + sql_output_field_list +\
                  " FROM " + self.database_table_name + \
                  " WHERE " + input_field + " in ({seq})".format(
                    seq=','.join(['?'] * len(query_data)))
            sql_output = self.db_cursor.execute(sql, query_data).fetchall()

        return sql_output

    def run_keyword_select_query(self, output_field_list,
                                 keyword_sql):

        # build SQL string from sanitized(!) field names
        sql_output_field_list = ",".join(output_field_list)

        # build SQL string
        sql = "SELECT " + sql_output_field_list +\
              " FROM " + self.database_table_name + \
              " WHERE " + keyword_sql + ";"

        sql_output = self.db_cursor.execute(sql).fetchall()

        return sql_output

    def run_circrna_query(self, circrna_id):

        # build SQL string
        sql_output = self.db_cursor.execute("SELECT "
                                            "* " +
                                            " FROM " + self.database_table_name +
                                            " INNER JOIN " + self.database_table_name + "_log " +
                                            " ON " + self.database_table_name + "_log.CircRNA_ID = " +
                                            self.database_table_name + ".CircRNA_ID " +
                                            " INNER JOIN " + self.database_table_name + "_db_info " +
                                            " ON " + self.database_table_name + "_db_info.DB_ID = " +
                                            self.database_table_name + "_log.DB_ID" +
                                            " WHERE " +
                                            " circhemy.CircRNA_ID == ? OR " +
                                            " CSNv1 == ? OR " +
                                            " circBase == ? OR " +
                                            " CircAtlas2 == ? OR " +
                                            " circRNADb == ? OR " +
                                            " circBank == ? OR " +
                                            " deepBase2 == ? OR " +
                                            " Circpedia2 == ? OR " +
                                            " riboCIRC == ? OR " +
                                            " exorBase2 == ? OR " +
                                            " ENSEMBL == ? OR " +
                                            " Gene == ? OR " +
                                            " Arraystar == ?;", (
                                                circrna_id, circrna_id, circrna_id,
                                                circrna_id, circrna_id, circrna_id,
                                                circrna_id, circrna_id, circrna_id,
                                                circrna_id, circrna_id, circrna_id, circrna_id)).fetchall()

        return sql_output

    def get_circrna_history_by_id(self, circrna_id):

        # build SQL string
        sql_output = self.db_cursor.execute("SELECT "
                                            "* " +
                                            " FROM " + self.database_table_name + "_log " +
                                            " INNER JOIN " + self.database_table_name + "_db_info " +
                                            " ON " + self.database_table_name + "_db_info.DB_ID = " +
                                            self.database_table_name + "_log.DB_ID" +
                                            " WHERE CircRNA_ID = ?", (circrna_id,)).fetchall()

        return sql_output