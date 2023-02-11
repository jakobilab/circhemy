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

import sqlite3


class Util(object):

    # global settings
    software_version = "0.0.2-dev"

    database_version = "DB-0.0.2-dev"

    program_name = "circhemy"

    support_email = program_name+"@jakobilab.org"

    support_web = "https://github.com/jakobilab/circconvert/ "+\
                  program_name +\
                  "/issues/new"

    news_url = "https://redmine.jakobilab.org/projects/circhemy/" \
               "news.atom?key=c616b9fb231445ac4ca65d94db1d3207382798e9"

    database_location = "../data/circconvert.sqlite3"

    database_table_name = "circrnadb"

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
        "circBase",
        "CircAtlas",
        "circRNADB",
        "circBank",
        "Deepbase2",
        "Circpedia2",
        "riboCIRC",
        "exorBase2",
        "Arraystar",
        "Pubmed"
    ]

    db_columns = ["Species"] + select_db_columns + ["Chr",
                                                    "Start",
                                                    "Stop",
                                                    "Strand",
                                                    "Genome"]

    external_db_urls = {"circBase": "http://www.circbase.org/cgi-bin/singlerecord.cgi?id=",
                        "CircAtlas": "http://159.226.67.237:8080/new/circ_detail.php?ID=",
                        "Circpedia2": "",
                        "circBank": "http://www.circbank.cn/infoCirc.html?id=",
                        "Deepbase2": "https://rna.sysu.edu.cn/deepbase3/subpages/ViewDetail_circRNA.php?spe=hg19&name=",
                        "Arraystar": "",
                        "circRNADB": "",
                        "riboCIRC": "http://www.ribocirc.com/rna_detail.php?dependent=Condition-independent&circ_id=",
                        "Circ2Disease": "http://bioinformatics.zju.edu.cn/Circ2Disease/browse_circRNA_result.php?circRNA=",
                        "Pubmed": "https://pubmed.ncbi.nlm.nih.gov/",
                        "exorBase2": "http://www.exorbase.org/exoRBaseV2/detail/detailInfo?kind=circRNA&id=",
                        "Chr": "",
                        "Start": "",
                        "Stop": "",
                        "Genome": "",
                        "Genome-Browser": "https://genome.ucsc.edu/cgi-bin/hgTracks?",
                        "Species": ""
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

    def check_input_return_found_circ_number(self, query_data, input_field):

        # build SQL string
        sql = "SELECT count(distinct("+input_field+")) FROM " + self.database_table_name + \
              " WHERE " + input_field + " in ({seq})".format(
                seq=','.join(['?'] * len(query_data)))

        sql_output = self.db_cursor.execute(sql, query_data).fetchall()

        # return ratio (0->1)
        return sql_output[0][0]/len(query_data),sql_output[0][0]

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

        # build SQL string
        sql = "SELECT " + sql_output_field_list +\
              " FROM " + self.database_table_name + \
              " WHERE " + input_field + " in ({seq})".format(
                seq=','.join(['?'] * len(query_data)))

        print(sql)
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
