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

# misc web functionality
from pathlib import Path
from typing import List, Any
from uuid import uuid4

# models for REST API endpoints
from pydantic import BaseModel, validator
from pygments.formatters import HtmlFormatter

# regex for circRNA parsing
import re
import os

# own util functions
import circhemy.common.util as common

from datetime import datetime


# core nicegui and web imports
from fastapi import Request, Response
from nicegui import Client, app, ui
from . import svg

# imports for parsing the Redmine news feed
import feedparser
from dateutil import tz

from dateutil.parser import parse

from io import StringIO
from html.parser import HTMLParser

# create util instance for the web app
util = common.Util

# add static files for fonts and favicon
app.add_static_files('/favicon', Path(__file__).parent / 'favicon')
app.add_static_files('/fonts', Path(__file__).parent / 'fonts')
app.add_static_files('/static', Path(__file__).parent / 'static')
app.add_static_files('/sitemap', 'web/sitemap/')

# set up global variables

# holds most variables for the convert form submission
ui_convert_form_values = dict()

# holds variables for query-based constraints that are added dynamically
ui_query_forms = list()

# setup SQLite connection
util.setup_database(util, util.database_location)

# initialize statistics chart on the righthand side
ui_convert_form_values['chart'], ui_convert_form_values['dbsize'],\
    ui_convert_form_values['chart2'] = util.database_stats(util)


def main():
    # run main application
    ui.run(title=util.program_name+" "+util.software_version,
           show=False,
           exclude="audio,chart colors,interactive_image,joystick,"
                   "keyboard,log,scene,video",
           binding_refresh_interval=0.1
           )

# util functions


class RedmineParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.reset()
        self.strict = False
        self.convert_charrefs= True
        self.text = StringIO()

    def handle_data(self, d):
        self.text.write(d)

    def get_data(self):
        return self.text.getvalue()


def redmine_remove_tags(html):
    s = RedmineParser()
    s.feed(html)
    return s.get_data()


def check_circrna_input_regex(strg, search=re.compile(r'[^A-Za-z0-9_\-|:\n\t]').search):
    return not bool(search(strg))


def check_text_field_input(upload_data) -> str:
    if not upload_data:
        circ_list = str(ui_convert_form_values['textfield'].value)
    else:
        circ_list = upload_data

    list_okay = check_circrna_input_regex(circ_list)

    db_selected = check_if_db_is_selected(ui_convert_form_values)
    ui_convert_form_values['circrna_found'].set_visibility(False)

    if circ_list and list_okay and db_selected:
        ui_convert_form_values['submit_button'].props(remove="disabled=true")
        ui_convert_form_values['submit_notification'].set_text(
            ui_update_found_circrnas(circ_list))

        return "Submit " + str(circ_list.count('\n')) + " circRNAs for ID conversion"

    elif not circ_list:
        ui_convert_form_values['submit_button'].props("disabled=true")
        return "Convert circRNA IDs"

    elif not list_okay:
        ui_convert_form_values['submit_button'].props("disabled=true")
        ui_convert_form_values['submit_notification']. \
            set_text("Allowed characters: A-Z,\\n,\\t,|,:,-,_ ")
        return "Unsupported characters detected in your circRNA list"

    elif not db_selected:
        ui_convert_form_values['submit_button'].props("disabled=true")

        return "No Input database selected"


def check_if_db_is_selected(form_values) -> List[Any]:
    checklist = []

    for item in form_values['db_checkboxes']:
        if item.value is True:
            checklist.append(item.text)

    return checklist


def check_query_text_field() -> None:

    if 'submit_query_button' in ui_convert_form_values:
        all_good = True

        for form in ui_query_forms:
            if not form['query'].value:
                all_good = False

        if all_good:
            ui_convert_form_values['submit_query_button'].props(remove="disabled=true")
        else:
            ui_convert_form_values['submit_query_button'].props("disabled=true")

        ui_convert_form_values['submit_query_button'].update()


def add_if_not_in_list(input_list=None, item_list=None):

    if item_list is None:
        item_list = []
    if input_list is None:
        input_list = []
    for item in item_list:
        if item not in input_list:
            input_list.append(item)

    return input_list


def ui_generate_readable_coordinate(value: int):

    return f'{value:,}'


def ui_generate_external_link(external_site: str, values: dict):
    # special case for genome browser links
    # not handled by normal external DB URL dict
    # since we need to build the link component first

    pos = values['Chr'] + ":" + \
          str(values['Start']) + \
          "-" + \
          str(values['Stop'])

    tmp_start = values['Start']
    tmp_stop = values['Stop']

    print_strand = ""
    if values['Strand']:
        print_strand = " ["+values['Strand']+"]"

    print_pos = values['Chr'] + ": " + f'{tmp_start:,}' + " - " + f'{tmp_stop:,}' + print_strand

    if external_site == "Coordinates" and values["Coordinates"]:
        # build pos format: chrXZY:1234-5789

        build = values['Genome']

        returnval = "<a style=\"text-decoration: underline;" \
                    "\" href=\"" + util.external_db_urls[
                        "Genome-Browser"] \
                    + "db=" + build + "&" \
                    + "position=" + pos + "\"" \
                    + " target=\"_blank\">" \
                    + print_pos + "</a>"

    elif external_site == "CSNv1" and values["CSNv1"]:

        pos = values['Chr'] + ":" + \
              str(values['Start']) + \
              "-" + \
              str(values['Stop'])

        build = values['Genome']

        returnval = "<a style=\"text-decoration: underline;\" href=\"/circrna/" +\
                    values['CSNv1'] + "\">" + values['CSNv1'] + \
                    "</a> <a style=\"text-decoration: underline;" \
                       "\" href=\"" + util.external_db_urls["CSNv1"] \
                    + "db=" + build + "&" \
                    + "position=" + pos + "&hgct_customText=https://static.jakobilab.org/circhemy/" +\
                    values['Genome'] + ".bed\"" \
                    + " target=\"_blank\">[view exons in Genome Browser]</a>"

    elif external_site == "Species" or external_site == "Gene":

        returnval = "<i>"+values[external_site].replace("_", " ")+"</i>"

    elif values[external_site] and util.external_db_urls[
                                external_site]:
        returnval = "<a style=\"text-decoration: underline;" \
                            "\" href=\"" + util.external_db_urls[
                                external_site] \
                            + str(values[external_site]) + "\" target=\"_blank\">" \
                            + str(values[external_site]) + "</a>"

    elif values[external_site] and not util.external_db_urls[
                                external_site]:
        returnval = str(values[external_site])
    else:
        returnval = ""

    return returnval


def ui_result_table_get_coordinates(input_list=None):

    if input_list is None:
        input_list = []
    return_dict = {"Chr": "", "Start": "", "Stop": "", "Genome": ""}

    for item in return_dict:
        index = 0
        for field in input_list:
            if field == item:
                return_dict[item] = index
            index = index + 1

    return return_dict


def ui_load_example_data() -> None:
    ui_convert_form_values['textfield'].value = \
                                        "chr1:100121447|100132793\n" \
                                        "chr1:100121550|100122379\n" \
                                        "chr1:100122380|100125942\n" \
                                        "chr1:100123210|100125942\n" \
                                        "chr1:1001237|1001530\n" \
                                        "chr1:100123872|100125942\n" \
                                        "chr1:100133404|100133583\n" \
                                        "chr1:100136882|100137085\n" \
                                        "chr1:100136882|100139711\n" \
                                        "chr1:100136882|100140514\n" \
                                     # "hsa-MYH9_0004\n" \
                                     # "hsa-MYH9_0005\n" \
                                     # "hsa-MYH10_0003\n" \
                                     # "hsa-MYH10_0016\n" \
                                     # "hsa-MYH10_0044\n" \
                                     # "hsa-MYH10_0075\n" \
                                     # "hsa-MYH10_0018\n" \
                                     # "hsa-MYH10_0002\n" \
                                     # "hsa-MYH14_0011\n" \
                                     # "hsa-MYH14_0013\n" \
                                     # "hsa-MYH14_0010\n" \
                                     # "hsa-MYH9_0116"


def ui_generate_result_table(input_id=None, output_ids=None, query_data=None):

    # initialize empty to allow for empty results
    output = ""

    output_fields = ""

    # function is called from REST API
    if input_id and type(input_id) != list:
        circrna_list = query_data

        output_fields = output_ids

        output = util.run_simple_select_query(util,
                                              output_ids,
                                              circrna_list,
                                              input_id
                                              )

    # REST API query gets input_id from type list
    # in this case the list holds the constraints for SQL query
    # construction
    elif input_id and type(input_id) is list:

        output_fields = output_ids

        sql_query = ""

        constraint_id = 0

        for constraint in input_id:

            if constraint_id > 0:
                sql_query += " " + constraint.operator1 + " "

            if constraint.operator2 == "LIKE":
                tmp = constraint.query.replace("*", "", 1)
                sql_query += constraint.field + " LIKE \"%" \
                             + tmp + "%\" "
            elif constraint.operator2 == "is":
                sql_query += constraint.field + " == \"" \
                             + constraint.query+ "\" "
            elif constraint.operator2 == ">":
                sql_query += constraint.field + " > \"" \
                             + constraint.query + "\" "
            elif constraint.operator2 == "<":
                sql_query += constraint.field + " < \"" \
                             + constraint.query + "\" "

            constraint_id = constraint_id + 1
        output = util.run_keyword_select_query(util,
                                               output_fields,
                                               sql_query)

    # function called from web convert module
    elif ui_convert_form_values['mode'] is "convert":

        output_fields = check_if_db_is_selected(ui_convert_form_values)

        # we always need these fields for genome browser links
        output_fields = add_if_not_in_list(output_fields, ["Chr",
                                                           "Start",
                                                           "Stop",
                                                           "Genome"
                                                           ])

        if "uploaded_data" in ui_convert_form_values:
            circrna_list = ui_convert_form_values['uploaded_data'].split('\n')
        else:
            circrna_list = ui_convert_form_values['textfield'].value.split('\n')

        output = util.run_simple_select_query(util,
                                              output_fields,
                                              circrna_list,
                                              ui_convert_form_values['db_checkbox'].value
                                              )

        if ui_convert_form_values['db_checkbox'].value not in output_fields \
                and ui_convert_form_values['db_checkbox'].value != "Coordinates":
            output_fields.insert(0, ui_convert_form_values['db_checkbox'].value)

    # function called from web query module
    elif ui_convert_form_values['mode'] is "query":

        output_fields = check_if_db_is_selected(ui_convert_form_values)

        # we always need these fields for genome browser links
        output_fields = add_if_not_in_list(output_fields, ["Chr",
                                                           "Start",
                                                           "Stop",
                                                           "Genome"
                                                           ])

        sql_query = ""

        for form in ui_query_forms:

            if 'operator1' in form:
                # this is an addon condition with two operators
                sql_query += " " + form['operator1'].value + " "

            if form['operator2'].value is "LIKE":
                tmp = form['query'].value.replace("*", "", 1)
                sql_query += form['field'].value + " LIKE \"%" \
                             + tmp + "%\" "
            elif form['operator2'].value is "is":
                sql_query += form['field'].value + " == \"" \
                             + form['query'].value + "\" "
            elif form['operator2'].value is ">":
                sql_query += form['field'].value + " > \"" \
                             + form['query'].value + "\" "
            elif form['operator2'].value is "<":
                sql_query += form['field'].value + " < \"" \
                             + form['query'].value + "\" "

        ui_query_forms.clear()

        output = util.run_keyword_select_query(util,
                                               output_fields,
                                               sql_query)

    full_list = list(output_fields)

    processed_output = ""

    # only add AG Grid structures if we are calling from web
    if not input_id:

        table_base_dict = {'defaultColDef': {
            'filter': True,
            'sortable': True,
            'resizable': True,
            'cellStyle': {'textAlign': 'left'},
            'headerClass': 'font-bold'
        },

            'columnDefs': [],
            'rowData': []
        }
    # REST API call, just return a more simple JSON-compatible table
    else:
        table_base_dict = {'columnDefs': [],
                           'rowData': []
                           }

    for item in full_list:
        table_base_dict['columnDefs'].append(
            {'headerName': item, 'field': item})

        # make sure we only call this if called from web,
        # otherwise this field is not initialized
        if not input_id:

            processed_output = processed_output + item \
                               + ui_convert_form_values['select2'].value

    # did we actually have SQL rows returned?
    if output:

        # holds index of coordinates + genome build in field list
        idx = ui_result_table_get_coordinates(full_list)

        for line in output:

            tmp_dict = dict()

            zipped = zip(line, full_list)

            # continue with integration of genome browser urls
            # we have fixed genome and potentially coordinates
            # we just have to tie all together in the special case of those
            # fields and manually get item information from the list
            # we probably have to add a counter to get to the correct position

            for item in zipped:

                if item[0] == "NA":
                    tmp_dict[item[1]] = ""

                # special case for genome browser links
                # not handled by normal external DB URL dict
                # since we need to build the link component first
                elif item[1] in ["Chr", "Start", "Stop"] and not input_id:

                    # build pos format: chrXZY:1234-5789
                    pos = line[idx['Chr']] + ":" + \
                          str(line[idx['Start']]) + \
                          "-" + \
                          str(line[idx['Stop']])

                    build = line[idx['Genome']]

                    tmp_dict[item[1]] = "<a style=\"text-decoration: underline;" \
                                        "\" href=\"" + util.external_db_urls[
                                            "Genome-Browser"] \
                                        + "db=" + build + "&" \
                                        + "pos=" + pos + "\"" \
                                        + " target=\"_blank\">" \
                                        + str(item[0]) + "</a>"

                elif item[1] == "CSNv1" and not input_id \
                        and line[idx['Genome']] in ["hg38", "mm10", "rn6"]:
                    # only if we use the new genome builds
                    # we re-use the genome browser integration but add the
                    # custom BED files

                    # build pos format: chrXZY:1234-5789
                    pos = line[idx['Chr']] + ":" + \
                          str(line[idx['Start']]) + \
                          "-" + \
                          str(line[idx['Stop']])

                    build = line[idx['Genome']]

                    tmp_dict[item[1]] = "<a style=\"text-decoration: underline;" \
                                        "\" href=\"" + util.external_db_urls[
                                            "Genome-Browser"] \
                                        + "db=" + build + "&" \
                                        + "pos=" + pos  \
                                        + "&hgct_customText=" +\
                                        util.external_db_urls["CSNv1"]\
                                        + build + ".bed\"" \
                                        + " target=\"_blank\">" \
                                        + str(item[0]) + "</a>"

                elif item[0] != "NA" \
                        and util.external_db_urls[item[1]] \
                        and not input_id:
                    tmp_dict[item[1]] = "<a style=\"text-decoration: underline;" \
                                        "\" href=\"" + util.external_db_urls[
                                            item[1]] \
                                        + str(item[0]) + "\" target=\"_blank\">" \
                                        + str(item[0]) + "</a>"
                else:
                    tmp_dict[item[1]] = item[0]

            table_base_dict['rowData'].append(tmp_dict)

        # remove last sep character
        processed_output = processed_output[:-1]

        # add new line for correct line break
        processed_output = processed_output + "\n"

        # form_values['table2'].style("height: 900px")

    # only set up web table if we are calling from web
    if not input_id:

        table = ui.aggrid(table_base_dict, html_columns=list(range(len(full_list))))

        table.style("width: 90%")
        table.style("text-align:center; "
                    "margin-left:auto; "
                    "margin-right:auto; ")
        table.update()

        processed_output = processed_output \
                           + util.process_sql_output(output,
                                                     seperator=
                                                     ui_convert_form_values[
                                                         'select2'].value,
                                                     empty_char=
                                                     ui_convert_form_values[
                                                         'select3'].value)

        # web return is the output and the nicegui table object
        return processed_output, table
    else:
        # return is the output and a simple table dictionary
        return processed_output, table_base_dict


def ui_update_found_circrnas(data) -> str:
    circrna_list = data.split('\n')
    circrna_list = list(filter(None, circrna_list))

    ratio, found = util.check_input_return_found_circ_number(util, input_field=
    ui_convert_form_values['db_checkbox'].value, query_data=circrna_list)

    if found > 0:
        ui_convert_form_values['circrna_found'].value = ratio
        ui_convert_form_values['circrna_found'].set_visibility(True)
    else:
        ui_convert_form_values['circrna_found'].value = ratio
        ui_convert_form_values['submit_button'].props("disabled=true")

    return str(found) + " of " + str(len(circrna_list)) + " CircRNA IDs found"


def ui_file_upload_handler(file) -> None:
    data = file.content.decode('UTF-8')
    check_text_field_input(data)
    ui_convert_form_values['uploaded_data'] = data


def ui_layout_add_left_drawer(convert=False) -> None:
    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4; '):

        ui_layout_generate_logo()

        ui.label('Select input ID type:').style("text-decoration: underline;")

        if convert:
            # Manually adding Coordinates here, as it's not a real DB field
            ui_convert_form_values['db_checkbox'] = ui.select(
                ["Coordinates"] + util.select_db_columns, value="Coordinates",
                label="ID format").style("width: 90%")

        with ui.column():
            ui.label('')
            ui.label('Select output fields:').style(
                "text-decoration: underline;")

            # ui.label("Genomic Coordinates")
            # with ui.row():
            #     tmp = [ui.checkbox('Chr', value=True),
            #            ui.checkbox('Start', value=True),
            #            ui.checkbox('Stop', value=True)]

            db_entries = []
            ui.label("Databases")

            with ui.row():
                with ui.column():
                    for item in range(int(len(util.select_db_columns) / 2)):
                        db_entries.append(
                            ui.checkbox(util.select_db_columns[item],
                                        value=util.active_db_columns[util.select_db_columns[item]]))

                with ui.column():
                    for item in range(int(len(util.select_db_columns) / 2),
                                      len(util.select_db_columns)):
                        db_entries.append(
                            ui.checkbox(util.select_db_columns[item],
                                        value=util.active_db_columns[util.select_db_columns[item]]))

            checkbox_list = db_entries
            ui_convert_form_values['db_checkboxes'] = checkbox_list
            ui.label('')

        with ui.column():
            ui.label('Select Output Format:').style(
                "text-decoration: underline;")

            ui_convert_form_values['select2'] = ui.select({"\t": "Tab-delimited [\\t]",
                                                ",": "Comma-delimited [,]",
                                                ";": "Semicolon-delimited [;]"},
                                                          value="\t",
                                                          label="Separator character") \
                .style("width: 90%")

            ui_convert_form_values['select3'] = ui.select({"NA": "NA",
                                                "\t": "Tab [\\t]",
                                                "": "Don't print anything"},
                                                          value="NA",
                                                          label="Placeholder"
                                                     " for unavailable "
                                                     "fields") \
                .style("width: 90%")


def ui_layout_add_head_html() -> None:
    ui.add_head_html(
        (Path(__file__).parent / 'static' / 'header.html').read_text())
    ui.add_head_html(
        f'<style>{HtmlFormatter(nobackground=True).get_style_defs(".codehilite")}</style>')
    ui.add_head_html(
        f"<style>{(Path(__file__).parent / 'static' / 'style.css').read_text()}</style>")


def ui_layout_generate_logo() -> None:
    # program name as link
    ui.html("<a href=\"/\">" + util.program_name + "</a>").style(
        'text-align: center; font-size: 26pt;')

    ui.image("/static/logo2.png"). \
        tooltip("al·che·my - noun - The medieval forerunner of chemistry, "
                "based on the supposed transformation of matter. "
                "\"A seemingly magical process of transformation, "
                "creation, or combination.\"")

    # subtitle
    ui.html("The alchemy of<br/>circular RNA ID conversion"). \
        tooltip('...are shown on mouse over').style(
        'text-align: center; font-size: 16pt;')

    ui.html("<br/>").style('text-align: center; font-size: 12t;')


def ui_layout_add_header() -> None:
    menu_items = {
        'What\'s new?': '/news',
        'CircRNA ID conversion': '/',
        'Browse databases': '/query',
        'CLI application': '/cli',
        'REST API access': '/rest',
        'About': '/about'
    }
    with ui.header() \
            .classes('items-center duration-400 p-0 px-4 no-wrap') \
            .style('box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1)'):

        with ui.row().classes('lg:hidden'):
            with ui.menu().classes('bg-primary text-white text-lg') as menu:
                for title, target in menu_items.items():
                    ui.menu_item(title,
                                 on_click=lambda _, target=target: ui.open(
                                     target))
            ui.button(on_click=menu.open).props('flat color=white icon=menu')
        with ui.row().classes('max-lg:hidden'):
            for title, target in menu_items.items():
                ui.link(title, target).classes(replace='text-lg text-white')
        with ui.link(target='https://github.com/jakobilab/'):
            svg.github().classes('fill-white scale-125 m-1')


def ui_layout_add_footer_and_right_drawer() -> None:
    with ui.right_drawer(fixed=False).style('background-color: #ebf1fa'):
        ui.label('Database Version ' + util.database_version + ' statistics')
        ui.label('')

        ui.label('Database by species/genome')

        chart = ui.chart(ui_convert_form_values['chart']).classes('w-full h-64') \
            .style("height: 350px")

        ui.label('Database by CircRNA ID')
        chart2 = ui.chart(ui_convert_form_values['chart2']).classes('w-full h-64') \
            .style("height: 350px")

    with ui.footer().style('background-color: #3874c8'):
        ui.label(util.program_name + " | software version " +
                 util.software_version + " | database version " +
                 util.database_version)
        ui.link('© 2023 Jakobi Lab', 'https://jakobilab.org')
        ui.link('Visit Jakobi Lab @ GitHub ', 'https://github.com/jakobilab/')


def ui_query_remove_conditions(container) -> None:
    if len(ui_query_forms) > 1:
        container.remove(-1)
        del ui_query_forms[-1]


def ui_query_add_conditions(container, new=False) -> None:
    query_values = dict()

    with container:
        with ui.row():
            if new:
                query_values['operator1'] = ui.select(["AND", "OR", "AND NOT"],
                                                      value="AND",
                                                      label="Query operator").style(
                    "width: 150px")

            query_values['field'] = ui.select(
                util.db_columns,
                value="circBase",
                label="Database field").style("width: 130px")

            query_values['operator2'] = ui.select(["is", "LIKE", ">", "<"],
                                                  value="is",
                                                  label="Query operator").style(
                "width: 130px")

            query_values['query'] = ui.input(label='Enter search term',
                                             placeholder='start typing',
                                             on_change=lambda e:
                                             check_query_text_field())

    ui_query_forms.append(query_values)


# application logic pages

# main landing page / also works as start page for converter function

@ui.page('/')
async def page_application_convert():
    ui_layout_add_head_html()
    ui_layout_add_header()

    with ui.card().style('width: 100%;') as card:
        ui.html("<strong>CircRNA ID conversion made easy</strong>")
        with ui.card_section():
            ui.markdown("""
Multiple databases gathered hundreds of thousands of circRNAs, however, most 
databases employ their own naming scheme, making it harder and harder to keep 
track of known circRNAs and their identifiers.
        """)

    ui_convert_form_values['mode'] = "convert"

    # ui.image('https://imgs.xkcd.com/comics/standards.png')

    ui_convert_form_values['textfield'] = ui.input(
        label='Please paste a list of circRNA IDs, one per line:',
        placeholder='start typing',
        on_change=lambda e: ui_convert_form_values['submit_button'].
        set_text(check_text_field_input(upload_data=None))). \
        props('type=textarea rows=30').style("width: 60%; ")

    ui_convert_form_values['or'] = ui.label('- OR -')

    ui_convert_form_values['upload'] = ui.upload(
        label="1) Click + to select file "
              "2) upload file via button to the right "
              "3) press 'convert circRNA IDs' button",
        on_upload=lambda e: ui_file_upload_handler(e.files[0])).style(
        "width: 60%")

    with ui.row():
        ui_convert_form_values['submit_button'] = \
            ui.button('Convert circRNA IDs', on_click=lambda:
            ui.open(page_application_display_results)).props("disabled=true")

        ui_convert_form_values['example_button'] = \
            ui.button('Load example data', on_click=lambda:
            ui_load_example_data())

    ui_convert_form_values['circrna_found'] = ui.linear_progress(show_value=False,
                                                                 value=0).style(
        "width: 60%; ")

    ui_convert_form_values['circrna_found'].set_visibility(False)
    ui_convert_form_values['submit_notification'] = ui.label('')

    ####################

    ui_layout_add_left_drawer(convert=True)

    ui_layout_add_footer_and_right_drawer()


# main landing page / also works as start page for converter function
@ui.page('/query')
async def page_application_query():
    ui_layout_add_head_html()
    ui_layout_add_header()

    ui_convert_form_values['mode'] = "query"

    ui_convert_form_values['chart'], ui_convert_form_values['dbsize'], ui_convert_form_values[
        'chart2'] = util.database_stats(util)

    ui_query_add_conditions(ui.column(), new=False)

    condition_row = ui.column()

    with ui.row():
        ui_convert_form_values['submit_query_button'] = \
            ui.button('Submit query', on_click=lambda:
            ui.open(page_application_display_results)).props("disabled=false")

        ui_convert_form_values['add_condition_button'] = \
            ui.button('Add condition', on_click=lambda:
            ui_query_add_conditions(condition_row, new=True))

        ui_convert_form_values['remove_condition_button'] = \
            ui.button('Remove condition', on_click=lambda:
            ui_query_remove_conditions(condition_row))

    ui_convert_form_values['circrna_found'] = ui.linear_progress(show_value=False, value=0).style("width: 60%; ")

    ui_convert_form_values['circrna_found'].set_visibility(False)

    ####################

    ui_layout_add_left_drawer()

    ui_layout_add_footer_and_right_drawer()


@ui.page('/results')
async def page_application_display_results():
    ui_layout_add_head_html()
    ui_layout_add_header()

    session_id = str(uuid4())

    # this just makes sure we built the landing page first
    if 'mode' in ui_convert_form_values:

        processed_output,ui_convert_form_values['table2'] = ui_generate_result_table()

        if not os.path.isdir("tmp/"):
            os.makedirs("tmp/")

        try:
            with open('tmp/' + session_id + ".csv", 'w') as f:
                f.write(processed_output)
        except FileNotFoundError:
            print("Output file could not be created")
            exit(-1)

        f.close()

        app.add_static_files('/download', 'tmp')

        with ui.row().classes('self-center'):
            ui.button('Download table',
                      on_click=lambda e: ui.open(
                          '/download/' + session_id + ".csv")) \
                .classes('self-center')

            ui.button('New query',
                      on_click=lambda e: ui.open('/')).classes('self-center')

    else:
        ui.open(page_application_convert)

        ui.html('<strong>Internal error encountered.</strong>'
                '<br/><a href=\"/\">Returning to main page</a>'
                ).style('text-align:center;')


@ui.page('/circrna/')
async def page_application_display_circrna_profile():
    ui_layout_add_head_html()
    ui_layout_add_header()

    with ui.left_drawer(top_corner=True, bottom_corner=False).style('background-color: #d7e3f4; '):
        ui_layout_generate_logo()

    # circRNA Name header
    ui.html("URL does not contain a known circRNA ID").style('font-size: 24pt;')

    ui.html('<br/><a href=\"/query/\">Click here to search the circhemy\'s circRNA database</a>')

    ui_layout_add_footer_and_right_drawer()


@ui.page('/circrna/{circ_id}', title=util.program_name_long)
async def page_application_display_circrna_profile(circ_id: str, client: Client):

    ui_layout_add_head_html()
    ui_layout_add_header()

    output = util.run_circrna_query(util, circrna_id=circ_id)

    if len(output) > 0:

        with ui.left_drawer(top_corner=True, bottom_corner=False).style('background-color: #d7e3f4; '):

            ui_layout_generate_logo()

            if len(output) > 1:
                ui.html("Multiple circRNAs matching "+"("+str(len(output))+")"+":").style('font-size: 14pt;')
                ui.html("&nbsp;").style('font-size: 14pt;')

            with ui.tabs() as tabs:
                tabs.style('background-color: #d7e3f4; height: 50%; border-radius: 0px;')
                tabs.props("vertical")
                tabs.props("dense")
                tabs.props("no-caps")
                tabs.props("inline-label")
                tabs.props("outside-arrows")
                tabs.props('align="left"')

                for hit in range(0, len(output)):
                    if output[hit][15]:
                        ui.tab(output[hit][16] + ":" + str(
                            output[hit][17]) + "-" + str(output[hit][18])+" ["+output[hit][20]+"]",
                               icon='change_circle').style("justify-content: initial; background-color: #d7e3f4; ")
                    else:
                        ui.tab(output[hit][16] + ":" + str(
                            output[hit][17]) + "-" + str(output[hit][18])+" ["+output[hit][20]+"]",
                               icon='donut_large').style("justify-content: initial; background-color: #d7e3f4; ")
                
                ui.html("&nbsp;").style('font-size: 14pt;')
                with ui.row():
                    ui.icon('change_circle').classes('text-2xl')
                    ui.label('CircRNA with new CSNv1 annotation')
                with ui.row():
                    ui.icon('donut_large').classes('text-2xl')
                    ui.html('CircRNA <u>without</u> CSNv1 annotation')

    # circRNA not found
    if len(output) == 0:

        with ui.left_drawer(top_corner=True, bottom_corner=False).style('background-color: #d7e3f4; '):

            ui_layout_generate_logo()

        # circRNA Name header
        ui.html("\""+circ_id+"\" is an unknown circRNA ID").style('font-size: 24pt;')
        ui.html('<br/><a href=\"/query/\">Click here to search the circhemy\'s circRNA database</a>')

    # at least one circRNA found
    if len(output) >= 1:

        # this is the name of the first tab, has to be known at time out panel setup
        first_tab = output[0][16] + ":" + str(output[0][17]) + "-" + str(output[0][18])+" ["+output[0][20]+"]"

        with ui.tab_panels(tabs, value=first_tab).style("background-color: #f8f8f8;"):

            # for each circRNA hit, e.g. different circRNAs or different genomes
            for hit in range(0, len(output)):

                output_dict = {}

                with ui.tab_panel(output[hit][16] + ":" + str(
                        output[hit][17]) + "-" + str(output[hit][18])+" ["+output[hit][20]+"]"):

                    local_columns = util.all_db_columns.copy()

                    local_columns.insert(0, "Stable circhemy database ID")

                    source_database_print = ""

                    for hit_id in range(0, len(local_columns)):
                        output_dict[local_columns[hit_id]] = output[hit][hit_id]
                        if str(output[hit][hit_id]) == circ_id:
                            source_database_print = local_columns[hit_id]

                    output_dict['Coordinates'] = 1
                    output_dict['Unspliced length'] = f"{output_dict['Stop'] - output_dict['Start']:,}" + " bp"

                    local_columns.append('Coordinates')
                    local_columns.append('Unspliced length')

                    with ui.row():
                        with ui.column():
                            with ui.column().style('width: 600px;'):
                                with ui.card().style('width: 600px;font-size: 12pt;'). \
                                        classes('column justify-between').style("background-color: #d7e3f4;") as card:
                                    with ui.row():
                                        with ui.column():
                                            ui.html(circ_id).style('font-size: 24pt;')
                                        with ui.column().classes("items-center"):
                                            ui.label(source_database_print)

                            if output_dict['CSNv1'] and source_database_print != 'CSNv1':
                                with ui.column().style('width: 600px;'):
                                    with ui.card().style('width: 600px;font-size: 12pt;'). \
                                            classes('column justify-between').\
                                            style("background-color: #d7e3f4;") as card:
                                        ui.html(output_dict['CSNv1']).style('font-size: 20pt;')

                            with ui.column().style('width: 600px;'):
                                with ui.card().style('width: 600px;font-size: 12pt;'). \
                                        classes('column justify-between'). \
                                        style("background-color: #d7e3f4;") as card:

                                    for hit_id in range(0, len(local_columns)):

                                        if local_columns[hit_id] not in ["Chr",
                                                                         "Start",
                                                                         "Stop",
                                                                         "Strand",
                                                                         "Stable circhemy database ID"]:
                                            with ui.row():
                                                with ui.column().style('width: 150px;'):
                                                    ui.html('<strong>' + local_columns[hit_id] + '</strong>')
                                                with ui.column():
                                                    ui.html(ui_generate_external_link(local_columns[hit_id], output_dict)).\
                                                       style('align:center')

                        with ui.column().style('width: 350px;'):
                            with ui.card().style('width: 350px;font-size: 12pt;'). \
                                    classes('column justify-between'). \
                                    style("background-color: #d7e3f4;") as card2:

                                output2 = util.get_circrna_history_by_id(util,
                                                                         circrna_id=output_dict['Stable circhemy database ID'])

                                with ui.row():
                                    with ui.column().style('width: 220px;'):
                                        ui.html("<strong>" + local_columns[0] + "<strong>")

                                    with ui.column().style('width: 30px;'):
                                        ui.html(str(output[hit][0]))

                                with ui.row():
                                    with ui.column().style('width: 80px;'):
                                        ui.html('<strong>DB Version</strong>')

                                    with ui.column().style('width: 80px;'):
                                        ui.html('<strong>Status</strong>')

                                    with ui.column().style('width: 80px;'):
                                        ui.html('<strong>Date</strong>')

                                for status_id in range(0, len(output2)):

                                    time_code = datetime.utcfromtimestamp(output2[status_id][5]).strftime('%m/%d/%Y')

                                    with ui.row():
                                        with ui.column().style('width: 80px;'):
                                            ui.html(str(output2[status_id][4]))

                                        with ui.column().style('width: 80px;'):
                                            ui.html(str(util.db_action_codes[output2[status_id][1]]))

                                        with ui.column().style('width: 80px;'):
                                            ui.html(str(time_code))

    ui_layout_add_footer_and_right_drawer()

    # update page title
    # https://github.com/zauberzeug/nicegui/discussions/830#discussioncomment-5714587
    await client.connected(timeout=5)
    await ui.run_javascript(f'document.title = "{"Circular RNA "+ circ_id + " - " + util.program_name_long or "nothing"}";',
                            respond=False)
# content pages


@ui.page('/about')
async def page_about():
    ui_layout_add_head_html()
    ui_layout_add_header()

    ui_layout_add_head_html()
    ui_layout_add_header()

    ui.html("About "+util.program_name).style('text-align: center; font-size: 26pt;')

    with ui.card().style('width: 100%;') as card:
        ui.html("Circular RNAs").style('text-align: center; font-size: 16pt;')
        with ui.card_section():
            ui.markdown("""
Circular RNAs (circRNAs) originate through back-splicing events from linear 
primary transcripts, are resistant to exonucleases, typically not 
polyadenylated, and have been shown to be highly specific for cell type and 
developmental stage. 

The prediction of circular RNAs is a multi-stage bioinformatics process starting 
with raw sequencing data and usually ending with a list of potential circRNA 
candidates which, depending on tissue and condition may contain hundreds to 
thousands of potential circRNAs. While there are a number of tools for the 
prediction process (e.g. circtools developed by our group) a unified naming 
convention for circRNA is not available.

Multiple databases gathered hundreds of thousands of circRNAs, however, most 
databases employ their own naming scheme, making it harder and harder to keep 
track of known circRNAs and their identifiers.
                """)

    with ui.card().style('width: 100%;') as card:
        ui.html("Circhemy").style('text-align: center; font-size: 16pt;')
        with ui.card_section():
            ui.markdown("""
We developed circhemy, a modular, Python3-based framework for circRNA ID 
conversion that unifies several functionalities in a single Python package. 
Three different routes are implemented within package to access more than 2 
million circRNA IDs:

Streamlined CLI application for direct access to the prepackaged local SQLite3 
database User-friendly Web frontend at circhemy.jakobilab.org A public REST API 
that enables direct access to the most recent ID database from HPC systems using 
curl or similar tools

Circhemy includes two different modes of action: convert and query. Convert 
allows the user to convert from one type of circRNA ID to a wide variety of 
other database identifiers, while query allows users to run direct queries on 
the circRNA database to extract circRNAs fulfilling a user-defined set of 
constraints.

Moreover, circhemy is the first circRNA resource that supports and integrates 
the first version of the **C**ircRNA **S**tandard **N**omenclature (abbreviated 
**CSNv1** in circhemy) as outlined in **["A guide to naming eukaryotic circular 
RNAs", Chen et al. 2023.](https://www.nature.com/articles/s41556-022-01066-9)**. 
Currently, circhemy contains CSNv1 names for nearly 1 million circRNAs of Human, 
mouse, and rat. 
            """)

    with ui.card().style('width: 100%;') as card:
        ui.html("Contact").style('text-align: center; font-size: 16pt;')
        with ui.card_section():
            ui.markdown("""
Circhemy is developed at the **[Jakobi Lab](https://jakobilab.org)**, part of 
the **[Translational Cardiovascular Research Center 
(TCRC)](https://phoenixmed.arizona.edu/tcrc)** in the Department of Internal 
Medicine at **[The University of Arizona College of Medicine
 – Phoenix](https://phoenixmed.arizona.edu/)**.

Contact: **<circhemy@jakobilab.org>**
    """)

    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4; '):
        ui_layout_generate_logo()

    ui_layout_add_footer_and_right_drawer()


@ui.page('/rest')
async def page_rest():
    ui_layout_add_head_html()
    ui_layout_add_header()

    ui.html("Representational State Transfer (REST) Interface").style(
        'text-align: center; font-size: 26pt;')

    with ui.card().style('width: 100%;') as card:
        ui.html("Introduction").style('text-align: center; font-size: 16pt;')
        with ui.card_section():
            ui.markdown("""
Representational State Transfer, or REST for short, allows users and software 
developers to easily access circhemy from within their own tools or pipelines. 
Circhemy's REST API uses JSON for input queries and returning output, making it 
easy to format queries from every programming language or even by hand. The REST
 API it publicly available and uses a fixed set of keywords to perform 
conversions or queries. Two examples for the two different modes of action are 
shown below.
    """)

    with ui.card().style('width: 100%;') as card:
        ui.html("Quick Start").style('text-align: center; font-size: 16pt;')
        with ui.card_section():
            ui.markdown("**Convert module**")
            ui.markdown("The convert module is able to convert from a range of "
                        "input circRNA ID into different one or more database"
                        "identifiers.")
            ui.markdown("""
Example: Convert a list of CircAtlas IDs into circBase and  
into CircPedia2 IDs, but also output CircAtlas IDs. 

                curl -X 'POST' 'https://circhemy.jakobilab.org/api/convert'
                  -H 'accept: application/json'
                  -H 'Content-Type: application/json'
                  -d '{
                      "input": "CircAtlas",
                      "output": ["CircPedia2","CircAtlas"],
                      "query": ["hsa-MYH9_0004","hsa-MYH9_0004"]
                      }'

""")
            ui.markdown("""

Output is returned as JSON-formatted string which can directly be used for AG 
Grid tables for any other postprocessing:

""")

            ui.markdown("""

                {
                  "columnDefs": [
                    {
                      "headerName": "circBase",
                      "field": "circBase"
                    },
                    {
                      "headerName": "Circpedia2",
                      "field": "Circpedia2"
                    }
                  ],
                  "rowData": [
                    {
                      "circBase": "hsa_circ_0004470",
                      "Circpedia2": "HSA_CIRCpedia_36582"
                    },
                    {
                      "circBase": "hsa_circ_0004470",
                      "Circpedia2": "HSA_CIRCpedia_36582"
                    }
                  ]
                }
""")

        with ui.card_section():
            ui.markdown("**Query module**")
            ui.markdown("The query module is able to retrieve circRNA IDs from "
                        "the internal database that fulfil a set of user-defined"
                        " constraints.")
            ui.markdown("""
Example: Retrieve all circRNAs with a CircAtlas ID containing *atf* or *xbp1*, 
return the IDs in circBase and circAtlas format:

                curl -X 'POST'
                  'https://circhemy.jakobilab.org/api/query'
                  -H 'accept: application/json'
                  -H 'Content-Type: application/json'
                  -d '{{
                  "input": [
                    {
                      "query": "pdia4",
                      "field": "CircAtlas",
                      "operator1": "AND",
                      "operator2": "LIKE"
                    },
                    {
                      "query": "rattus_norvegicus",
                      "field": "species",
                      "operator1": "AND",
                      "operator2": "IS"
                    }
                  ],
                  "output": [
                    "circBase",
                    "CircAtlas"
                  ]
                }}'


            """)
            ui.markdown("""

Output is returned as JSON-formatted string which can directly be used for AG 
Grid tables for any other postprocessing:

            """)

            ui.markdown("""

                {
                  "columnDefs": [
                    {
                      "headerName": "circBase",
                      "field": "circBase"
                    },
                    {
                      "headerName": "CircAtlas",
                      "field": "CircAtlas"
                    }
                  ],
                  "rowData": [
                    {
                      "circBase": "",
                      "CircAtlas": "hsa-NPPA_0001"
                    },
                    {
                      "circBase": "",
                      "CircAtlas": "hsa-NPPA_0002"
                    },
                    {
                      "circBase": "",
                      "CircAtlas": "hsa-NPPA-AS1_0001"
                    },
                    {
                      "circBase": "hsa_circ_0009871",
                      "CircAtlas": "hsa-NPPA-AS1_0004"
                    },
                    {
                      "circBase": "",
                      "CircAtlas": "hsa-NPPA-AS1_0002"
                    },
                    {
                      "circBase": "",
                      "CircAtlas": "hsa-NPPA-AS1_0003"
                    },
                    {
                      "circBase": "",
                      "CircAtlas": "hsa-NPPA_0001"
                    },
                    {
                      "circBase": "",
                      "CircAtlas": "hsa-NPPA_0002"
                    },
                    {
                      "circBase": "",
                      "CircAtlas": "hsa-NPPA-AS1_0001"
                    },
                    {
                      "circBase": "hsa_circ_0009871",
                      "CircAtlas": "hsa-NPPA-AS1_0004"
                    },
                    {
                      "circBase": "",
                      "CircAtlas": "hsa-NPPA-AS1_0002"
                    },
                    {
                      "circBase": "",
                      "CircAtlas": "hsa-NPPA-AS1_0003"
                    }
                  ]
                }
            """)

    with ui.card().style('width: 100%;') as card:
        ui.html("Documentation").style('text-align: center; font-size: 16pt;')
        with ui.card_section():
            ui.markdown("The full documentation is available on "
                        "***[GitHub](https://)***.")

    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4; '):
        ui_layout_generate_logo()

    ui_layout_add_footer_and_right_drawer()


@ui.page('/cli')
async def page_cli():
    ui_layout_add_head_html()
    ui_layout_add_header()

    ui.html("Command Line Interface (CLI)").style('text-align: center; font-size: 26pt;')

    with ui.card().style('width: 100%;') as card:
        ui.html("Setup").style('text-align: center; font-size: 16pt;')
        with ui.card_section():
            ui.markdown("""
The circhemy package is written in Python3 (>=3.7) and consists of two 
core modules, namely convert and query. The command line version requires 
only one external dependency, ``sqlite3``, for access to the internal SQLite3
database with circRNA ID data

Installation is managed through ``pip3 install circhemy`` or ``python3 setup.py
install`` when installed from the cloned GitHub repository. No sudo access is
required if the installation is executed with ``--user`` which will install the
package in a user-writeable folder. The binaries should be installed
to ``/home/$user/.local/bin/`` in case of Debian-based systems.

circhemy was developed and tested on Debian Buster, but should run with
any other distribution.

The installation requires running python on the command line:

    git clone https://github.com/dieterich-lab/circtools.git
    cd circtools
    python3 setup.py install --verbose --user
""")

    with ui.card().style('width: 100%;') as card:
        ui.html("Quick Start").style('text-align: center; font-size: 16pt;')
        with ui.card_section():
            ui.markdown("**Convert module**")
            ui.markdown("The convert module is able to convert from a range of "
                        "input circRNA ID into different one or more database"
                        "identifiers.")
            ui.markdown("""
Example: Convert a list of circatlas IDs read via STDIN from file `input.csv` 
into circpedia2 IDs, but also output circatlas IDs, while writing the output to
`/tmp/output.csv`:

``cat input.csv | circhemy convert  -q STDIN -i circatlas -o circpedia2 
circatlas -O /tmp/output.csv``""")

        with ui.card_section():
            ui.markdown("**Query module**")
            ui.markdown("The query module is able to retrieve circRNA IDs from "
                        "the internal database that fulfil a set of user-defined"
                        " constraints.")
            ui.markdown("""
Example: Retrieve a list of circbase and circatlas circRNA IDs that are located
on chromosome 3 of the species rattus norvegicus; only print out circRNAs from 
the rn6 genome build.

``circhemy query -o circbase circatlas -C chr3 
-s rattus_norvegicus -g rn6``
""")

    with ui.card().style('width: 100%;') as card:
        ui.html("Documentation").style('text-align: center; font-size: 16pt;')
        with ui.card_section():
            ui.markdown("The full documentation is available on "
                        "***[GitHub](https://)***.")

    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4; '):

        ui_layout_generate_logo()

    ui_layout_add_footer_and_right_drawer()


@ui.page('/news')
async def page_news():
    ui_layout_add_head_html()
    ui_layout_add_header()

    ui.html("<img style='text-align: center' src=\"static/news_small.png\">").style(
        'text-align: center; padding:10px;')

    # get feed
    d = feedparser.parse(util.news_url)

    max_posts = 20

    news_count = 1

    for post in d.entries:
        date_obj = parse(post.updated)

        with ui.card().style('width: 100%;') as card:
            ui.html("<strong>"+date_obj.astimezone(tz.tzlocal()).strftime(
            "%d %B, %Y")+"</strong>")
            ui.html(post.title).style('font-size: 1.25rem; line-height: 1.75rem;')
            with ui.card_section():
                ui.html(post.description)

        news_count += 1
        if news_count > max_posts:
            break

    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4; '):

        ui_layout_generate_logo()

    ui_layout_add_footer_and_right_drawer()


# error pages for error 404 and 500
@app.exception_handler(404)
async def exception_handler_404(request: Request, exception: Exception) -> Response:
    with Client(ui.page('')) as client:
        with ui.column().style('width: 100%; padding: 5rem 0; align-items: center; gap: 0'):
            ui.label('Error 404').style('font-size: 3.75rem; line-height: 1; padding: 1.25rem 0')
            ui.html("<strong>Page not found: </strong>" + str(request.url.path)).\
                style('font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
            ui.html(str(Request)).style('font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
            ui.html("<a href=\"/\"><img src=\"/static/error.png\"></a>").style('text-align: center; padding:10px;')
            ui.label('Please contact the server administrator at ' + util.support_email+' for additional support.').\
                style('font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
            ui.link('Click here to report an issue on GitHub.', util.support_web + "?title=circhemy 404 error" +
                "&body=URL: " + str(request.url) + "%0A%0APlease describe how the error occurred").\
                style('font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
    return client.build_response(request, 404)


@app.exception_handler(500)
async def exception_handler_500(request: Request, exception: Exception) -> Response:
    with Client(ui.page('')) as client:
        with ui.column().style('width: 100%; padding: 5rem 0; align-items: center; gap: 0'):
            ui.label('Internal server error').style('font-size: 3.75rem; line-height: 1; padding: 1.25rem 0')
            ui.html("<strong>Exception</strong>").style('font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
            ui.html(str(exception)).style('font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
            ui.html("<a href=\"/\"><img src=\"/static/error.png\"></a>").style('text-align: center; padding:10px;')
            ui.label('Please contact the server administrator at ' + util.support_email + ' for additional support.').\
                style('font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
            ui.link('Click here to report an issue on GitHub.', util.support_web + "?title=circhemy exception: " +
                str(exception) + "&body=URL: " + str(request.url) + "%0A%0APlease describe how the error occurred").\
                style('font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
    return client.build_response(request, 500)


# Below: REST API data definitions and endpoint functions for
# convert and query endpoint

# Data class for REST API calls from the convert module
class ConvertModel(BaseModel):
    input: str
    output: List[str]
    query: List[str]

    @validator('query', each_item=True)
    def circrna_id_pattern_check(cls, v):
        if not check_circrna_input_regex(v):
            raise ValueError('CircRNA IDs contains illegal characters.')
        if v == "":
            raise ValueError('No CircRNA IDs for conversion provided.')
        return v

    @validator('query')
    def circrna_id_check_length(cls, v):
        if len(v) == 0:
            raise ValueError('No CircRNA IDs for conversion provided.')
        return v

    @validator('input', allow_reuse=True)
    def database_name_check(cls, v):

        fields_allowed = list(util.external_db_urls.keys())

        if v not in fields_allowed:
            raise ValueError("Unsupported input field provided."
                             " Supported fields are: "
                             + ', '.join(fields_allowed)+
                             ". Field names are case-sensitive.")
        return v

    @validator('output', each_item=True)
    def database_name_check(cls, v):

        fields_allowed = list(util.external_db_urls.keys())

        if v not in fields_allowed:
            raise ValueError("Unsupported input field provided."
                             " Supported fields are: "
                             + ', '.join(fields_allowed)+
                             ". Field names are case-sensitive.")
        return v


# Data subclass for REST API calls from the query module
class ConstraintModel(BaseModel):
    query: str
    field: str
    operator1: str
    operator2: str

    @validator('query')
    def circrna_id_pattern_check(cls, v):
        if not check_circrna_input_regex(v):
            raise ValueError('Query contains illegal characters.')
        return v

    @validator('field')
    def field_name_check(cls, v):

        fields_allowed = list(util.external_db_urls.keys())

        if v not in fields_allowed:
            raise ValueError("Unsupported input field provided."
                             " Supported fields are: "
                             + ', '.join(fields_allowed)+
                             ". Field names are case-sensitive.")
        return v

    @validator('operator1')
    def operator1_check(cls, v):

        fields_allowed = ["AND", "OR", "AND NOT"]

        if v not in fields_allowed:
            raise ValueError("Unsupported operator1 provided."
                             " Supported fields are: "
                             + ', '.join(fields_allowed)+
                             ". Field names are case-sensitive.")
        return v

    @validator('operator2')
    def operator2_check(cls, v):

        fields_allowed = ["is", "LIKE", ">", "<"]

        if v not in fields_allowed:
            raise ValueError("Unsupported operator2 provided."
                             " Supported fields are: "
                             + ', '.join(fields_allowed)+
                             ". Field names are case-sensitive.")
        return v


# Data class for REST API calls from the query module
class QueryModel(BaseModel):
    input: List[ConstraintModel]
    output: List[str]


@app.post("/api/convert")
async def process_api_convert_call(data: ConvertModel):
    data, table = ui_generate_result_table(data.input, data.output, data.query)
    return table


@app.post("/api/query")
async def process_api_query_call(data: QueryModel):
    out, table = ui_generate_result_table(data.input, data.output)
    return table
