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

from pathlib import Path
from typing import List, Any, Union
from uuid import uuid4

from pydantic import BaseModel
from pygments.formatters import HtmlFormatter

import re

# own util functions
import common.util

# nicegui and web imports
from nicegui import ui
from nicegui import app
from web import svg

# create util instance for the web app
util = common.util.Util

# add static files for fonts and favicon
ui.add_static_files('/favicon', Path(__file__).parent / 'web' / 'favicon')
ui.add_static_files('/fonts', Path(__file__).parent / 'web' / 'fonts')

# set up global variables

# holds most variables for the convert form submission
form_values = dict()

# holds variables for query-based constraints that are added dynamically
query_forms = list()

# setup SQLite connection
util.setup_database(util, util.database_location)

# initialize statistics chart on the righthand side
form_values['chart'], form_values['dbsize'], form_values[
    'chart2'] = util.database_stats(util)

# run main application
ui.run(title=util.program_name, show=False)


def special_match(strg, search=re.compile(r'[^A-Za-z0-9_\-|:\n\t]').search):
    return not bool(search(strg))


def check_text_field_input(upload_data) -> str:
    if not upload_data:
        circ_list = str(form_values['textfield'].value)
    else:
        circ_list = upload_data
        print(upload_data)

    print(circ_list)

    list_okay = special_match(circ_list)

    db_selected = db_is_selected(form_values)
    form_values['circrna_found'].set_visibility(False)

    if circ_list and list_okay and db_selected:
        form_values['submit_button'].props(remove="disabled=true")
        form_values['submit_notification'].set_text(
            update_found_circrnas(circ_list))

        return "Submit " + \
            str(circ_list.count('\n') + 1) + \
            " circRNAs for ID conversion"

    elif not circ_list:
        form_values['submit_button'].props("disabled=true")
        return "Convert circRNA IDs"

    elif not list_okay:
        form_values['submit_button'].props("disabled=true")
        form_values['submit_notification']. \
            set_text("Allowed characters: A-Z,\\n,\\t,|,:,-,_ ")
        return "Unsupported characters detected in your circRNA list"

    elif not db_selected:
        form_values['submit_button'].props("disabled=true")

        return "No Input database selected"


def db_is_selected(form_values) -> List[Any]:
    checklist = []

    for item in form_values['db_checkboxes']:
        if item.value is True:
            checklist.append(item.text)

    return checklist


def load_example_data() -> None:
    form_values['textfield'].value = "hsa-MYH9_0004\n" \
                                     "hsa-MYH9_0005\n" \
                                     "hsa-MYH10_0003\n" \
                                     "hsa-MYH10_0016\n" \
                                     "hsa-MYH10_0044\n" \
                                     "hsa-MYH10_0075\n" \
                                     "hsa-MYH10_0018\n" \
                                     "hsa-MYH10_0002\n" \
                                     "hsa-MYH14_0011\n" \
                                     "hsa-MYH14_0013\n" \
                                     "hsa-MYH14_0010\n" \
                                     "hsa-MYH9_0116"


def generate_result_table(input_id=None, output_ids=None, query_data=None):

    # initialize empty to allow for empty results
    output = ""

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

            print(constraint.operator1)
            print(constraint.operator2)
            print(constraint.query)
            print(constraint.field)

            if constraint_id > 0:
                sql_query += " " + constraint.operator1 + " "

            if constraint.operator2 == "LIKE":
                tmp = constraint.query.replace("*", "", 1)
                sql_query += constraint.field + " LIKE \"%" \
                             + tmp + "%\" "
            elif constraint.operator2 is "is":
                sql_query += constraint.field + " == \"" \
                             + constraint.query+ "\" "
            elif constraint.operator2 is ">":
                sql_query += constraint.field + " > \"" \
                             + constraint.query + "\" "
            elif constraint.operator2 is "<":
                sql_query += constraint.field + " < \"" \
                             + constraint.query + "\" "

            constraint_id = constraint_id + 1

        print(sql_query)

        output = util.run_keyword_select_query(util,
                                               output_fields,
                                               sql_query)

    # function called from web convert module
    elif form_values['mode'] is "convert":

        output_fields = db_is_selected(form_values)

        if "uploaded_data" in form_values:
            circrna_list = form_values['uploaded_data'].split('\n')
        else:
            circrna_list = form_values['textfield'].value.split('\n')

        output = util.run_simple_select_query(util,
                                              output_fields,
                                              circrna_list,
                                              form_values['db_checkbox'].value
                                              )

        if form_values['db_checkbox'].value not in output_fields:
            output_fields.insert(0, form_values['db_checkbox'].value)

    # function called from web query module
    elif form_values['mode'] is "query":

        output_fields = db_is_selected(form_values)

        sql_query = ""
        print(sql_query+" after reset")

        for form in query_forms:
            print("form")

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

        query_forms.clear()

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

            processed_output = processed_output + item + form_values[
                'select2'].value

    # did we actually have SQL rows returned?
    if output:

        print("output")

        for line in output:

            tmp_dict = dict()

            for item in zip(line, full_list):

                if item[0] and util.external_db_urls[item[1]] and not input_id:
                    tmp_dict[item[1]] = "<a style=\"text-decoration: underline;" \
                                        "\" href="+util.external_db_urls[item[1]]\
                                        +item[0]+" target=\"_blank\">"\
                                        +item[0]+"</a>"
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

        table = ui.table(table_base_dict, html_columns=list(range(len(full_list))))

        table.style("width: 75%")
        table.style("text-align:center; "
                    "margin-left:auto; "
                    "margin-right:auto; ")
        table.update()

        processed_output = processed_output \
                           + util.process_sql_output(output,
                                                     seperator=
                                                     form_values[
                                                         'select2'].value,
                                                     empty_char=
                                                     form_values[
                                                         'select3'].value)

        # web return is the output and the nicegui table object
        return processed_output, table
    else:
        # return is the output and a simple table dictionary
        return processed_output, table_base_dict


def update_found_circrnas(data) -> str:
    circrna_list = data.split('\n')

    ratio, found = util.check_input_return_found_circ_number(util, input_field=
    form_values[
        'db_checkbox'].value, query_data=circrna_list)

    form_values['circrna_found'].value = ratio
    form_values['circrna_found'].set_visibility(True)

    return str(found) + " of " + str(len(circrna_list)) + " CircRNA IDs found"


def check_query_text_field() -> None:

    if 'submit_query_button' in form_values:
        all_good = True

        for form in query_forms:
            if not form['query'].value:
                all_good = False
        # print(form['query'].value)
        if all_good:
            form_values['submit_query_button'].props(remove="disabled=true")
        else:
            form_values['submit_query_button'].props("disabled=true")
        form_values['submit_query_button'].update()


def file_upload_handler(file) -> None:
    data = file.content.decode('UTF-8')
    check_text_field_input(data)
    form_values['uploaded_data'] = data


def add_left_drawer() -> None:
    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4; '):
        ui.image(
            'https://docs.circ.tools/en/latest/_static/circtools_150px.png')
        ui.label('Select input ID type:')

        form_values['db_checkbox'] = ui.select(
            ["circBase", "CircAtlas", "Circpedia2", "CircBank", "Deepbase2",
             "Arraystar", "CircRNADB"], value="CircAtlas",
            label="ID format").style("width: 90%")

        with ui.column():
            ui.label('')
            ui.label('Select output fields:')

            with ui.row():
                tmp = [ui.checkbox('Chr', value=True),
                ui.checkbox('Start', value=True),
                ui.checkbox('Stop', value=True)]

            checkbox_list = tmp + [
                ui.checkbox('circBase', value=True),
                ui.checkbox('CircAtlas', value=True),
                ui.checkbox('Circpedia2', value=True),
                ui.checkbox('CircBank', value=True),
                ui.checkbox('Deepbase2', value=True),
                ui.checkbox('Arraystar', value=True),
                ui.checkbox('CircRNADB', value=True)
            ]
            form_values['db_checkboxes'] = checkbox_list
            ui.label('')

        with ui.column():
            ui.label('Output formatting:')

            form_values['select2'] = ui.select({"\t": "Tab-delimited [\\t]",
                                                ",": "Comma-delimited [,]",
                                                ";": "Semicolon-delimited [;]"},
                                               value="\t",
                                               label="Separator character") \
                .style("width: 90%")

            form_values['select3'] = ui.select({"NA": "NA",
                                                "\t": "Tab [\\t]",
                                                "": "Don't print anything"},
                                               value="NA",
                                               label="Placeholder"
                                                     " for unavailable "
                                                     "fields") \
                .style("width: 90%")


def add_head_html() -> None:
    ui.add_head_html(
        (Path(__file__).parent / 'web' / 'static' / 'header.html').read_text())
    ui.add_head_html(
        f'<style>{HtmlFormatter(nobackground=True).get_style_defs(".codehilite")}</style>')
    ui.add_head_html(
        f"<style>{(Path(__file__).parent / 'web' / 'static' / 'style.css').read_text()}</style>")


def add_header() -> None:
    menu_items = {
        'Convert circRNA IDs': '/',
        'Search circRNA database': '/query',
        'CLI application': '/#cli',
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


def add_footer_and_right_drawer() -> None:
    with ui.right_drawer(fixed=False).style('background-color: #ebf1fa').props(
            'bordered'):
        ui.label('Database Version ' + util.database_version + ' statistics')
        ui.label('')

        ui.label('Database by species/genome')

        chart = ui.chart(form_values['chart']).classes('w-full h-64') \
            .style("height: 350px")

        ui.label('Database by CircRNA ID')
        chart2 = ui.chart(form_values['chart2']).classes('w-full h-64') \
            .style("height: 350px")

    with ui.footer().style('background-color: #3874c8'):
        ui.label(util.program_name + " | software version" +
                 util.software_version + " | database version" +
                 util.database_version)
        ui.link('| © 2022 Jakobi Lab |', 'https://jakobilab.org')
        ui.link('Visit Jakobi Lab @ GitHub ',
                'https://github.com/jakobilab/')


def remove_conditions(container) -> None:
    if len(query_forms) > 1:
        container.remove(-1)
        del query_forms[-1]


def add_conditions(container, new=False) -> None:
    query_values = dict()

    with container:
        with ui.row():
            if new:
                query_values['operator1'] = ui.select(["AND", "OR", "AND NOT"],
                                                      value="AND",
                                                      label="Query operator").style(
                    "width: 150px")

            query_values['field'] = ui.select(
                ["Chr","Start","Stop","circBase", "CircAtlas", "Circpedia2", "CircBank", "Deepbase2",
                 "Arraystar", "CircRNADB", "Chr", "Start", "Stop", "Genome",
                 "Species"],
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

    query_forms.append(query_values)


# main landing page / also works as start page for converter function
@ui.page('/')
async def landing_page():
    add_head_html()
    add_header()

    form_values['mode'] = "convert"

    # ui.image('https://imgs.xkcd.com/comics/standards.png')

    form_values['textfield'] = ui.input(
        label='Please paste a list of circRNA IDs, one per line:',
        placeholder='start typing',
        on_change=lambda e: form_values['submit_button'].
        set_text(check_text_field_input(upload_data=None))). \
        props('type=textarea rows=30').style("width: 60%; ")

    form_values['or'] = ui.label('- OR -')

    form_values['upload'] = ui.upload(
        file_picker_label="Click to select a file; upload via button to the right",
        on_upload=lambda e: file_upload_handler(e.files[0])).style(
        "width: 100%")

    with ui.row():
        form_values['submit_button'] = \
            ui.button('Convert circRNA IDs', on_click=lambda:
            ui.open(display_results_page)).props("disabled=true")

        form_values['example_button'] = \
            ui.button('Load example data', on_click=lambda:
            load_example_data())

    form_values['circrna_found'] = ui.linear_progress(show_value=False,
                                                      value=0).style(
        "width: 60%; ")

    form_values['circrna_found'].set_visibility(False)
    form_values['submit_notification'] = ui.label('')

    ####################

    add_left_drawer()

    add_footer_and_right_drawer()


# main landing page / also works as start page for converter function
@ui.page('/query')
async def query_page():
    add_head_html()
    add_header()

    form_values['mode'] = "query"

    form_values['chart'], form_values['dbsize'], form_values[
        'chart2'] = util.database_stats(util)

    add_conditions(ui.column(), new=False)

    condition_row = ui.column()

    with ui.row():
        form_values['submit_query_button'] = \
            ui.button('Submit query', on_click=lambda:
            ui.open(display_results_page)).props("disabled=false")

        form_values['add_condition_button'] = \
            ui.button('Add condition', on_click=lambda:
            add_conditions(condition_row, new=True))

        form_values['remove_condition_button'] = \
            ui.button('Remove condition', on_click=lambda:
            remove_conditions(condition_row))

    form_values['circrna_found'] = ui.linear_progress(show_value=False,
                                                      value=0).style(
        "width: 60%; ")

    form_values['circrna_found'].set_visibility(False)

    ####################

    add_left_drawer()

    add_footer_and_right_drawer()


@ui.page('/results')
async def display_results_page():
    add_head_html()
    add_header()

    session_id = str(uuid4())

    # this just makes sure we built the landing page first
    if 'mode' in form_values:

        processed_output,form_values['table2'] = generate_result_table()

        try:
            with open('tmp/' + session_id + ".csv", 'w') as f:
                f.write(processed_output)
        except FileNotFoundError:
            print("Output file could not be created")
            exit(-1)

        f.close()

        ui.add_static_files('/download', 'tmp')

        with ui.row().classes('self-center'):
            ui.button('Download table',
                      on_click=lambda e: ui.open(
                          '/download/' + session_id + ".csv")) \
                .classes('self-center')

            ui.button('New query',
                      on_click=lambda e: ui.open('/')).classes('self-center')

    else:
        ui.open(landing_page)

        ui.html('<strong>Internal error encountered.</strong>'
                '<br/><a href=\"/\">Returning to main page</a>'
                ).style('text-align:center;')


@ui.page('/about')
async def landing_page():
    add_head_html()
    add_header()

    ui.html('<strong>About page</strong>'
            '<br/><a href=\"/\">Returning to main page</a>'
            ).style('text-align:center;')

    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4; '):
        ui.image(
            'https://docs.circ.tools/en/latest/_static/circtools_150px.png')

    add_footer_and_right_drawer()


# Data class for REST API calls from the convert module
class ConvertModel(BaseModel):
    input: str
    output: List[str]
    query: List[str]


# Data subclass for REST API calls from the query module
class ConstraintModel(BaseModel):
    query: str
    field: str
    operator1: str
    operator2: str


# Data class for REST API calls from the query module
class QueryModel(BaseModel):
    input: List[ConstraintModel]
    output: List[str]


@app.post("/api/convert")
async def process_api_convert_call(data: ConvertModel):
    data, table = generate_result_table(data.input, data.output, data.query)
    return table


@app.post("/api/query")
async def process_api_query_call(data: QueryModel):
    print(data)
    out, table = generate_result_table(data.input, data.output)
    return table
