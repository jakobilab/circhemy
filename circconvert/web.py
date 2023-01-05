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
from typing import List, Any
from uuid import uuid4
from nicegui import ui
from web import svg
from pygments.formatters import HtmlFormatter
import re
import common.util

util = common.util.Util

ui.add_static_files('/favicon', Path(__file__).parent / 'web' / 'favicon')
ui.add_static_files('/fonts', Path(__file__).parent / 'web' / 'fonts')

form_values = dict()
query_forms = list()


def special_match(strg, search=re.compile(r'[^A-Za-z0-9_\-|:\n\t]').search):
    return not bool(search(strg))


def prepare_application() -> None:
    util.setup_database(util, util.database_location)


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


def submit_query_button_function():

    # work with input data
    output_fields = db_is_selected(form_values)

    sql_query = ""

    for form in query_forms:
        print("form")

        # this is the base condition, no additional operator1
        if 'operator1' not in form:

            print(form['query'].value)
            print(form['operator2'].value)
            print(form['field'].value)

        else:
            # this is an addon condition with two operators

            # since we already have a start query, we need to attach this one
            # via the operator1 supplied value, e.g. AND

            sql_query += " " + form['operator1'].value + " "

            print(form['query'].value)
            print(form['operator1'].value)
            print(form['operator2'].value)
            print(form['field'].value)

        if form['operator2'].value is "like":
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
        print(sql_query)

    # print(sql_query)



    # input_dict = dict(circbase=args.circbase_query,
    #                   circatlas=args.circatlas_query,
    #                   deepbase2=args.deepbase2_query,
    #                   circpedia2=args.circpedia2_query,
    #                   circbank=args.circbank_query,
    #                   arraystar=args.arraystar_query,
    #                   circrnadb=args.circrnadb_query,
    #                   species=args.species_query,
    #                   genome=args.genome_query,
    #                   chr=args.chr_query,
    #                   start=args.start_query,
    #                   stop=args.stop_query,
    #                   gene=args.gene_query
    #                   )
    #
    # # print(input_dict)
    #
    #
    #
    #
    #

    #
    # # remove last AND from query
    # sql_query = sql_query[:-4]
    #
    #
    #
    # output = util.run_simple_select_query(util,
    #                                       output_fields,
    #                                       circrna_list,
    #                                       form_values['db_checkbox'].value
    #                                       )
    #
    # output_fields.insert(0, form_values['db_checkbox'].value)
    # full_list = list((output_fields))
    #
    # processed_output = ""
    #
    # for item in full_list:
    #     form_values['table2'].options['columnDefs'].append(
    #         {'headerName': item, 'field': item})
    #
    #     processed_output = processed_output+item+form_values[
    #                                                'select2'].value
    #
    # for line in output:
    #
    #     tmp_dict = dict()
    #
    #     for item in zip(line, full_list):
    #         tmp_dict[item[1]] = item[0]
    #         # if item[0]:
    #         #     tmp_dict[item[1]] = "<strong>"+item[0]+"</strong>"
    #         # else:
    #         #     tmp_dict[item[1]] = item[0]
    #
    #
    #     form_values['table2'].options['rowData'].append(tmp_dict)
    #
    # # remove last sep character
    # processed_output = processed_output[:-1]
    #
    # # add new line for correct line break
    # processed_output = processed_output+"\n"
    #
    # # form_values['table2'].style("height: 900px")
    #
    # form_values['table2'].style("width: 60%")
    # form_values['table2'].style("text-align:center; "
    #                             "margin-left:auto; "
    #                             "margin-right:auto; ")
    # form_values['table2'].update()
    #
    # processed_output = processed_output + util.process_sql_output(output,
    #                                            seperator=form_values[
    #                                                'select2'].value,
    #                                            empty_char=form_values[
    #                                                'select3'].value)
    #
    # return processed_output


def submit_convert_button_function():
    # work with input data
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

    full_list = list(output_fields)

    processed_output = ""

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

    for item in full_list:
        table_base_dict['columnDefs'].append(
            {'headerName': item, 'field': item})

        processed_output = processed_output + item + form_values[
            'select2'].value

    for line in output:

        tmp_dict = dict()

        for item in zip(line, full_list):

            if item[0] and util.external_db_urls[item[1]]:
                tmp_dict[item[1]] = "<a style=\"text-decoration: underline;\" href="+util.external_db_urls[item[1]]+item[0]+" target=\"_blank\">"+item[0]+"</a>"
            else:
                tmp_dict[item[1]] = item[0]

        table_base_dict['rowData'].append(tmp_dict)

    # remove last sep character
    processed_output = processed_output[:-1]

    # add new line for correct line break
    processed_output = processed_output + "\n"

    # form_values['table2'].style("height: 900px")

    table = ui.table(table_base_dict, html_columns=list(range(len(full_list))))

    table.style("width: 75%")
    table.style("text-align:center; "
                                "margin-left:auto; "
                                "margin-right:auto; ")
    table.update()

    processed_output = processed_output + util.process_sql_output(output,
                                                                  seperator=
                                                                  form_values[
                                                                      'select2'].value,
                                                                  empty_char=
                                                                  form_values[
                                                                      'select3'].value)

    return processed_output, table


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
        'About': '/#about'
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
                ["circBase", "CircAtlas", "Circpedia2", "CircBank", "Deepbase2",
                 "Arraystar", "CircRNADB", "Chr", "Start", "Stop", "Genome",
                 "Species"],
                value="circBase",
                label="Database field").style("width: 130px")

            query_values['operator2'] = ui.select(["is", "like", ">", "<"],
                                                  value="is",
                                                  label="Query operator").style(
                "width: 130px")

            query_values['query'] = ui.input(label='Enter search term',
                                             placeholder='start typing', on_change=lambda e: check_query_text_field())

    query_forms.append(query_values)


# main landing page / also works as start page for converter function
@ui.page('/')
async def landing_page():
    add_head_html()
    add_header()

    # ui.image('https://imgs.xkcd.com/comics/standards.png')

    prepare_application()

    form_values['chart'], form_values['dbsize'], form_values[
        'chart2'] = util.database_stats(util)

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

            checkbox_list = [
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

    add_drawers()


# main landing page / also works as start page for converter function
@ui.page('/query')
async def query_page():
    add_head_html()
    add_header()

    prepare_application()

    form_values['chart'], form_values['dbsize'], form_values[
        'chart2'] = util.database_stats(util)

    add_conditions(ui.column(), new=False)

    condition_row = ui.column()

    with ui.row():
        form_values['submit_query_button'] = \
            ui.button('Submit query', on_click=lambda:
            ui.open(display_query_results_page)).props("disabled=false")

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

    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4; '):
        ui.image(
            'https://docs.circ.tools/en/latest/_static/circtools_150px.png')

        with ui.column():
            ui.label('')
            ui.label('Select output fields:')

            checkbox_list = [
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

    add_drawers()


@ui.page('/convert_results')
async def display_results_page():
    add_head_html()
    add_header()

    session_id = str(uuid4())

    # this just makes sure we built the landing page first
    if 'or' in form_values:
        # form_values['table2'] = ui.table({'defaultColDef': {
        #     'filter': True,
        #     'sortable': True,
        #     'resizable': True,
        #     'cellStyle': {'textAlign': 'left'},
        #     'headerClass': 'font-bold'
        # },
        #
        #     'columnDefs': [],
        #     'rowData': []
        # }, html_columns=[0])

        processed_output,form_values['table2'] = submit_convert_button_function()

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


@ui.page('/results')
async def display_query_results_page():
    add_head_html()
    add_header()

    session_id = str(uuid4())

    submit_query_button_function()

    # # this just makes sure we built the landing page first
    # if 'or' in form_values:
    #     form_values['table2'] = ui.table({'defaultColDef': {
    #         'filter': True,
    #         'sortable': True,
    #         'resizable': True,
    #         'cellStyle': {'textAlign': 'left'},
    #         'headerClass': 'font-bold'
    #     },
    #
    #         'columnDefs': [],
    #         'rowData': []
    #     }, html_columns=[0])
    #
    #     processed_output = smubit_convert_button_function()
    #
    #     try:
    #         with open('tmp/' + session_id + ".csv", 'w') as f:
    #             f.write(processed_output)
    #     except FileNotFoundError:
    #         print("Output file could not be created")
    #         exit(-1)
    #
    #     f.close()
    #
    #     ui.add_static_files('/download', 'tmp')
    #
    #     with ui.row().classes('self-center'):
    #         ui.button('Download table',
    #                   on_click=lambda e: ui.open(
    #                       '/download/' + session_id + ".csv")) \
    #             .classes('self-center')
    #
    #         ui.button('New query',
    #                   on_click=lambda e: ui.open('/')).classes('self-center')
    #
    # else:
    #     ui.open(landing_page)
    #
    #     ui.html('<strong>Internal error encountered.</strong>'
    #             '<br/><a href=\"/\">Returning to main page</a>'
    #             ).style('text-align:center;')


def add_drawers() -> None:
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


# run main application
ui.run(title=util.program_name, show=False)
