#!/usr/bin/env python3
# Copyright (C) 2024 Tobias Jakobi
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
from datetime import datetime

from pydantic import BaseModel, field_validator
from pygments.formatters.html import HtmlFormatter

import re

import circhemy.common.util as common

from fastapi import Request, Response
from nicegui import Client, app, ui
from . import svg

from html import escape
from urllib.parse import quote

# ========================================================
# MODULE-LEVEL READ-ONLY CONFIGURATION
# These are safe to share across all users because they
# are never mutated after initialization.
# ========================================================

util = common.Util

download_dir = Path("tmp")
download_dir.mkdir(exist_ok=True)

# Static file mounts
app.add_static_files('/favicon', Path(__file__).parent / 'favicon')
app.add_static_files('/fonts', Path(__file__).parent / 'fonts')
app.add_static_files('/static', Path(__file__).parent / 'static')
app.add_static_files('/sitemap', 'web/sitemap/')
app.add_static_files('/download', download_dir)

# Setup SQLite connection (read-only after init)
util.setup_database(util, util.database_location)

# Pre-compute database statistics once at startup (read-only)
_chart_data, _dbsize_data, _chart2_data = util.database_stats(util)


def main():
    ui.run(
        title=util.program_name_long + " - Release " + util.database_version,
        show=False,
        binding_refresh_interval=0.1,
        storage_secret="circhemy-secret-change-in-production",
    )


# ========================================================
# PURE UTILITY FUNCTIONS (stateless, safe for all users)
# ========================================================

def html_link(url: str, text: str, external: bool = True) -> str:
    """Generate an HTML link string."""
    target = ' target="_blank" rel="noopener noreferrer"' if external else ''
    return (
        f'<a style="text-decoration: underline;" '
        f'href="{escape(url, quote=True)}"{target}>'
        f'{escape(text)}'
        f'</a>'
    )


def check_circrna_input_regex(strg, search=re.compile(
    r'[^A-Za-z0-9_\-|:\n\t]').search):
    """Validate that input contains only allowed characters."""
    return not bool(search(strg))


def add_if_not_in_list(input_list=None, item_list=None):
    """Add items to list if not already present, maintaining order."""
    if item_list is None:
        item_list = []
    if input_list is None:
        input_list = []
    for item in item_list[::-1]:
        if item not in input_list:
            input_list.insert(0, item)
    return input_list


def ui_generate_readable_coordinate(value: int):
    return f'{value:,}'


def ui_result_table_get_coordinates(input_list=None):
    """Get indices of coordinate columns in the field list."""
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


def ui_generate_external_link(external_site: str, values: dict):
    """Generate external link HTML for a given database field."""
    pos = values['Chr'] + ":" + \
          str(values['Start']) + "-" + \
          str(values['Stop'])

    print_strand = ""
    if values['Strand']:
        print_strand = " [" + values['Strand'] + "]"

    print_pos = values['Chr'] + ": " + \
                f'{values["Start"]:,}' + " - " + \
                f'{values["Stop"]:,}' + print_strand

    if external_site == "Coordinates" and values["Coordinates"]:
        return html_link(
            util.external_db_urls["Genome-Browser"] +
            "?db=" + values["Genome"] +
            "&position=" + quote(pos, safe=''),
            print_pos)
    elif external_site == "CSNv1" and values["CSNv1"]:
        pos = values['Chr'] + ":" + \
              str(values['Start']) + "-" + \
              str(values['Stop'])
        build = values['Genome']
        return (f"{html_link('/', values['CSNv1'])} "
                f"<a style='text-decoration: underline;' "
                f"href='{util.external_db_urls['CSNv1']}db={build}"
                f"&position={pos}"
                f"&hgct_customText={util.external_db_urls['CSNv1']}{build}.bed' "
                f"target='_blank'>[view exons in Genome Browser]</a>")
    elif external_site == "Species" or external_site == "Gene" and values[external_site]:
        return f"<i>{values[external_site].replace('_', ' ')}</i>"
    elif values[external_site] and util.external_db_urls[external_site]:
        return html_link(
            util.external_db_urls[external_site] + quote(str(values[external_site]), safe=''),
            str(values[external_site])
        )
    elif values[external_site] and not util.external_db_urls[external_site]:
        return escape(str(values[external_site]))
    else:
        return ""


def ui_layout_circrna_header_meta(input_dict, circrna_id):
    """Generate meta description for a circRNA profile page."""
    if 'CSNv1' in input_dict and input_dict['CSNv1'] is not None:
        csnv1 = ", " + input_dict['CSNv1']
    else:
        csnv1 = ""

    pos = (input_dict['Chr'] + ":" +
           str(input_dict['Start']) + "-" + str(input_dict['Stop']))
    species = input_dict['Species']
    build = input_dict['Genome']

    if 'Gene' in input_dict and input_dict['Gene'] is not None:
        gene = ", originating from " + input_dict['Gene']
    else:
        gene = ""

    meta_string = ("CircRNA " + circrna_id + csnv1 +
                   ", " + pos + " (" + build + ")" + gene +
                   " in circhemy, the circRNA ID database")
    return meta_string


# ========================================================
# SHARED UI LAYOUT FUNCTIONS (create UI elements per-call,
# read only from module-level constants)
# ========================================================

def ui_layout_add_head_html(input_dict=None, circrna_id=None) -> None:
    """Add HTML head elements including meta tags and styles."""
    if input_dict:
        ui.add_head_html(
            "<meta name=\"description\" content=\" " +
            ui_layout_circrna_header_meta(input_dict, circrna_id) +
            " \"/>")
    else:
        ui.add_head_html(
            "<meta name=\"description\" content=\"Circhemy is an "
            "easy-to-use web service and command line tool helping"
            " researchers to convert between different"
            " circular RNA ID formats.\"/>")

    ui.add_head_html(
        (Path(__file__).parent / 'static' / 'header.html').read_text())
    ui.add_head_html(
        f'<style>{HtmlFormatter(nobackground=True).get_style_defs(".codehilite")}</style>')
    ui.add_head_html(
        f"<style>{(Path(__file__).parent / 'static' / 'style.css').read_text()}</style>")


def ui_layout_generate_logo() -> None:
    """Generate the circhemy logo and subtitle."""
    ui.image('/static/logo_small_2.png').classes('w-[320px]').tooltip(
        "al·che·my - noun - The medieval forerunner of chemistry, "
        "based on the supposed transformation of matter. "
        "\"A seemingly magical process of transformation, "
        "creation, or combination.\"")
    ui.html("<i>circhemy</i> - the alchemy of<br/>circular RNA ID conversion"). \
        tooltip('...are shown on mouse over').style(
        'text-align: center; font-size: 16pt;').style("width: 320px").classes(
        'q-py-none')


def ui_layout_add_header() -> None:
    """Generate the top navigation header bar."""
    menu_items = {
        'What\'s new?': '/news',
        'CircRNA ID conversion': '/',
        'Browse databases': '/query',
        'CLI application': '/cli',
        'REST API access': '/rest',
        'About': '/about'
    }

    with ui.header().style("height: 40px").classes('q-pa-xs') \
            .classes('items-center duration-400 p-0 px-4 no-wrap') \
            .style('box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1)'):
        with ui.row().classes('lg:hidden'):
            with ui.menu().classes('bg-primary text-white text-lg') as menu:
                for title, target in menu_items.items():
                    ui.menu_item(title,
                                 on_click=lambda _, target=target:
                                 ui.navigate.to(target))
            ui.button(on_click=menu.open).props('flat color=white icon=menu')
        with ui.row().classes('max-lg:hidden'):
            for title, target in menu_items.items():
                ui.link(title, target).classes(replace='text-lg text-white')
        with ui.link(target='https://github.com/jakobilab/'):
            svg.github().classes('fill-white scale-125 m-1')


def ui_layout_add_footer_and_right_drawer() -> None:
    """Add the right drawer with stats and the footer.
    Uses module-level read-only chart data — safe for all users."""
    with ui.right_drawer(fixed=False).style('background-color: #ebf1fa'):
        ui.label('Database Version ' + util.database_version + ' statistics')
        ui.label('')
        ui.label('Database by species/genome')
        ui.highchart(_chart_data).classes('w-full h-64').style("height: 330px")
        ui.label('Database by CircRNA ID')
        ui.highchart(_chart2_data).classes('w-full h-64').style("height: 330px")

    with ui.footer().style('background-color: #3874c8').style(
            "height: 30px").classes('q-pa-xs'):
        ui.label(" circhemy version " +
                 util.software_version + " | database version " +
                 util.database_version)
        ui.link('| Jakobi Lab @ GitHub ', 'https://github.com/jakobilab/')
        ui.link(' | © 2024 Jakobi Lab', 'https://jakobilab.org')


def build_result_table_for_web(output_fields, output, mode,
                               separator='\t', empty_char='NA'):
    """Build an AG Grid table for the web UI.

    This is a pure function — all data is passed explicitly.
    No global state is read or written.

    Args:
        output_fields: List of column names to display.
        output: SQL query result rows.
        mode: 'convert' or 'query'.
        separator: Character to separate fields in download output.
        empty_char: Placeholder for missing values in download output.

    Returns:
        Tuple of (processed_output_string, aggrid_table_element).
    """
    full_list = list(output_fields)
    processed_output = ""

    table_base_dict = {
        'defaultColDef': {
            'filter': True,
            'sortable': True,
            'resizable': True,
            'cellStyle': {'textAlign': 'left'},
            'headerClass': 'font-bold',
        },
        'columnDefs': [],
        'rowData': []
    }

    for item in full_list:
        table_base_dict['columnDefs'].append(
            {'headerName': item, 'field': item})
        if mode == "convert":
            processed_output = processed_output + item + separator

    if output:
        idx = ui_result_table_get_coordinates(full_list)

        for line in output:
            tmp_dict = dict()
            zipped = zip(line, full_list)

            for item in zipped:
                if item[0] == "NA":
                    tmp_dict[item[1]] = ""
                elif item[1] in ["Chr", "Start", "Stop"]:
                    pos = (str(line[idx['Chr']]) + ":" +
                           str(line[idx['Start']]) + "-" +
                           str(line[idx['Stop']]))
                    build = line[idx['Genome']]
                    tmp_dict[item[1]] = (
                        f'<a style="text-decoration: underline;" '
                        f'href="{util.external_db_urls["Genome-Browser"]}'
                        f'db={build}&pos={pos}" '
                        f'target="_blank">{item[0]}</a>')
                elif item[1] == "CSNv1" and line[idx['Genome']] in ["hg38", "mm10", "rn6"]:
                    pos = (str(line[idx['Chr']]) + ":" +
                           str(line[idx['Start']]) + "-" +
                           str(line[idx['Stop']]))
                    build = line[idx['Genome']]
                    tmp_dict[item[1]] = (
                        f'<a style="text-decoration: underline;" '
                        f'href="{util.external_db_urls["Genome-Browser"]}'
                        f'db={build}&pos={pos}'
                        f'&hgct_customText={util.external_db_urls["CSNv1"]}'
                        f'{build}.bed" target="_blank">{item[0]}</a>')
                elif item[1] == "CircRNA_ID":
                    tmp_dict[item[1]] = (
                        f'<a style="text-decoration: underline;" '
                        f'href="{util.external_db_urls["CircRNA_ID"]}'
                        f'{item[0]}" target="_blank">Circhemy profile</a>')
                elif (item[0] is not None and
                      util.external_db_urls.get(item[1])):
                    tmp_dict[item[1]] = html_link(
                        util.external_db_urls[item[1]] + quote(str(item[0]), safe=''),
                        str(item[0]))
                else:
                    tmp_dict[item[1]] = item[0]

            table_base_dict['rowData'].append(tmp_dict)

        # Remove last separator character and add newline
        processed_output = processed_output[:-1] + "\n"

    # Create the AG Grid UI element with html_columns support [8]
    table = ui.aggrid(
        table_base_dict,
        html_columns=list(range(len(full_list)))
    )
    table.style(add='height: calc(100vh - 220px)')
    table.props("domLayout='autoHeight'")
    table.on('firstDataRendered',
             lambda: table.run_grid_method('autoSizeAllColumns'))
    table.update()

    # Build downloadable text output
    if mode == "convert":
        processed_output = processed_output + util.process_sql_output(
            output, seperator=separator, empty_char=empty_char)
    else:
        processed_output = processed_output + util.process_sql_output(output)

    return processed_output, table


def build_result_table_for_api(input_id, output_ids, query_data=None):
    """Build a result table for REST API calls.

    Pure function — no global state accessed.

    Args:
        input_id: Input field name (str) or list of ConstraintModel for query.
        output_ids: List of output field names.
        query_data: List of circRNA IDs (for convert mode).

    Returns:
        Tuple of (processed_output_string, table_dict).
    """
    output_fields = output_ids

    if not isinstance(input_id, list):
        # Convert mode API call
        output = util.run_simple_select_query(
            util, output_ids, query_data, input_id)
    else:
        # Query mode API call
        sql_query = ""
        constraint_id = 0
        for constraint in input_id:
            if constraint_id > 0:
                sql_query += " " + constraint.operator1 + " "
            if constraint.operator2 == "LIKE":
                tmp = constraint.query.replace("*", "", 1)
                sql_query += constraint.field + ' LIKE "%' + tmp + '%" '
            elif constraint.operator2 == "IS":
                sql_query += constraint.field + ' == "' + constraint.query + '" '
            elif constraint.operator2 == ">":
                sql_query += constraint.field + ' > "' + constraint.query + '" '
            elif constraint.operator2 == "<":
                sql_query += constraint.field + ' < "' + constraint.query + '" '
            constraint_id += 1

        output = util.run_keyword_select_query(util, output_fields, sql_query)

    full_list = list(output_fields)
    processed_output = ""

    table_base_dict = {'columnDefs': [], 'rowData': []}

    for item in full_list:
        table_base_dict['columnDefs'].append(
            {'headerName': item, 'field': item})

    if output:
        for line in output:
            tmp_dict = dict()
            zipped = zip(line, full_list)
            for item in zipped:
                tmp_dict[item[1]] = item[0]
            table_base_dict['rowData'].append(tmp_dict)

    return processed_output, table_base_dict


# ========================================================
# APPLICATION PAGES — ALL STATE IS LOCAL PER-VISIT
# ========================================================

@ui.page('/')
async def page_application_convert():
    """Main landing page / circRNA ID conversion form.

    All mutable state is local to this function call.
    Each visitor gets their own isolated state."""

    # Must await connection BEFORE accessing app.storage.tab
    await ui.context.client.connected()

    ui_layout_add_head_html()
    ui_layout_add_header()

    # ---- Per-visit local state (replaces global ui_convert_form_values) ----
    page_state = {
        'mode': 'convert',
        'textfield': None,
        'db_checkbox': None,
        'db_checkboxes': [],
        'select2': None,
        'select3': None,
        'submit_button': None,
        'submit_notification': None,
        'circrna_found': None,
        'uploaded_data': None,
    }

    # ---- Closures that capture page_state ----

    def check_if_db_is_selected_local() -> List[Any]:
        """Return list of selected output database checkbox labels."""
        checklist = []
        for item in page_state['db_checkboxes']:
            if item.value is True:
                checklist.append(item.text)
        return checklist

    def ui_update_found_circrnas_local(data: str) -> str:
        """Update progress indicator with found circRNA count."""
        circrna_list = data.split('\n')
        circrna_list = list(filter(None, circrna_list))
        ratio, found = util.check_input_return_found_circ_number(
            util,
            input_field=page_state['db_checkbox'].value,
            query_data=circrna_list
        )
        if found > 0:
            page_state['circrna_found'].value = ratio
            page_state['circrna_found'].set_visibility(True)
        else:
            page_state['circrna_found'].value = ratio
            page_state['submit_button'].props("disabled=true")
        return f"{found} of {len(circrna_list)} CircRNA IDs found"

    def check_text_field_input_local(upload_data=None) -> str:
        """Validate input and update UI elements accordingly."""
        if not upload_data:
            circ_list = str(page_state['textfield'].value)
        else:
            circ_list = upload_data

        list_okay = check_circrna_input_regex(circ_list)
        db_selected = check_if_db_is_selected_local()
        page_state['circrna_found'].set_visibility(False)

        if circ_list and list_okay and db_selected:
            page_state['submit_button'].props(remove="disabled=true")
            page_state['submit_notification'].set_text(
                ui_update_found_circrnas_local(circ_list))
            return ("Submit " + str(circ_list.count('\n')) +
                    " circRNAs for ID conversion")
        elif not circ_list:
            page_state['submit_button'].props("disabled=true")
            return "Convert circRNA IDs"
        elif not list_okay:
            page_state['submit_button'].props("disabled=true")
            page_state['submit_notification'].set_text(
                "Allowed characters: A-Z,\\n,\\t,|,:,-,_ ")
            return "Unsupported characters detected in your circRNA list"
        elif not db_selected:
            page_state['submit_button'].props("disabled=true")
            return "No Input database selected"

    def handle_file_upload(file) -> None:
        """Handle uploaded file content."""
        data = file.content.decode('UTF-8')
        check_text_field_input_local(data)
        page_state['uploaded_data'] = data

    def load_example_data() -> None:
        """Load example circRNA data into the text field."""
        page_state['textfield'].value = (
            "chr1:100121447|100132793\n"
            "chr1:100121550|100122379\n"
            "chr1:100122380|100125942\n"
            "chr1:100123210|100125942\n"
            "chr1:1001237|1001530\n"
            "chr1:100123872|100125942\n"
            "chr1:100133404|100133583\n"
            "chr1:100136882|100137085\n"
            "chr1:100136882|100139711\n"
            "chr1:100136882|100140514\n"
        )

    def handle_submit() -> None:
        """Serialize local state into app.storage.tab, then navigate to results."""
        app.storage.tab['convert_state'] = {
            'mode': 'convert',
            'textfield_value': page_state['textfield'].value or '',
            'uploaded_data': page_state.get('uploaded_data'),
            'db_checkbox_value': page_state['db_checkbox'].value,
            'db_checkboxes_selected': check_if_db_is_selected_local(),
            'select2_value': page_state['select2'].value,
            'select3_value': page_state['select3'].value,
        }
        ui.navigate.to('/results')


    # ---- Build the UI ----

    with ui.card().style('width: 100%;'):
        ui.html("<strong>CircRNA ID conversion made easy</strong>")
        with ui.card_section():
            ui.markdown("""
Multiple databases gathered hundreds of thousands of circRNAs, however, most
databases employ their own naming scheme, making it harder and harder to keep
track of known circRNAs and their identifiers.
            """)

    page_state['textfield'] = ui.input(
        label='Please paste a list of circRNA IDs, one per line:',
        placeholder='start typing',
        on_change=lambda e: page_state['submit_button'].set_text(
            check_text_field_input_local(upload_data=None))
    ).props('type=textarea rows=18').style(
        "width: 100%; background-color: #ffffff;").classes('q-pa-md')

    ui.label('- OR -')

    ui.upload(
        label="1) Click + to select file "
              "2) upload file via button to the right "
              "3) press 'convert circRNA IDs' button",
        on_upload=lambda e: handle_file_upload(e.files[0])
    ).style("width: 100%")

    with ui.row():
        page_state['submit_button'] = ui.button(
            'Convert circRNA IDs',
            on_click=handle_submit
        ).props("disabled=true")
        ui.button('Load example data', on_click=load_example_data)

    page_state['circrna_found'] = ui.linear_progress(
        show_value=False, value=0).style("width: 60%;")
    page_state['circrna_found'].set_visibility(False)
    page_state['submit_notification'] = ui.label('')

    # ---- Left drawer with convert options (local to this page) ----
    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4;').props('width=390').classes(
        'q-py-none').classes('q-my-none'):

        with ui.column().classes('q-py-none').classes('q-my-none'):
            ui_layout_generate_logo()

        with ui.column().classes('q-py-none').classes('q-my-none'):
            ui.label('Step 1: select input ID type:').style(
                "text-decoration: underline;")
            page_state['db_checkbox'] = ui.select(
                ["Coordinates"] + util.select_db_columns,
                value="Coordinates",
                label="ID format").style("width: 320px")

            with ui.column().classes('q-py-none').classes('q-my-none'):
                ui.label('Step 2: select output fields:').style(
                    "text-decoration: underline;")

                db_entries = []
                ui.label("CircRNA databases & external IDs")
                third = len(util.select_db_columns) // 3
                two_thirds = (len(util.select_db_columns) * 2) // 3

                with ui.row().classes('q-py-none'):
                    with ui.column():
                        for item in range(third):
                            db_entries.append(
                                ui.checkbox(
                                    util.select_db_columns[item],
                                    value=util.active_db_columns[
                                        util.select_db_columns[item]]
                                ).props('size=xs'))
                    with ui.column():
                        for item in range(third, two_thirds):
                            db_entries.append(
                                ui.checkbox(
                                    util.select_db_columns[item],
                                    value=util.active_db_columns[
                                        util.select_db_columns[item]]
                                ).props('size=xs'))
                    with ui.column():
                        for item in range(two_thirds, len(util.select_db_columns)):
                            db_entries.append(
                                ui.checkbox(
                                    util.select_db_columns[item],
                                    value=util.active_db_columns[
                                        util.select_db_columns[item]]
                                ).props('size=xs'))

                page_state['db_checkboxes'] = db_entries

            with ui.column():
                ui.label('Step 3: select output format:').style(
                    "text-decoration: underline;")
                page_state['select2'] = ui.select(
                    {"\t": "Tab-delimited '\\t'",
                     ",": "Comma-delimited ','",
                     ";": "Semicolon-delimited ';'"},
                    value="\t",
                    label="Separator character").style("width: 320px")
                page_state['select3'] = ui.select(
                    {"NA": "NA",
                     "\t": "Tab [\\t]",
                     "": "Don't print anything"},
                    value="NA",
                    label="Placeholder for unavailable fields"
                ).style("width: 320px")

    ui_layout_add_footer_and_right_drawer()


@ui.page('/query')
async def page_application_query():
    """Query page for browsing databases.

    All mutable state is local to this function call.
    Each visitor gets their own isolated query form state."""

    # Must await connection BEFORE accessing app.storage.tab
    await ui.context.client.connected()

    ui_layout_add_head_html()
    ui_layout_add_header()

    # ---- Per-visit local state ----
    page_state = {
        'mode': 'query',
        'db_checkboxes': [],
        'query_forms': [],  # Replaces the old global ui_query_forms
        'submit_query_button': None,
    }

    # ---- Closures ----

    def check_query_text_field_local() -> None:
        """Enable/disable submit button based on query field completeness."""
        if page_state['submit_query_button'] is None:
            return

        all_good = all(
            form['query'].value
            for form in page_state['query_forms']
        )

        if all_good:
            page_state['submit_query_button'].props(remove="disabled=true")
        else:
            page_state['submit_query_button'].props("disabled=true")

        page_state['submit_query_button'].update()

    def add_conditions_local(container, new=False) -> None:
        """Add a new row of query condition inputs to the form."""
        query_values = dict()

        with container:
            with ui.row():
                if new:
                    query_values['operator1'] = ui.select(
                        ["AND", "OR", "AND NOT"],
                        value="AND",
                        label="Query operator"
                    ).style("width: 150px")

                query_values['field'] = ui.select(
                    util.db_columns,
                    value="circBase",
                    label="Database field"
                ).style("width: 130px")

                query_values['operator2'] = ui.select(
                    ["is", "LIKE", ">", "<"],
                    value="LIKE",
                    label="Query operator"
                ).style("width: 130px")

                query_values['query'] = ui.input(
                    label='Enter search term',
                    placeholder='start typing',
                    on_change=lambda e: check_query_text_field_local()
                )

        page_state['query_forms'].append(query_values)

    def remove_conditions_local(container) -> None:
        """Remove the last added query condition from the form."""
        if len(page_state['query_forms']) > 1:
            container.remove(-1)
            del page_state['query_forms'][-1]

    def check_if_db_is_selected_local() -> list:
        """Return list of selected database checkbox labels."""
        checklist = []
        for item in page_state['db_checkboxes']:
            if item.value is True:
                checklist.append(item.text)
        return checklist

    def handle_query_submit() -> None:
        """Serialize query form data into per-tab storage and navigate."""
        constraints = []
        for form in page_state['query_forms']:
            constraint = {
                'field': form['field'].value,
                'operator2': form['operator2'].value,
                'query': form['query'].value,
            }
            if 'operator1' in form:
                constraint['operator1'] = form['operator1'].value
            constraints.append(constraint)

        app.storage.tab['convert_state'] = {
            'mode': 'query',
            'constraints': constraints,
            'db_checkboxes_selected': check_if_db_is_selected_local(),
        }
        ui.navigate.to('/results')

    # ---- Build the UI ----

    initial_container = ui.column()
    add_conditions_local(initial_container, new=False)

    condition_row = ui.column()

    with ui.row():
        page_state['submit_query_button'] = ui.button(
            'Submit query',
            on_click=handle_query_submit
        ).props("disabled=false")

        ui.button(
            'Add condition',
            on_click=lambda: add_conditions_local(condition_row, new=True)
        )

        ui.button(
            'Remove condition',
            on_click=lambda: remove_conditions_local(condition_row)
        )

    circrna_found = ui.linear_progress(
        show_value=False, value=0).style("width: 60%;")
    circrna_found.set_visibility(False)

    # ---- Left drawer with output field checkboxes ----
    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4;').props('width=390').classes(
        'q-py-none').classes('q-my-none'):

        with ui.column().classes('q-py-none').classes('q-my-none'):
            ui_layout_generate_logo()

        with ui.column().classes('q-py-none').classes('q-my-none'):
            ui.label('Select output fields:').style(
                "text-decoration: underline;")

            db_entries = []
            ui.label("CircRNA databases & external IDs")
            third = len(util.select_db_columns) // 3
            two_thirds = (len(util.select_db_columns) * 2) // 3

            with ui.row().classes('q-py-none'):
                with ui.column():
                    for item in range(third):
                        db_entries.append(
                            ui.checkbox(
                                util.select_db_columns[item],
                                value=util.active_db_columns[
                                    util.select_db_columns[item]]
                            ).props('size=xs'))
                with ui.column():
                    for item in range(third, two_thirds):
                        db_entries.append(
                            ui.checkbox(
                                util.select_db_columns[item],
                                value=util.active_db_columns[
                                    util.select_db_columns[item]]
                            ).props('size=xs'))
                with ui.column():
                    for item in range(two_thirds, len(util.select_db_columns)):
                        db_entries.append(
                            ui.checkbox(
                                util.select_db_columns[item],
                                value=util.active_db_columns[
                                    util.select_db_columns[item]]
                            ).props('size=xs'))

            page_state['db_checkboxes'] = db_entries

    ui_layout_add_footer_and_right_drawer()


@ui.page('/results')
async def page_application_display_results():
    # Must await connection BEFORE accessing app.storage.tab
    await ui.context.client.connected()

    ui_layout_add_head_html()
    ui_layout_add_header()

    session_id = str(uuid4())

    # Read from per-tab storage instead of per-user storage
    user_state = app.storage.tab.get('convert_state', None)

    if user_state is None:
        ui.navigate.to('/')
        return

    mode = user_state.get('mode', '')
    output = None
    output_fields = []

    if mode == 'convert':
        output_fields = user_state.get('db_checkboxes_selected', [])
        output_fields = add_if_not_in_list(output_fields, [
            "CircRNA_ID", "Chr", "Start", "Stop", "Genome"
        ])

        if user_state.get('uploaded_data'):
            circrna_list = user_state['uploaded_data'].split('\n')
        else:
            circrna_list = user_state.get('textfield_value', '').split('\n')

        db_checkbox_value = user_state.get('db_checkbox_value', 'Coordinates')

        output = util.run_simple_select_query(
            util, output_fields, circrna_list, db_checkbox_value
        )

        if (db_checkbox_value not in output_fields and
                db_checkbox_value != "Coordinates"):
            output_fields.insert(0, db_checkbox_value)

    elif mode == 'query':
        output_fields = user_state.get('db_checkboxes_selected', [])
        output_fields = add_if_not_in_list(output_fields, [
            "CircRNA_ID", "Chr", "Start", "Stop", "Genome"
        ])

        constraints = user_state.get('constraints', [])
        sql_query = ""

        for idx, constraint in enumerate(constraints):
            if idx > 0 and 'operator1' in constraint:
                sql_query += " " + constraint['operator1'] + " "

            if constraint['operator2'] == "LIKE":
                tmp = constraint['query'].replace("*", "", 1)
                sql_query += constraint['field'] + ' LIKE "%' + tmp + '%" '
            elif constraint['operator2'] == "is":
                sql_query += constraint['field'] + ' == "' + constraint['query'] + '" '
            elif constraint['operator2'] == ">":
                sql_query += constraint['field'] + ' > "' + constraint['query'] + '" '
            elif constraint['operator2'] == "<":
                sql_query += constraint['field'] + ' < "' + constraint['query'] + '" '

        output = util.run_keyword_select_query(util, output_fields, sql_query)

    else:
        ui.navigate.to('/')
        return

    # Build result table using pure function
    separator = user_state.get('select2_value', '\t')
    empty_char = user_state.get('select3_value', 'NA')

    processed_output, table = build_result_table_for_web(
        output_fields=output_fields,
        output=output,
        mode=mode,
        separator=separator,
        empty_char=empty_char,
    )

    # Write per-session download file
    output_file = download_dir / f"{session_id}.csv"
    try:
        with output_file.open('w') as f:
            f.write(processed_output)
    except OSError:
        ui.html('<strong>Output file could not be created.</strong>')
        return

    with ui.row().classes('self-center'):
        ui.button(
            'Download table',
            on_click=lambda: ui.navigate.to(f'/download/{session_id}.csv')
        ).classes('self-center')
        ui.button(
            'New query',
            on_click=lambda: ui.navigate.to('/')
        ).classes('self-center')


@ui.page('/circrna/')
async def page_application_display_circrna_profile_missing_id():
    """Handle missing circRNA ID in URL."""
    ui_layout_add_head_html()
    ui_layout_add_header()

    with ui.left_drawer(top_corner=True, bottom_corner=False).style(
            'background-color: #d7e3f4;').props('width=390').classes(
        'q-py-none').classes('q-my-none'):
        ui.html("URL does not contain a known circRNA ID").style(
            'font-size: 24pt;')
        ui.html(
            '<br/><a href="/query/">Click here to search '
            'circhemy\'s circRNA database</a>')

    ui_layout_add_footer_and_right_drawer()


@ui.page('/circrna/{circ_id}', title=util.program_name_long)
async def page_application_display_circrna_profile(circ_id: str,
                                                   client: Client):
    """Display a circRNA profile page.

    All state is local — no globals accessed. Each visitor
    viewing a circRNA profile gets their own isolated UI context."""

    ui_layout_add_header()

    output = util.run_circrna_query(util, circrna_id=circ_id)
    output_dict = None

    if len(output) > 0:
        with ui.left_drawer(top_corner=True, bottom_corner=False).style(
                'background-color: #d7e3f4;').props('width=390').classes(
            'q-py-none').classes('q-my-none'):
            ui_layout_generate_logo()

            if len(output) > 1:
                ui.html("Multiple circRNAs matching " + "(" + str(
                    len(output)) + ")" + ":").style('font-size: 14pt;')
                ui.html("&nbsp;").style('font-size: 14pt;')

            with ui.tabs() as tabs:
                tabs.style(
                    'background-color: #d7e3f4; height: 50%; border-radius: 0px;')
                tabs.props("vertical")
                tabs.props("dense")
                tabs.props("no-caps")
                tabs.props("inline-label")
                tabs.props("outside-arrows")
                tabs.props('align="left"')

                for hit in range(0, len(output)):
                    tab_label = (str(output[hit][17]) + ":" +
                                 str(output[hit][18]) + "-" +
                                 str(output[hit][19]) + " [" +
                                 output[hit][21] + "]")
                    if output[hit][16]:
                        ui.tab(tab_label, icon='change_circle').style(
                            "justify-content: initial; "
                            "background-color: #d7e3f4;")
                    else:
                        ui.tab(tab_label, icon='donut_large').style(
                            "justify-content: initial; "
                            "background-color: #d7e3f4;")

                ui.html("&nbsp;").style('font-size: 14pt;')
                with ui.row():
                    ui.icon('change_circle').classes('text-2xl')
                    ui.label('CircRNA with new CSNv1 annotation')
                with ui.row():
                    ui.icon('donut_large').classes('text-2xl')
                    ui.html('CircRNA <u>without</u> CSNv1 annotation')

    # circRNA not found
    if len(output) == 0:
        with ui.left_drawer(top_corner=True, bottom_corner=False).style(
                'background-color: #d7e3f4;').props('width=390').classes(
            'q-py-none').classes('q-my-none'):
            ui_layout_generate_logo()

        ui.html("\"" + circ_id + "\" is an unknown circRNA ID").style(
            'font-size: 24pt;')
        ui.html(
            '<br/><a href="/query/">Click here to '
            'search circhemy\'s circRNA database</a>')

    # at least one circRNA found
    if len(output) >= 1:
        first_tab = (output[0][17] + ":" + str(output[0][18]) + "-" +
                     str(output[0][19]) + " [" + output[0][21] + "]")

        with ui.tab_panels(tabs, value=first_tab).style(
                "background-color: #f8f8f8;").classes('q-ma-none'):

            for hit in range(0, len(output)):
                output_dict = {}
                tab_label = (output[hit][17] + ":" +
                             str(output[hit][18]) + "-" +
                             str(output[hit][19]) + " [" +
                             output[hit][21] + "]")

                with ui.tab_panel(tab_label):
                    local_columns = util.all_db_columns.copy()
                    local_columns.insert(0, "Stable circhemy database ID")
                    source_database_print = ""

                    for hit_id in range(0, len(local_columns)):
                        output_dict[local_columns[hit_id]] = output[hit][hit_id]
                        if str(output[hit][hit_id]) == circ_id:
                            source_database_print = local_columns[hit_id]

                    output_dict['Coordinates'] = 1
                    output_dict['Unspliced length'] = (
                        f"{output_dict['Stop'] - output_dict['Start']:,} bp")
                    local_columns.append('Coordinates')
                    local_columns.append('Unspliced length')

                    with ui.row():
                        with ui.column():
                            with ui.column().style('width: 700px;'):
                                with ui.card().style(
                                        'width: 700px;font-size: 12pt;'). \
                                        classes('column justify-between').style(
                                    "background-color: #d7e3f4;"):
                                    with ui.row().style('height: 34px;'):
                                        with ui.column():
                                            ui.html(circ_id).style(
                                                'font-size: 24pt;')
                                        with ui.column().classes("items-center"):
                                            ui.label(source_database_print)

                            if (output_dict['CSNv1'] and
                                    source_database_print != 'CSNv1'):
                                with ui.column().style('width: 700px;'):
                                    with ui.card().style(
                                            'width: 700px;font-size: 12pt;'). \
                                            classes('column justify-between'). \
                                            style("background-color: #d7e3f4;"):
                                        with ui.row().style('height: 30px;').classes('q-pa-none'):
                                            with ui.column():
                                                ui.html(output_dict['CSNv1']).style(
                                                    'font-size: 20pt;')
                                            with ui.column().classes("items-center"):
                                                ui.label('CSNv1')

                            with ui.column().style('width: 700px;'):
                                with ui.card().style(
                                        'width: 700px;font-size: 12pt;'). \
                                        classes('column justify-between'). \
                                        style("background-color: #d7e3f4;"):
                                    for hit_id in range(0, len(local_columns)):
                                        if local_columns[hit_id] not in [
                                            "Chr", "Start", "Stop", "Strand",
                                            "Stable circhemy database ID"
                                        ]:
                                            with ui.row().style('height: 15px;'):
                                                with ui.column().style('width: 150px;'):
                                                    ui.html(
                                                        '<strong>' +
                                                        local_columns[hit_id] +
                                                        '</strong>'
                                                    ).classes('q-pa-none')
                                                with ui.column():
                                                    ui.html(
                                                        ui_generate_external_link(
                                                            local_columns[hit_id],
                                                            output_dict)
                                                    ).style('align:center').classes('q-pa-none')

                        with ui.column().style('width: 350px;'):
                            with ui.card().style(
                                    'width: 350px;font-size: 12pt;'). \
                                    classes('column justify-between'). \
                                    style("background-color: #d7e3f4;"):
                                output2 = util.get_circrna_history_by_id(
                                    util,
                                    circrna_id=output_dict['Stable circhemy database ID'])

                                with ui.column().style('width: 350px'):
                                    with ui.row():
                                        with ui.column().style('width: 220px;'):
                                            ui.html("<strong>" +
                                                    local_columns[0] +
                                                    "<strong>")
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
                                        time_code = datetime.utcfromtimestamp(
                                            output2[status_id][5]).strftime(
                                            '%m/%d/%Y')
                                        with ui.row():
                                            with ui.column().style('width: 80px;'):
                                                ui.html(str(output2[status_id][4]))
                                            with ui.column().style('width: 80px;'):
                                                ui.html(str(util.db_action_codes[
                                                                output2[status_id][1]]))
                                            with ui.column().style('width: 80px;'):
                                                ui.html(str(time_code))

    ui_layout_add_footer_and_right_drawer()

    if output_dict:
        ui_layout_add_head_html(input_dict=output_dict, circrna_id=circ_id)
    else:
        ui_layout_add_head_html()

    ui.page_title(circ_id + " - " + util.program_name_long)


# ========================================================
# ERROR PAGES
# ========================================================

@app.exception_handler(404)
async def exception_handler_404(request: Request,
                                exception: Exception) -> Response:
    with Client(ui.page(''), request=request) as client:
        with ui.column().style('width: 100%;align-items: center; gap: 0'):
            ui.label('Error 404').style(
                'font-size: 2.75rem; line-height: 1;')
            ui.html(
                "<strong>Page not found: </strong>" +
                str(request.url.path)).style(
                'font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
            ui.html("<a href=\"/\"><img src=\"/static/error.png\"></a>").style(
                'text-align: center; padding:10px;')
            ui.label(
                'Please contact the server administrator at ' +
                util.support_email + ' for additional support.').style(
                'font-size: 1.25rem; line-height: 1.75rem;')
            ui.link(
                'Click here to report an issue on GitHub.',
                util.support_web + "?title=circhemy 404 error" +
                "&body=URL: " + str(request.url) +
                "%0A%0APlease describe how the error occurred").style(
                'font-size: 1.25rem; line-height: 1.75rem;')
    return client.build_response(request, 404)


@app.exception_handler(500)
async def exception_handler_500(request: Request,
                                exception: Exception) -> Response:
    with Client(ui.page(''), request=request) as client:
        with ui.column().style('width: 100%;align-items: center; gap: 0'):
            ui.label('Internal server error').style(
                'font-size: 2.75rem; line-height: 1;')
            ui.html("<strong>Exception:</strong>").style(
                'font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
            ui.html(str(exception)).style(
                'font-size: 1.25rem; line-height: 1.75rem;')
            ui.html("<a href=\"/\"><img src=\"/static/error.png\"></a>").style(
                'text-align: center; padding:10px;')
            ui.label(
                'Please contact the server administrator at ' +
                util.support_email + ' for additional support.').style(
                'font-size: 1.25rem; line-height: 1.75rem;')
            ui.link(
                'Click here to report an issue on GitHub.',
                util.support_web + "?title=circhemy exception: " +
                str(exception) + "&body=URL: " + str(request.url) +
                "%0A%0APlease describe how the error occurred").style(
                'font-size: 1.25rem; line-height: 1.75rem; padding: 1.25rem 0')
    return client.build_response(request, 500)


# ========================================================
# REST API — Stateless by design, already multi-user safe
# ========================================================

class ConvertModel(BaseModel):
    input: str
    output: List[str]
    query: List[str]

    @field_validator('query')
    @classmethod
    def circrna_id_pattern_check(cls, v):
        if len(v) == 0:
            raise ValueError('No CircRNA IDs for conversion provided.')
        for item in v:
            if not check_circrna_input_regex(item):
                raise ValueError('CircRNA IDs contains illegal characters.')
            if item == "":
                raise ValueError('No CircRNA IDs for conversion provided.')
        return v

    @field_validator('input')
    @classmethod
    def database_name_check_input(cls, v):
        fields_allowed = list(util.external_db_urls.keys())
        if v not in fields_allowed:
            raise ValueError("Unsupported input field provided."
                             " Supported fields are: "
                             + ', '.join(fields_allowed) +
                             ". Field names are case-sensitive.")
        return v

    @field_validator('output')
    @classmethod
    def database_name_check_output(cls, v):
        fields_allowed = list(util.external_db_urls.keys())
        for field in v:
            if field not in fields_allowed:
                raise ValueError("Unsupported output field provided."
                                 " Supported fields are: "
                                 + ', '.join(fields_allowed) +
                                 ". Field names are case-sensitive.")
        return v


class ConstraintModel(BaseModel):
    query: str
    field: str
    operator1: str
    operator2: str

    @field_validator('query')
    @classmethod
    def circrna_id_pattern_check(cls, v):
        if not check_circrna_input_regex(v):
            raise ValueError('Query contains illegal characters.')
        return v

    @field_validator('field')
    @classmethod
    def field_name_check(cls, v):
        fields_allowed = list(util.external_db_urls.keys())
        if v not in fields_allowed:
            raise ValueError("Unsupported input field provided."
                             " Supported fields are: "
                             + ', '.join(fields_allowed) +
                             ". Field names are case-sensitive.")
        return v

    @field_validator('operator1')
    @classmethod
    def operator1_check(cls, v):
        fields_allowed = ["AND", "OR", "AND NOT"]
        if v not in fields_allowed:
            raise ValueError("Unsupported operator1 provided."
                             " Supported fields are: "
                             + ', '.join(fields_allowed) +
                             ". Field names are case-sensitive.")
        return v

    @field_validator('operator2')
    @classmethod
    def operator2_check(cls, v):
        fields_allowed = ["IS", "LIKE", ">", "<"]
        if v not in fields_allowed:
            raise ValueError("Unsupported operator2 provided."
                             " Supported fields are: "
                             + ', '.join(fields_allowed) +
                             ". Field names are case-sensitive.")
        return v


class QueryModel(BaseModel):
    input: List[ConstraintModel]
    output: List[str]

    @field_validator('output')
    @classmethod
    def database_name_check_output(cls, v):
        fields_allowed = list(util.external_db_urls.keys())
        for field in v:
            if field not in fields_allowed:
                raise ValueError("Unsupported output field provided."
                                 " Supported fields are: "
                                 + ', '.join(fields_allowed) +
                                 ". Field names are case-sensitive.")
        return v


@app.post("/api/convert")
async def process_api_convert_call(data: ConvertModel):
    """REST API convert endpoint — stateless, no globals."""
    _, table = build_result_table_for_api(data.input, data.output, data.query)
    return table


@app.post("/api/query")
async def process_api_query_call(data: QueryModel):
    """REST API query endpoint — stateless, no globals."""
    _, table = build_result_table_for_api(data.input, data.output)
    return table
