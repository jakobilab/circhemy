#!/bin/bash

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

if ! command -v sqlite3 &> /dev/null
then
    echo "SQLite3 is required, but not installed or in \$PATH."
    (exit 1)
fi


echo "Loading base SQL table schema."

sqlite3 ../data/circconvert.sqlite3 < ../data/circrnadb_schema.sql

echo "SQL database schema ready, starting data import."

sqlite3 ../data/circconvert.sqlite3 << EOF
.mode csv
.separator "\t"
.import ../data/circrnadb_data.csv circrnadb
VACUUM;
EOF

echo "Data import finished."
