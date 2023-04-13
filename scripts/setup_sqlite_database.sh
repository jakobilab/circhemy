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

sqlite3 ../circhemy/data/circhemy.sqlite3 < ../circhemy/data/circhemy_schema.sql

echo "SQL database schema ready, starting data import."

commandfile=$(mktemp)

# create temporary init script
cat <<EOF > $commandfile
.mode csv
.separator "\t"
.import /dev/stdin circhemy
UPDATE circhemy SET Species= NULL WHERE Species= 'NA';
UPDATE circhemy SET Gene= NULL WHERE Gene= 'NA';
UPDATE circhemy SET Entrez= NULL WHERE Entrez= 'NA';
UPDATE circhemy SET ENSEMBL= NULL WHERE ENSEMBL= 'NA';
UPDATE circhemy SET Description= NULL WHERE Description= '';
UPDATE circhemy SET Description= NULL WHERE Description= 'NA';
UPDATE circhemy SET circBase= NULL WHERE circBase= 'NA';
UPDATE circhemy SET CircAtlas2= NULL WHERE CircAtlas2= 'NA';
UPDATE circhemy SET circRNADb= NULL WHERE circRNADb= 'NA';
UPDATE circhemy SET deepBase2= NULL WHERE deepBase2= 'NA';
UPDATE circhemy SET Circpedia2= NULL WHERE Circpedia2= 'NA';
UPDATE circhemy SET circBank= NULL WHERE circBank= 'NA';
UPDATE circhemy SET riboCIRC= NULL WHERE riboCIRC= 'NA';
UPDATE circhemy SET exoRBase2= NULL WHERE exoRBase2= 'NA';
UPDATE circhemy SET Arraystar= NULL WHERE Arraystar= 'NA';
UPDATE circhemy SET CSNv1= NULL WHERE CSNv1= 'NA';
UPDATE circhemy SET Chr=NULL WHERE Chr= 'NA';
UPDATE circhemy SET Start=NULL WHERE Start= 'NA';
UPDATE circhemy SET Stop=NULL WHERE Stop= 'NA';
UPDATE circhemy SET Strand=NULL WHERE Strand= 'NA';
UPDATE circhemy SET Genome=NULL WHERE Genome= 'NA';
UPDATE circhemy SET Pubmed=NULL WHERE Pubmed= 'NA';
VACUUM;
EOF

bzip2 -d -c ../circhemy/data/circhemy_data.csv.bz2 | sqlite3 --init "$commandfile" ../circhemy/data/circhemy.sqlite3

echo "Data import finished, creating bzipped2 file for deployment"
bzip2 -f -k --best ../circhemy/data/circhemy.sqlite3
