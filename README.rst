**circhemy**
======================================================================

**The alchemy of circular RNA ID conversion**

.. image:: https://github.com/jakobilab/circhemy/raw/main/circhemy/web/static/logo2.png
   :alt: circhemy - The alchemy of circular RNA ID conversion

|downloads| |pypi|

Introduction
-------------

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

Circhemy
-------------

We developed circhemy, a modular, Python3-based framework for circRNA ID
conversion that unifies several functionalities in a single Python package.
Three different routes are implemented within package to access more than 2
million circRNA IDs:

* User-friendly web application at `circhemy.jakobilab.org <https://circhemy.jakobilab.org>`__
* Streamlined CLI application for direct access to the prepackaged local SQLite3 database
* A public `REST API <https://circhemy.jakobilab.org/rest/>`__ that enables direct access to the most recent ID database from HPC systems using curl or similar tools

Circhemy includes two different modes of action: ``convert`` and ``query``. Convert
allows the user to convert from one type of circRNA ID to a wide variety of
other database identifiers, while query allows users to run direct queries on
the circRNA database to extract circRNAs fulfilling a user-defined set of
constraints.

Installation
-------------

The circhemy CLI package is written in Python3 (>=3.7) and consists of two
core modules, namely ``convert`` and ``query``. The command line version requires
only one external dependency, ``sqlite3``, for access to the internal SQLite3
database with circRNA ID data

Installation is managed through ``pip3 install circhemy`` or ``python3 setup.py
install`` when installed from the cloned GitHub repository. No sudo access is
required if the installation is executed with ``--user`` which will install the
package in a user-writeable folder. The binaries should be installed
to ``/home/$user/.local/bin/`` in case of Debian-based systems.

circhemy was developed and tested on Debian Buster, but should run with
any other distribution.

The latest release version of circhemy can be installed via pip:

.. code-block:: console

    pip3 install circhemy

Additionally, this repository offers the latest development version:

.. code-block:: console

    pip3 install git+https://github.com/jakobilab/circhemy.git



Command Line Interface
-----------------------

Circhemy currently offers two modules:

Convert module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The convert module is able to convert from a range of input circRNA ID into different one or more database identifiers.

Example: Convert a list of CircAtlas IDs read via STDIN from file input.csv into circPedia2 IDs, but also output  CircAtlas IDs, while writing the output to /tmp/output.csv:

.. code-block:: console

    cat input.csv | circhemy convert -q STDIN -i circatlas -o circpedia2 circatlas -O /tmp/output.csv

Query module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The query module is able to retrieve circRNA IDs from the internal database that fulfil a set of user-defined constraints.

Example: Retrieve a list of circbase and circatlas circRNA IDs that are located on chromosome 3 of the species rattus norvegicus; only print out circRNAs from the rn6 genome build.

.. code-block:: console

    circhemy query -o circbase circatlas -C chr3 -s rattus_norvegicus -g rn6


Representational State Transfer Interface (REST)
-------------------------------------------------

Representational State Transfer, or REST for short, allows users and software
developers to easily access circhemy from within their own tools or pipelines.
Circhemy's REST API uses JSON for input queries and returning output, making it
easy to format queries from every programming language or even by hand.

The REST API it publicly available and uses a fixed set of keywords to perform
conversions or queries. Two examples for the two different modes of action are
shown below.

Convert module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The convert module is able to convert from a range of input circRNA ID into
different one or more database identifiers.

Example: Convert a list of CircAtlas IDs into circBase and
into CircPedia2 IDs, but also output CircAtlas IDs.

.. code-block:: console

    curl -X 'POST' 'https://circhemy.jakobilab.org/api/convert'
      -H 'accept: application/json'
      -H 'Content-Type: application/json'
      -d '{
          "input": "CircAtlas",
          "output": ["CircPedia2","CircAtlas"],
          "query": ["hsa-MYH9_0004","hsa-MYH9_0004"]
          }'

Output is returned as JSON-formatted string which can directly be used for AG
Grid tables for any other postprocessing:

.. code-block:: json

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

Query module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The query module is able to retrieve circRNA IDs from the internal database that
fulfil a set of user-defined constraints.

Example: Retrieve all circRNAs with a CircAtlas ID containing *atf* or *xbp1*,
return the IDs in circBase and circAtlas format:

.. code-block:: console

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

Output is returned as JSON-formatted string which can directly be used for AG
Grid tables for any other postprocessing:

.. code-block:: json

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


.. |downloads| image:: https://pepy.tech/badge/circtools
    :alt: Python Package Index Downloads
    :scale: 100%
    :target: https://pepy.tech/project/circhemy

.. |pypi| image:: https://badge.fury.io/py/circtools.svg
    :alt: Python package version
    :scale: 100%
    :target: https://badge.fury.io/py/circhemy


About
-------------
Circhemy is developed at the `Jakobi Lab <https://jakobilab.org/>`__, part of
the `Translational Cardiovascular Research Center (TCRC) <https://phoenixmed.arizona.edu/tcrc/>`__, in the Department of Internal Medicine at `The University of Arizona College of Medicine – Phoenix <https://phoenixmed.arizona.edu/>`__.

Contact: **circhemy@jakobilab.org**
