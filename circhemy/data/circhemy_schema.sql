BEGIN TRANSACTION;
DROP TABLE IF EXISTS `circhemy`;
CREATE TABLE IF NOT EXISTS `circhemy` (
	`Species`	TEXT,
    `Gene`	TEXT,
    `Description` TEXT,
    `ENSEMBL` TEXT,
    `Entrez` INTEGER,
  	`circBase`	TEXT,
	`CircAtlas2`	TEXT,
	`circRNADb`	TEXT,
	`deepBase2`	TEXT,
	`Circpedia2`	TEXT,
	`circBank`	TEXT,
	`riboCIRC`	TEXT,
	`exoRBase2`	TEXT,
	`Arraystar`	TEXT,
	`CSNv1`	TEXT,
	`Chr`	TEXT,
	`Start`	INTEGER,
	`Stop`	INTEGER,
	`Strand`	TEXT,
	`Genome`	TEXT,
	`Pubmed`	INTEGER
);
COMMIT;
