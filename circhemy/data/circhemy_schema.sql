BEGIN TRANSACTION;
DROP TABLE IF EXISTS `circhemy`;
CREATE TABLE IF NOT EXISTS `circhemy` (
	`Species`	TEXT,
	`circBase`	TEXT,
	`CircAtlas`	TEXT,
	`circRNADB`	TEXT,
	`Deepbase2`	TEXT,
	`Circpedia2`	TEXT,
	`circBank`	TEXT,
	`riboCIRC`	TEXT,
	`exorBase2`	TEXT,
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