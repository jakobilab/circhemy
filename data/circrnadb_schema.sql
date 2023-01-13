BEGIN TRANSACTION;
DROP TABLE IF EXISTS `circrnadb`;
CREATE TABLE IF NOT EXISTS `circrnadb` (
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
	`Chr`	TEXT,
	`Start`	INTEGER,
	`Stop`	INTEGER,
	`Strand`	TEXT,
	`Genome`	TEXT,
	`Pubmed`	INTEGER
);
COMMIT;
