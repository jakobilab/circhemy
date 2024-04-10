BEGIN TRANSACTION;
DROP TABLE IF EXISTS `circhemy`;
CREATE TABLE IF NOT EXISTS `circhemy` (
    `CircRNA_ID` INTEGER PRIMARY KEY,
	`Species` TEXT,
    `Gene` TEXT,
    `Description` TEXT,
    `ENSEMBL` TEXT,
    `Entrez` INTEGER,
  	`circBase` TEXT,
  	`circBase_alt` TEXT,
	`CircAtlas2` TEXT,
	`circRNADb`	TEXT,
	`deepBase2`	TEXT,
	`Circpedia2` TEXT,
	`circBank` TEXT,
	`riboCIRC` TEXT,
	`exoRBase2` TEXT,
	`Arraystar` TEXT,
	`CSNv1` TEXT,
	`Chr` TEXT NOT NULL,
	`Start` INTEGER NOT NULL,
	`Stop` INTEGER NOT NULL,
	`Strand` TEXT,
	`Genome` TEXT NOT NULL,
	`Pubmed` INTEGER
);
CREATE TABLE IF NOT EXISTS `circhemy_db_info` (
    `DB_ID` INTEGER PRIMARY KEY,
	`Version` TEXT NOT NULL UNIQUE,
    `Date` INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS `circhemy_log` (
    `CircRNA_ID` INTEGER PRIMARY KEY,
    `Action` INTEGER,
    `DB_ID`	INTEGER
);
COMMIT;
