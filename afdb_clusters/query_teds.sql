-- Select Foldseek AFDB clusters where >=20 members match the
-- CATH domain architecture (up to T level) of the representative.
-- Expects connection to db with ted_domain table (from ted_domains_first_split.tsv)
-- then attach the afdb clusters db.


-- Match on CAT, mismatch on CATH
CREATE TABLE IF NOT EXISTS cath_annotations AS
	SELECT
		accession,
		COUNT(1) AS cnt,
		GROUP_CONCAT(c || '.' || a || '.' || t || '.' || h) AS cath,
		GROUP_CONCAT(c || '.' || a || '.' || t) AS cat
	FROM ted_domain
	GROUP BY accession
	HAVING cnt >= 2;

CREATE INDEX IF NOT EXISTS idx_cath_annotations_accession ON cath_annotations(accession);
CREATE INDEX IF NOT EXISTS idx_cath_annotations_cat ON cath_annotations(cat);

CREATE TEMP TABLE IF NOT EXISTS filtered_am AS
	SELECT accession, rep_accession
	FROM afdb.member
	WHERE flag = 2;

CREATE TABLE IF NOT EXISTS unique_members AS
	SELECT MIN(filtered_am.accession) as accession
	FROM filtered_am
	INNER JOIN cath_annotations ca1 on filtered_am.accession == ca1.accession
	INNER JOIN cath_annotations ca2 on filtered_am.rep_accession == ca2.accession
	WHERE ca1.cat = ca2.cat
	GROUP BY filtered_am.rep_accession, ca1.cath;

CREATE INDEX IF NOT EXISTS idx_unique_members_accession ON unique_members(accession);

CREATE TEMP TABLE IF NOT EXISTS unique_am AS
	SELECT *
	FROM afdb.member
	WHERE accession IN unique_members;

CREATE TABLE IF NOT EXISTS cath_clusters AS 
	SELECT
		unique_am.accession AS mem_accession,
		unique_am.rep_accession AS rep_accession,
		ca2.cnt as num_domains,
		ca2.cath AS rep_cath,
		ca1.cath AS mem_cath
	FROM unique_am
	INNER JOIN cath_annotations ca1 on unique_am.accession == ca1.accession
	INNER JOIN cath_annotations ca2 on unique_am.rep_accession == ca2.accession;

-- SELECT am.rep_accession, ca2.cath, COUNT(1) as num_members, ca2.cnt as num_domains, GROUP_CONCAT(am.accession)
-- FROM afdb.member am
-- INNER JOIN cath_annotations ca1 on am.accession == ca1.accession
-- INNER JOIN cath_annotations ca2 on am.rep_accession == ca2.accession
-- WHERE am.flag = 2 AND ca1.cat = ca2.cat AND ca1.cath != ca2.cath
-- GROUP BY am.rep_accession
-- HAVING num_members >= 20
-- ORDER BY num_domains DESC;
