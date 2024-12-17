-- Select Foldseek AFDB clusters where >=20 members match the
-- CATH domain architecture (up to T level) of the representative.
-- Expects connection to db with ted_domain table (from ted_domains_first_split.tsv)
-- then attach the afdb clusters db.


-- Table with CAT and CATH domain strings
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

-- Foldseek-clustered members
CREATE TEMP TABLE IF NOT EXISTS filtered_am AS
	SELECT accession, rep_accession
	FROM afdb.member
	WHERE flag = 2;

-- Table with structures with unique domain architectures (i.e. unique H)
-- Group by CATH, then take the first
CREATE TABLE IF NOT EXISTS unique_members AS
	SELECT MIN(filtered_am.accession) as accession
	FROM filtered_am
	INNER JOIN cath_annotations ca1 on filtered_am.accession == ca1.accession
	INNER JOIN cath_annotations ca2 on filtered_am.rep_accession == ca2.accession
	WHERE ca1.cat = ca2.cat
	GROUP BY filtered_am.rep_accession, ca1.cath;
CREATE INDEX IF NOT EXISTS idx_unique_members_accession ON unique_members(accession);

-- Filter AFDB database for unique structures only
CREATE TEMP TABLE IF NOT EXISTS unique_am AS
	SELECT *
	FROM afdb.member
	WHERE accession IN unique_members;

-- Re-join CATH annotation data
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

-- Select clusters with >=2 domains and >=20 members
SELECT rep_accession, rep_cath, COUNT(mem_accession) AS num_members, num_domains, GROUP_CONCAT(mem_accession)
FROM cath_clusters
GROUP BY rep_accession
HAVING num_domains >= 2 AND num_members >=20
ORDER BY num_domains DESC;