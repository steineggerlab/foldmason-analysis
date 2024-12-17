-- Select clusters with >=20 members
-- Requires database built with query_teds.sql
SELECT am.rep_accession, ca2.cath, COUNT(1) as num_members, ca2.cnt as num_domains, GROUP_CONCAT(am.accession)
FROM afdb.member am
INNER JOIN cath_annotations ca1 on am.accession == ca1.accession
INNER JOIN cath_annotations ca2 on am.rep_accession == ca2.accession
WHERE am.flag = 2 AND ca1.cat = ca2.cat AND ca1.cath != ca2.cath
GROUP BY am.rep_accession
HAVING num_members >= 20
ORDER BY num_domains DESC;
