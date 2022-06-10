CREATE TABLE papers AS
SELECT u.id, u.institution, u.url_stem, u.date_created, u.current, u.url_scraped, p.authors, p.orcids, p.author_urls, p.subjects
FROM papers_old p LEFT JOIN urls u
ON p.id=u.id;
