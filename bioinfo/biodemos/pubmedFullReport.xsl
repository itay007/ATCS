<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:variable name="article" select="PubmedArticleSet/PubmedArticle" />
<xsl:variable name="totalNumArticles" select="count($article)" />

<xsl:template match="/">

<html>
<head>
<title>PubMed Query</title>
<link href="/css/default.css" type="text/css" rel="stylesheet" />
</head>
<body>
<h2>PubMed Query</h2>
<h3>Total Number of Unique Articles: <xsl:value-of select="$totalNumArticles" /></h3>
<table border="1" width="100%">
<tr>
<th>Publication Date</th>
<th>Journal</th>
<th>Title</th>
<th>Link to Abstract</th>
</tr>
<xsl:for-each select="PubmedArticleSet/PubmedArticle">
<xsl:sort select="MedlineCitation/DateCreated/Year"/>
<tr>

<!--Publication Date-->
<td valign="middle">
<xsl:value-of select="MedlineCitation/DateCreated/Year"/>
</td>

<!--Journal-->
<td valign="middle">
<xsl:value-of select="MedlineCitation/MedlineJournalInfo/MedlineTA"/>
</td>
<!-- Title -->
<td>
<xsl:value-of select="MedlineCitation/Article/ArticleTitle"/>
</td>

<!-- Link to abstract -->
<td valign="middle">
<a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&#38;cmd=Retrieve&#38;dopt=AbstractPlus&#38;query_hl=12&#38;itool=pubmed_docsum&#38;list_uids={MedlineCitation/PMID}"><xsl:value-of select="MedlineCitation/PMID"/></a>
</td>

</tr>
</xsl:for-each>
</table>

</body>

</html>

</xsl:template>
</xsl:stylesheet>