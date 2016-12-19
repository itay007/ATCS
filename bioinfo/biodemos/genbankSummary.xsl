<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:template match="/">

<html>
<head>
<title>GenBank Query</title>
<link href="/css/default.css" type="text/css" rel="stylesheet" />
</head>

<body>
<h2>GenBank Query Results</h2>

<table border="3">
<tr>
<th>GenBank Accession</th>
<th>Sequence Description</th>
<th>ID</th>
<th>Link to Full Report</th>
</tr>


<xsl:for-each select="eSummaryResult/DocSum">
<xsl:sort select="Item[@Name='Caption']"/>
<tr>
<!--GenBank Accession-->
<td valign="middle">
<xsl:value-of select="Item[@Name='Caption']"/>
</td>
<!-- Sequence Description -->
<td>
<xsl:value-of select="Item[@Name='Title']"/>
</td>
<!-- ID -->
<td>
<xsl:value-of select="Item[@Name='Extra']"/>
</td>

<!-- Link to Full Report -->
<td valign="middle">
<a href="http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&#38;val={Item[@Name='Gi']}"><xsl:value-of select="Item[@Name='Gi']"/></a>
</td>

</tr>
</xsl:for-each>
</table>
</body>

</html>

</xsl:template>
</xsl:stylesheet>

