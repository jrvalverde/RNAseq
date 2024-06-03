#!/bin/bash

#org='cjo'
#org='gga'

# Start with GO

cat > index.html << END
<HTML>
<head>
   <title>ClusterProfiler GO and KEGG Analysis</title>
   <link rel="stylesheet" href="./style.css" />
</head>
<BODY>
<center><h1>ClusterProfiler GSEA Analysis</h1></center>

<p>Index</p>
<ul>
    <li><a href="#go">Enrichment analysis uning Gene Ontologies</a></li>
    <ul>
        <li><a href="#gographs">General graphs</a></li>
		<li><a href="#gotable">detailed table</a></li>
    </ul>
    <li><a href="#kegg">Enrichment analysis using KEGG Pathways</a></li>
    <ul>
        <li><a href="#kegggraphs">General graphs</a></li>
		<li><a href="#keggtable">Detailed table</a></li>
    </ul>
</ul>


<a id="go"></a>
<center><h2>Enrichment analysis using Gene Ontologies</h2></center>

END


if [ -e cProf.topGO.tab ] ; then
  base=cProf
else
  base=cProf.raw_p
  echo "<h3>NOTE: raw probabilities were used (not p.adjusted)!</h3>" >> index.html
  echo "<p>This probably means that no significantly regulated gene <b>sets</b> were found.</p>" >> index.html
fi


cat >> index.html << END

<h3>Summary plots</h3>
<center>
<a id="gographs"></a>
<a href="${base}.GOdotplot.png"><img src="${base}.GOdotplot.png" width="256" height="256" alt="dotplot"></img></a>
<a href="${base}.GOemmaplot.png"><img src="${base}.GOemmaplot.png" width="256" height="256" alt="emmaplot"></img></a>
<a href="${base}.GOridgeplot.png"><img src="${base}.GOridgeplot.png" width="256" height="256" alt="ridgeplot"></img></a>
</center>
<p>(click to enlarge)</p>

<a id="gotable"></a>
<center><h3>Top GOs</h3></center>
<p>Click on Graph to enlarge</p>
<table border="1">
  <colgroup>
      <col span="1" style="width: 5%">
      <col span="1" style="width: 5%">
      <col span="1" style="width: 10%">
      <col span="1" style="width: 30%">
      <col span="1" style="width: 5%">
      <col span="1" style="width: 5%">
      <col span="1" style="width: 5%">
      <col span="1" style="width: 30%">
  </colgroup>
<tr>
	<td>ID</td>
	<td>Graph</td>
	<td>TYPE</td>
	<td>DESCRIPTION</td>
	<td>N.members</td>
	<td>ES</td>
	<td>NES</td>
	<td>p.adjust</td>
	<td>core enrichment</td>
</tr>
END

i=0
cat ${base}.topGO.tab \
| tail -n +2 \
| while IFS=$'\t' read \
	name ontology ID description setSize ES NES \
	pvalue padjust qvalue rank leading_edge core_enrichment ; do
    
    if [ "$ontology" == '"BP"' ] ; then ontology="Biological Process"
    elif [ "$ontology" == '"MF"' ] ; then ontology="Molecular Function"
    elif [ "$ontology" == '"CC"' ] ; then ontology="Cellular Component"
    fi
    
    ID=`echo $ID | tr -d '"'`	# remove quotes
    
    (( i++ ))
    graph="./${base}.GOcnetplot.$i.$ID.png"
    echo " <tr>" \
         "<td>$ID</td>" \
         "<td><a href=\"$graph\"><img src=\"$graph\" width=\"64\" height=\"64\"></img></a></td>" \
	 "<td>$ontology</td>" \
	 "<td>$description</td>" \
	 "<td align=\"right\">$setSize</td>" \
	 "<td align=\"right\">" `/usr/bin/printf "%0.4g" "$ES"` "</td>" \
	 "<td align=\"right\">" `/usr/bin/printf "%0.4g" "$NES"` "</td>" \
	 "<td align=\"right\">" `/usr/bin/printf "%0.4g" "$padjust"` "</td>" \
	 "<td>$core_enrichment</td>" \
	 "</tr>"
done >> index.html

#-----------------------------------------------------------------------------#

# Proceed to KEGG

cat >> index.html <<END
  </table>

<a id="kegg"></a>
<center><h2>Enrichment analysis using KEGG pathways</h2></center>
END

if [ -e cProf.topKEGG.tab ] ; then
  base=cProf
else
  base=cProf.raw_p
  echo "<h3>NOTE: raw probabilities were used (not p.adjusted)!</h3>" >> index.html
  echo "<p>This probably means that no significantly regulated gene <b>sets</b> were found.</p>" >> index.html
fi

cat >> index.html <<END

<center>
<a id="kegggraphs"></a>
<a href="${base}.KEGGdotplot.png"><img src="${base}.KEGGdotplot.png" width="256" height="256" alt="dotplot"></img></a>
<a href="${base}.KEGGemmaplot.png"><img src="${base}.KEGGemmaplot.png" width="256" height="256" alt="emmaplot"></img></a>
<a href="${base}.KEGGridgeplot.png"><img src="${base}.KEGGridgeplot.png" width="256" height="256" alt="ridgeplot"></img></a>
<a href="${base}.KEGGcnetplot.png"><img src="${base}.KEGGcnetplot.png" width="256" height="256" alt="cnetplot"></img></a>
<p>(click to enlarge)</p>
</center>

<a id="keggtable"></a>
<center><h3>Top KEGG Pathways</h3></center>

<p> Click on graph to enlarge. Click on pathway name to see as a PDF file</p>

<table border="1">
  <colgroup>
      <col span="1" style="width: 5%">
      <col span="1" style="width: 5%">
      <col span="1" style="width: 30%">
      <col span="1" style="width: 5%">
      <col span="1" style="width: 5%">
      <col span="1" style="width: 5%">
      <col span="1" style="width: 5%">
      <col span="1" style="width: 40%">
  </colgroup>
	<tr>
  		<td>ID</td>
		<td>Graph</td>
		<td>DESCRIPTION</td>
		<td>N.members</td>
		<td>ES</td>
		<td>NES</td>
		<td>p.adjust</td>
		<td>core enrichment</td>
	</tr>
END

cat ${base}.topKEGG.tab \
| tail -n +2 \
| while IFS=$'\t' read \
	name ID description setSize ES NES \
	pvalue padjust qvalue rank leading_edge core_enrichment ; do
    
    if [ "$ontology" == "BP" ] ; then ontology="Biological Process"
    elif [ "$ontology" == "MF" ] ; then ontology="Molecular Function"
    elif [ "$ontology" == "CC" ] ; then ontology="Cellular Component"
    fi
    
    ID=`echo $ID | tr -d '"'`	# remove quotes

    echo -n " <tr>" 
    if [ -e "$ID.pathview.pdf" ] ; then
        echo -n "<td><a href=\"$ID.pathview.pdf\">$ID</a></td>" 
    else
        echo -n "<td>$ID</td>"
    fi
    echo "<td><a href=\"$ID.pathview.png\"><img src=\"$ID.pathview.png\" width=\"64\" height=\"64\"></img></a></td>" \
	 "<td>$description</td>" \
	 "<td align=\"right\">$setSize</td>" \
	 "<td align=\"right\">" `/usr/bin/printf "%0.4g" "$ES"` "</td>" \
	 "<td align=\"right\">" `/usr/bin/printf "%0.4g" "$NES"` "</td>" \
	 "<td align=\"right\">" `/usr/bin/printf "%0.4g" "$padjust"` "</td>" \
	 "<td>$core_enrichment</td>" \
	 "</tr>"
done >> index.html



cat >> index.html <<END
  </table>
END
exit



cat >> index.html <<END
  </table>
</td></tr>
</table>
</BODY>
</HTML>
END

