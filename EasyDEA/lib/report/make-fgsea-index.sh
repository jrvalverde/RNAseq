topUp=`cat topUp* | wc -l`
topDn=`cat topDn* | wc -l`


cat > index.html <<END
<HTML>
<head>
   <title>FGSEA GO Analysis</title>
   <link rel="stylesheet" href="./style.css" />
</head>
<BODY>
<center><h1>Gene Set Enrichment Analysis with FGSEA</h1></center>

<H2>Barplot</H2>
<center><a href="fgsea.barplot.png"><img src="fgsea.barplot.png" width="128" height="128"></img></a>
<p>(click to enlarge)</p></center>

<h2>Citations of the topmost up/down regulated ontologies in PubMed</h2>

<center><a href="fgsea.EUPMC.png"><img src="fgsea.EUPMC.png" width="128" height="128"></img></a>
<p>(click to enlarge)</p></center>

<H2>Graphic Summary</H2>

<center><a href="fgsea.GSEAtable.png"><img src="fgsea.GSEAtable.png" width="256" height="256"></img></a>
<p>(click to enlarge)</p></center>

<p>The leading edge subset of a gene set is the subset of members that
is most affected. For a positive ES,
the leading edge subset is the set of members that appear in the ranked list
prior to the peak score. For a negative ES, it is the set of members that
appear subsequent to the peak score. The plots help you get an idea of how
many genes in the pathway/ontology group are the most affected.</p>

<p><a href="http://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Enrichment_Score_(ES)">Help on interpreting the graphic</a></p>

<table>
<tr><td>

  <center><H2>Top $topUp up-regulated GOs</H2></center>
  <table border="1">
  <tr>
    <th>No.</th>
    <th>Description</th>
    <th>Graph<br>(Click)</th>
  </tr>
END

if [ ! -s topUp.*.txt ] ; then
    echo "ERROR:"
    pwd
    echo "topUp$topUp.txt is empty or does not exist"
    exit
fi

cat topUp.*.txt \
| while read n go ; do
    echo "    <tr>"
    echo "      <td>$n</td>"
    echo "      <td>$go</td>"
    echo -n '      <td><a href="fgsea.topUp.'
      printf "%03d" $n
      echo -n '.enrichment.png">'
    echo -n '<img src="fgsea.topUp.'
      printf "%03d" $n
      echo '.enrichment.png" width="32" height="32" alt="enrichment plot"></img></a></td>'
    echo "    </tr>"
done >> index.html


cat >> index.html <<END
  </table>
</td>
<td>
  <center><H2>Top $topDn Down-regulated GOs</H2></center>
  <table border="1">
  <tr>
    <th>No.</th>
    <th>Description</th>
    <th>Graph<br>(Click)</th>
  </tr>
END


if [ ! -s topDn.*.txt ] ; then
    echo "ERROR:"
    pwd
    echo "topDn.$topDn.txt is empty or does not exist"
    exit
fi

cat topDn.*.txt \
| while read n go ; do
    echo "    <tr>"
    echo "      <td>$n</td>"
    echo "      <td>$go</td>"
    echo -n '      <td><a href="fgsea.topDn.'
      printf "%03d" $n
      echo -n '.enrichment.png">'
    echo -n '<img src="fgsea.topDn.'
      printf "%03d" $n
      echo '.enrichment.png" width="32" height="32" alt="enrichment plot"></img></a></td>'
    echo "    </tr>"
done >> index.html



cat >> index.html <<END
  </table>
</td></tr>
</table>
</BODY>
</HTML>
END

