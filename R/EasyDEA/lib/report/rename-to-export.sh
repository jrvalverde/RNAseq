prefix="sample_"
OLDWAY="FALSE"

cp -Ru max.size=500/* max_grp_size_500/

cd max_grp_size_500

for i in sample* ; do
    if [ ! -d $i ] ; then continue ; fi
    newdir=${i#sample_}
    cd $newdir
    for f in * ; do
        name=${f%.*}
        ext=${f##*.}
        nname=${f//./_}
        mv $name.$ext $nname.$ext
        pwd
        if [ -x ../../mk_index_export.sh ] ; then
        ../../mk_index_export.sh
        fi
        echo "$i done"
    done
    cd ..
done


if [ "$OLDWAY" = "TRUE" ] ; then
    mv max.size=500 max_grp_size_500
    for i in max_grp_size_500/* ; do mv $i ${i//./_}; done

    for i in */* ; do mv $i ${i//./_}; done

    for i in */*_png ; do mv $i ${i/%_png/.png}; done
    for i in */*_txt ; do mv $i ${i/%_txt/.txt}; done
    for i in */*_xml ; do mv $i ${i/%_xml/.xml}; done
    for i in */*_tab ; do mv $i ${i/%_tab/.tab}; done
    for i in */*_RData ; do mv $i ${i/%_RData/.RData}; done
    for i in */*_html ; do mv $i ${i/%_html/.html}; done

    for i in ${prefix}* ; do cd $i ; ../../mk_index_export.sh ; echo $i done ; cd - ; done

    if [ "$GGALLUS" = "TRUE" ] ; then
        for i in ${prefix}* ; do 
	    entry=${i##$prefix}
	    #entry=${entry/pfu_/pfu-vs-}
	    mv $i $entry
        done
    #    for i in P_* ; do echo mv $i ${i/#_/-vs-} done
    #    for i in PC_* ; do echo mv $i ${i/#_/-vs-} done
    #    for i in wt_* ; do echo mv $i ${i/#_/-vs-} done
    else
        for i in ${prefix}* ; do 
	    entry=${i##$prefix}
	    #entry=${entry/_pfu/-vs-pfu}
	    echo mv $i $entry
        done
        #for i in *_wt ; do mv $i ${i/_wt/-vs-wt}
    fi
fi

cp  ../style.css .
cat >index.html<<END
<html>
<head>
  <title>FGSEA results</title>
  <link rel="stylesheet" href="./style.css" />
</head>
</html><body>
<center><h1>FGSEA results</H1></center>
<table>
  <tr><th>Sample</th></tr>
END

for i in sample* ; do
    echo "<tr><td><a href=\"$i/index.html\">${i}</a><br /></td></tr>" ; 
done >> index.html

echo "</body></html>" >> index.html

cd ..

zip -r max_grp_size_500.zip max_grp_size_500

cp max_grp_size_500.zip /home/jr/

exit

if [ "$OLDWAY" = "TRUE" ] ; then
    ifile=$1

    echo '#---------------------------------------------------------------'
    echo "#$ifile"

    fnam=`basename "$ifile"`
    dnam=`dirname "$ifile"`

    echo "#$dnam / $fnam"

    ext=${fnam##*.}
    nam=${fnam%.*}

    newname=${nam//./_}
    echo "#$nam => $newname . $ext"

    if [ "$nam" != "$newname" ] ; then
        echo "mv $ifile $dnam/$newname.$ext" 
    else
        echo "#$ifile unchanged"
    fi
fi
