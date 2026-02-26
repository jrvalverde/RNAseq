echo -e "file\tsample\tstatus\ttype\treplica\tpersistent\trecovered\twt\tresistant\tsensible\tko" 
ls -1 *R1*.bam \
| while read line ; do 
    file=${line%.*}
    sample=${line%_S*}
    status=${sample##*_}
    type=${sample%_*}
    replica=${file#*_S}
    replica=${replica%%_*}
    if echo "$file" | grep -q P_ ; then persistent=Y ; else persistent=N ; fi
    if echo "$file" | grep -q PC_ ; then recovered=Y ; else recovered=N ; fi
    if echo "$file" | grep -q DF-1_ ; then wt=Y ; else wt=N ; fi
    if [ "$persistent" = 'Y' ] ; then resistant=Y ; else resistant=N ; fi
    if [ "$persistent" = 'Y' ] ; then sensible=N ; else sensible=Y ; fi
    if echo "$file" | grep -q KO ; then ko=Y ; else ko=N ; fi
    echo -e "${file}\t${sample}\t${status}\t${type}\t${replica}\t${persistent}\t${recovered}\t${wt}\t${resistant}\t${sensible}\t${ko}" 
done
