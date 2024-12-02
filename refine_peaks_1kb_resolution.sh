#!/usr/bin/bash
#
#requirements:
#1.bedtools
#2.bigWigAverageOverBed (https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
#
echo
echo "Give me the called peaks in bed format"
read peaks
echo
echo "Give the bigwig file of peak calling"
read bw
echo
echo "Set number of threads"
read threads
echo 
bedtools makewindows -b $peaks -w 1000 | bedtools intersect -a - -b $peaks -wa -wb  | awk -F '\t' 'BEGIN { OFS=FS } $7 != save { counter = 1; save = $7 } { print $0, counter++ }' | cut -f1-3,7,8 | sed 's/\(.*\)\t/\1_/' > ${peaks%.bed}.1kb_annotated.bed
echo
echo
ls *1kb_annotated.bed | sed "s/.bed//g" | xargs -P$threads -I{} sh -c "bigWigAverageOverBed $bw {}.bed {}.tab" -- {}
echo
echo
for f in *.1kb_annotated.tab
 do paste ${f%.1kb_annotated.tab}.1kb_annotated.bed $f | cut -f1-4,9 | sed "s/_/\t/g" | awk -v OFS="\t" '$6>max[$4]{max[$4]=$6; row[$4]=$0} END{for (i in row) print row[i]}' > ${f%.1kb_annotated.tab}.summits.bed
done
echo
rm *.1kb_annotated.tab *.1kb_annotated.bed
echo "DONE!!!"
echo
exit
