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
echo 
bedtools makewindows -b $peaks -w 1000 | bedtools intersect -a - -b $peaks -wa -wb  | awk -F '\t' 'BEGIN { OFS=FS } $7 != save { counter = 1; save = $7 } { print $0, counter++ }' | cut -f1-3,7,8 | sed 's/\(.*\)\t/\1_/' > $peaks.1kb_annotated.bed
echo
echo
for f in *.1kb_annotated.bed
 do bigWigAverageOverBed $bw $f ${f%.bed.1kb_annotated.bed}.1kb_annotated.tab
done
echo
echo
for f in *.1kb_annotated.tab
 do paste ${f%.1kb_annotated.tab}.bed.1kb_annotated.bed $f | cut -f1-4,9 | sed "s/_/\t/g;s/\tregion/_region/g" | awk -v OFS="\t" '{print $4, $6, $0}' | sort -k1,1 -k2,2nr | awk -v OFS="\t" '!seen[$1]++ {print $3,$4,$5,$6,$8}' > ${f%.1kb_annotated.tab}.summits.bed
done

