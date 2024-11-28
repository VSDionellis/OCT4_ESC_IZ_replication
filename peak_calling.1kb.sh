#!/usr/bin/bash
#
#requirements:
#1.bedtools
#2.bigWigAverageOverBed (https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
#3.datamash 
#Adjust sensitivity-specificity in line 37 by modifying the perc value of "datamash perc:95". This is practically the background threshold. Values 1-100. The lower the value the more peaks will be called. It is suggested to lower this value for broad peak calling. 
#To split broader peaks with multiple summits modify perc value in lines 55,78. Values 1-100. Higher values split the broader peaks in multiple summits. For peak calling of proader peaks low values are suggested.
#If peaks are too many and close enough (e.g from potential sequencing gaps or artifacts, or broad peak calling), please merge the detected peaks by bedtools merge -d values in lines 37,60,83, where -d applies for the distance in bases (e.g. -d 10000 stands for merging peaks found in a distance of less than 10kb).
echo
echo "Give me the reference genome coordinates in bed format [chr start end]"
read reference_bed
echo
echo "Give the bigwig file for peak calling"
read bw
echo
echo "Set number of threads"
read threads
echo
echo "Give output name"
read name
echo
echo "STEP1: binning reference genome in 1kb windows"
bedtools makewindows -b $reference_bed -w 1000 | grep -v chrM > reference.1kb.bed
echo
echo "STEP2: split per chromosome and annotate bins"
cut -f1 reference.1kb.bed | sort -V | uniq | xargs -P$threads -I{} sh -c 'grep -w {} reference.1kb.bed | awk -v OFS="\t" '\''{print $1,$2,$3,"bin"NR}'\'' > reference.{}.1kb.bed' -- {}
echo
echo "STEP3: extract bigwig values for every 1kb bin per chromosome"
ls *chr*.1kb.bed | sed "s/.bed//g" | xargs -P$threads -I{} sh -c "bigWigAverageOverBed $bw {}.bed {}.tab" -- {}
echo
echo "STEP4: prepare signal input per chromosome"
ls *chr*.1kb.bed | sed "s/.bed//g" | xargs -P$threads -I{} sh -c 'paste {}.bed {}.tab | cut -f 1-4,9 > {}.signal.bed' -- {}
echo
echo "STEP5: get global first pass peaks per chromosome"
for f in *.1kb.signal.bed
 do echo $(cat $f | datamash perc:95 5) | paste $f - | awk -v OFS="\t" '{if(NF==5)$6=p FS $6; else p=$6 FS $6}1' | cut -d ' ' -f1 | awk '$5 > $6' | bedtools merge -d 10000 -i - | awk -v OFS="\t" '{print $1,$2,$3,$1"_region"NR}' > ${f%.1kb.signal.bed}.first_pass.peaks.bed
done 
echo
echo "STEP6: binning first pass peaks in 1kb windows and annotate"
for f in *.first_pass.peaks.bed
 do bedtools makewindows -b $f -w 1000 | bedtools intersect -a - -b $f -wa -wb | awk -F '\t' 'BEGIN { OFS=FS } $7 != save { counter = 1; save = $7 } { print $0, counter++ }' | cut -f1-3,7,8 | sed 's/\(.*\)\t/\1_/' > ${f%.bed}.1kb.annotated.bed
done
echo
echo "STEP7: extract bigwig values for every 1kb bin of first pass peaks per chromosome"
ls *.first_pass.peaks.1kb.annotated.bed | sed "s/.bed//g" | xargs -P$threads -I{} sh -c "bigWigAverageOverBed $bw {}.bed {}.tab" -- {}
echo
echo "STEP8: prepare signal input of first pass peaks"
for f in *.first_pass.peaks.1kb.annotated.bed
 do paste $f ${f%.bed}.tab | cut -f 1-4,9 > ${f%.bed}.signal.bed
done
echo
echo "STEP9: calculate medians or desired percentile of first pass peaks to break broad peaks and define multiple summits"
for f in *.first_pass.peaks.1kb.annotated.signal.bed
 do sed 's/\(.*\)_/\1\t/' $f | cut -f4,6 | datamash -W groupby 1 perc:50 2 > ${f%.bed}.percentile.txt
done
echo
echo "STEP10: get local second peaks per chromosome"
for f in *.first_pass.peaks.1kb.annotated.signal.bed
 do sed "s/\(.*\)_/\1\t/" $f | cut -f1-4,6 | join -1 4 -2 1 - ${f%.bed}.percentile.txt | awk '$5 > $6' | cut -d ' ' -f 2-4 | sed "s/ /\t/g" | bedtools merge -d 10000 -i - | awk -v OFS="\t" '{print $1,$2,$3,$1"_region"NR}' > ${f%.first_pass.peaks.1kb.annotated.signal.bed}.second_pass.peaks.bed
done
echo
echo "STEP11: binning second pass peaks in 1kb windows and annotate"
for f in *.second_pass.peaks.bed
 do bedtools makewindows -b $f -w 1000 | bedtools intersect -a - -b $f -wa -wb | awk -F '\t' 'BEGIN { OFS=FS } $7 != save { counter = 1; save = $7 } { print $0, counter++ }' | cut -f1-3,7,8 | sed 's/\(.*\)\t/\1_/' > ${f%.bed}.1kb.annotated.bed
done
echo
echo "STEP12: extract bigwig values for every 1kb bin of second pass peaks per chromosome"
ls *.second_pass.peaks.1kb.annotated.bed | sed "s/.bed//g" | xargs -P$threads -I{} sh -c "bigWigAverageOverBed $bw {}.bed {}.tab" -- {}
echo
echo "STEP13: prepare signal input of second pass peaks"
for f in *.second_pass.peaks.1kb.annotated.bed
 do paste $f ${f%.bed}.tab | cut -f 1-4,9 > ${f%.bed}.signal.bed
done
echo
echo "STEP14: calculate medians or desired percentile of second pass peaks to break broad peaks and define multiple summits"
for f in *.second_pass.peaks.1kb.annotated.signal.bed
 do sed 's/\(.*\)_/\1\t/' $f | cut -f4,6 | datamash -W groupby 1 perc:50 2 > ${f%.bed}.percentile.txt
done
echo
echo "STEP15: get local third pass peaks per chromosome"
for f in *.second_pass.peaks.1kb.annotated.signal.bed
 do sed "s/\(.*\)_/\1\t/" $f | cut -f1-4,6 | join -1 4 -2 1 - ${f%.bed}.percentile.txt | awk '$5 > $6' | cut -d ' ' -f 2-4 | sed "s/ /\t/g" | bedtools merge -d 10000 -i - | awk -v OFS="\t" '{print $1,$2,$3,$1"_region"NR}' > ${f%.second_pass.peaks.1kb.annotated.signal.bed}.third_pass.peaks.bed
done
echo
echo "STEP16: merge peaks all chromosomes in one file"
cat *third_pass.peaks.bed | sort -V -k1,1 -k2,2 -k3,3 > $name.1kb_resolution.peaks.bed
echo
echo "STEP17: remove files"
rm *.1kb.bed *tab *.1kb.signal.bed *pass.peaks.bed *pass.peaks.1kb.annotated.bed *pass.peaks.1kb.annotated.signal.bed *pass.peaks.1kb.annotated.signal.percentile.txt
echo
echo "DONE!!!"
echo
exit











