#!/bin/bash
#
echo
echo
echo "EdU or EU?"
read run
echo
echo
echo "give me the reference genome"
read reference
echo
echo
echo "how many threads would you like to use?"
read threads
#
#
#
#################################### multi mv fastq.gz files ###################################################################################
#
#
#
find *.fastq.gz -type f -print0 | xargs -0 -n 1 | sed "s/_L[0-9].*//g" | xargs -n 1 -P 2 -I {} sh -c "mv {}*.fastq.gz {}.fastq.gz" -- {}
wait
#
#
#
#################################### unzip fastq.gz in parallel: argument -P specifies numbers of cores ########################################
#
#
#
find *.fastq.gz -type f -print0 | sed "s/.fastq.gz//g" | xargs -0 -n 1 -P 2 -I {} sh -c " gunzip -c {}.fastq.gz > {}.fastq" -- {}
wait
#
#
#
#################################### prepare sai files for alignment ###########################################################################
#
#
#
for f in *.fastq
do 
    bwa aln -t $threads $reference $f > ${f%.fastq}.sai
done
wait
#
#
#
#################################### align sai files ###########################################################################################
#
#
#
find *.sai -type f -print0 | sed "s/.sai//g" | xargs -0 -n 1 -P $threads -I {} sh -c "bwa samse -n1 $reference {}.sai {}.fastq > {}.sam" -- {}
wait
#
#
#
#################################### remove unnecessary fastq files ############################################################################
#
#
#
rm *fastq
#
#
#
#################################### run calculate and addSD scripts ###########################################################################
#
#
#
genome=$(echo $reference | sed "s/\//\t/g" | awk '{print $NF}')
#
#
#
if [ "$genome" == "allchr_masked.fasta" ] || [ "$genome" == "allchr.fasta" ] && [ "$run" == "EdU" ]
then
    find *.sam -type f -print0 | xargs -0 -n 1 -P $threads -I {} sh -c "./calculate_number_hits_1files_qual37_adjust-with-noEdU_checkchr_quality_clean_multichr_v1_0b.pl 10000   37  {}" -- {}  &&  find *_bin-size_10000_quality_37_chr1-X_adjbin_0b.csv -type f -print0 | xargs -0 -n 1 -P $threads -I {} sh -c "./add_SD_Back_multichr_v8relB3_0b.pl {}" -- {}
#
#
#
elif [ "$genome" == "allmchr_masked.fasta" ] || [ "$genome" == "allmchr.fasta" ] && [ "$run" == "EdU" ]
then
    find *.sam -type f -print0 | xargs -0 -n 1 -P $threads -I {} sh -c "./calculate_number_hits_1files_qual37_adjust-with-noEdU_checkchr_quality_clean_multichr_v1_0b_mouse.pl 10000   37  {}" -- {}  &&  find *_bin-size_10000_quality_37_mchr1-XY_adjbin_0b.csv -type f -print0 | xargs -0 -n 1 -P $threads -I {} sh -c "./add_SD_Back_multichr_v8relB_0b_mouse.pl {}" -- {}
#
#
#
#elif [ "$genome" == "hg19.fasta" ] && [ "$run" == "EU" ]
#then
#    find *.sam -type f -print0 | xargs -0 -n 1 -P $threads -I {} sh -c "./?????????????????????????? 10000   37  {}" -- {}  &&  &&  find *noHU_bin-size_10000_quality_37_mchr1-XY_adjbin_0b.csv -type f -print0 | xargs -0 -n 1 -P $threads -I {} sh -c "./add_SD_Back_multichr_v8relB_0b_mouse.pl {}" -- {}
#
#
#
elif [ "$genome" == "allmchr_masked.fasta" ] || [ "$genome" == "allmchr.fasta" ] && [ "$run" == "EU" ]
then
    find *.sam -type f -print0 | xargs -0 -n 1 -P $threads -I {} sh -c "./calculate_number_hits_1files_qual37_adjust-with-noEdU_checkchr_quality_clean_multichr_EU_v1_0b_mouse.pl 10000   37  {}" -- {}  &&  find *EU*_bin-size_10000_quality_37_mchr1-X_adjbin_0b.csv -type f -print0 | xargs -0 -n 1 -P $threads -I {} sh -c "./add_SD_Back_multichr_v8relB_0b_mouse.pl {}" -- {}
#
#
#
else
    echo "no correspondence between reference and scripts"
fi
#
#
#
exit














