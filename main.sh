


echo "On sample: $sample"
    
cutadapt -g CCAGCAGCYGCGGTAAN\
    -G  CCGTCAATTCNTTTRAGT \
    -o ${sample}R1_trimmed.fq.gz -p ${sample}R2_trimmed.fq.gz \
    ${sample}R1.fastq.gz ${sample}R2.fastq.gz \
    >> cutadapt_primer_trimming_stats.txt 2>&1

done
