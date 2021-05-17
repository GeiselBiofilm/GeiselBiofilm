# Author: Owen Wilkins (Data Analytics Core)



conda activate geiselbiof

bash main.sh \
	samples.txt \
	CCAGCAGCYGCGGTAAN \
	CCGTCAATTCNTTTRAGT 



mkdir trim 
mkdir dada2
mkdir phyloseq

samples_in=$1
samples_in="samples.txt"
#R1_five_prime_adapter="$2"
#R2_five_prime_adapter="$3"
R1_five_prime_adapter="CCAGCAGCYGCGGTAAN"
R2_five_prime_adapter="CCGTCAATTCNTTTRAGT"

samples=`cut -f1 $samples_in | tail -n +11`

# run cutadapt for all samples 
for sample in $samples; do
cutadapt \
    -g $R1_five_prime_adapter \
    -G $R2_five_prime_adapter \
    -o trim/${sample}.R1.trim.fastq.gz \
    -p trim/${sample}.R2.trim.fastq.gz \
    -m 1 \
    -j 4 \
    raw-data/${sample}_R1.fastq.gz raw-data/${sample}_R2.fastq.gz \
    > trim/${sample}_cutadapt.log
done





Rscript pipeline.R 




Rscript phyloseq.R 

# add fastqc & multiqc 
# write conditional for if pipeline succeeded (all output files exist)

