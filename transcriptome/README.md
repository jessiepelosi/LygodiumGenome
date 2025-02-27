# Transcriptome Analyses 

Raw reads were trimmed with Trimmomatic: 
```
trimmomatic PE "$file"_1.fq.gz "$file"_2.fq.gz "$file"_1P.fq.gz "$file"_1U.fq.gz "$file"_2P.fq.gz "$file"_2U.fq.gz HEADCROP:10 ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35
```

Hisat2 v20230203 (Kim et al. 2019) was used to map RNASeq reads to the genome and read counts generated with HTSeq v2.0.3 (Putri et al. 2022). 
```
hisat2 -x ${INDEX} -1 trimmed_reads/"$file"/"$file"_1P.fq.gz -2 trimmed_reads/"$file"/"$file"_2P.fq.gz --threads 24 --rna-strandness RF | samtools view -b | samtools sort -o "$file"_2v1.1.9.bam
htseq-count -f bam -m union -r pos -s reverse -t mRNA -i ID -o $file-stranded-exon-counts $file $gff > $file-stranded-exon_summary.txt
```

Read counts were used in DESeq2 v1.44.0 (Love et al.2014) for several different R scripts (in this directory) depending on the analysis. 

For alternative splicing, we used the custom scripts in Chamala et al. (2015) and Mei et al. (2017) after running PASA. 

