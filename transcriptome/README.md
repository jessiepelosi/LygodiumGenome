# Transcriptome analysis of <i>Lygodium microphyllum</i>

Raw transcriptomic reads are available through NCBI under project XXXXXXX. 

INSERT TABLE HERE WITH ACCESSION NUMBERS

Assess read quality with FASTQC. 
```
for file in \*.fq.gz; do fastqc "$file";done
```

Trim reads with trimmomatic. 
```
trimmomatic PE "$file"_1.fq.gz "$file"_2.fq.gz "$file"_1P.fq.gz "$file"_1U.fq.gz "$file"_2P.fq.gz \
"$file"_2U.fq.gz HEADCROP:10 SLIDINGWINDOW:4:20 MINLEN:36 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10
```
Reads were mapped to the reference genome using HISAT2. 
```

```
Raw counts were generated with HTSeq 
```
for file in *.bam; do samtools index "$file";done
htseq-count -f bam -m intersection-nonempty -s reverse -t exon -i geneID -o $file-counts $file $gff > $file-gene_summary.txt
```
