# Methylation Analyses

Genome assemblies and raw bisulfite sequencing reads were downloaded for the following taxa: 

| Taxon                    | Lineage   | Annotation     | Accession (Gametophyte/Sporophyte) |References                             | 
| ------------------------ | --------- | -------------- |----------------------------------- | --------------------------------------|
|<i>Chlamydomonas reinhardtii</i>| Green Algae| v5.6    | SRR1045037-8/SRR2051061-2          | Lopez et al. 2015 |
| <i>Marchantia polymorpha | Liverwort | v3.1           | SRR5314027-38/SRR5314015-20        | Schmid et al. 2018                    |
| <i>Arabidopsis thaliana  | Eudicot   | TAIR10.1       | SRR516180/SRR12737894-5            | Hsieh et al. 2016, Halter et al. 2021 |
| <i>Oryza sativa          | Monocot   | v7.0           | SRR8594859/SRR9637033-6            | Kim et al. 2019, Cui et al. 2021      |
  
Raw reads were trimmed with Trimmomatic v0.39. Most reads were PE, but some were SE and the command was adjusted as needed. 
```
trimmomatic PE [accession number]_1.fastq [accession number]_2.fastq [accession number]_1P.fastq [accession number]_1U.fastq [accession number]_2P.fastq [accession number]_2U.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 MINLEN:35 HEADCROP:5
```
## Bismark   
Trimmed reads were mapped to the reference genome using Bismark v0.22.3. The genome was first prepared: 
```
bismark_genome_preparation  --verbose /path/to/genome/
```
Reads were then mapped to the prepared genome: 
```
bismark -p 4 -N 1 -L 0,-0.4 --genome /path/to/genome/ -1 ${R1} -2 ${R2}
```
PCR duplicates were removed from the bismark output: 
```
deduplicate_bismark ${sample}_1P_bismark_bt2_pe.bam
```
And methylation calls were performed and reports were generated with samtools v1.15: 
```
bismark_methylation_extractor --bedGraph --CX ${sample}.bam
coverage2cytosine -CX -o $[sample}.bis_rep.cov --genome_folder path/to/genome/ ${sample}.bismark.cov
bgzip $[sample}.bis_rep.cov
tabix -C -p vcf ${sample}.bis_rep.cov.gz
```
Chloroplast, mitochondrial, and sex-chromosomes were removed prior to downstream analyses as necessary. 

## Alternative Methylation Pipeline 
Given the low mapping efficiency from Bismark for several samples, we ran an alternative pipeline to call methylation. 

Trimmed reads were mapped to the genome with bwa-meth v0.2.7: 
```
bwameth.py --reference ${REF} -t 12 ${R1} ${R2}| samtools view -b | samtools sort > ${sample}.bam
```
Alingment statistics were calculated with samtools v1.15. 
```
samtools flagstat ${sample}.bam
```
PCR and optical duplicates were marked with removed with Picard in GATK 4.3.0.0.
```
gatk MarkDuplicates -I ${sample}.bam -M ${sample}_deduplication_stats.txt -O ${sample}_deduplicated.bam --REMOVE_DUPLICATES true
```
Methylation data was extracted with MethylDackel v0.6.1 and summarized with a custom python script. 
```
MethylDackel extract -@ 8 --CHG --CHH ${REF} ${sample}_deduplicated.bam -o ${sample}_deduplicated_methylDackel
MethylDackel extract -@ 8 --CHH --CHG ${REF} ${sample}_deduplicated.bam -o ${sample}_deduplicated_methylDackel --cytosine_report
MethylDackel extract -@ 8 --CHG --mergeContext ${REF} ${sample}_deduplicated.bam -o ${sample}_deduplicated_methylDackel_merged
bgzip ${sample}_deduplicated_methylDackel.cytosine_report.txt
tabix -C -p vcf ${sample}_deduplicated_methylDackel.cytosine_report.txt.gz
python summarize_methylation.py -i ${sample}_deduplicated_methylDackel_merged_CpG.bedGraph
python summarize_methylation.py -i ${sample}_deduplicated_methylDackel_merged_CHG.bedGraph
```

## Summarizing and Analysis of Methylation Calls
Deduplicated results were read into methylKit v1.20 in R v4.1.3 (see `DiffMethylation.R`) to call differentially methylated bases between gametophyte and sporophytes. 
Methylation reports were used to visualize methylation across contexts in viewBS v0.1.10 (see `VisualizeMethylation.R`). 

# Methylation from Oxford Nanopore Reads 
Use modified basecalling models with dorado v0.3.3.
```
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.2.0 ${POD_DIR} --modified-bases 5mC_5hmC > ${OUTBAM}
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.2.0 ${POD_DIR} --modified-bases 6mA > ${OUTBAM}
```
Convert bam files to fastq files for mapping and retain modification tags with samtools v0.1.18 then map to genome with minimap2.  
```
samtools fastq -T '*' ${OUTBAM} > ${FASTQ}
minimap2 -ax map-ont -t 64 -y --secondary no ${REF} ${IN} --split-prefix temp | samtools view -b | samtools sort > ${BAM}
samtools flagstat ${BAM}
```
Extract methylation information with modkit v. 
```
${MODKIT} pileup --cpg ${BAM} 5mC_5hmC_CHG.bedMethyl --ref ${REF} --log-filepath 5mC_5hmC_CHG.log --threads 16
${MODKIT} pileup --motif CHG 0 ${BAM} 5mC_5hmC_CHG.bedMethyl --ref ${REF} --log-filepath 5mC_5hmC_CHG.log --threads 16
${MODKIT} pileup --motif CHG 0 ${BAM} 5mC_5hmC_CHH.bedMethyl --ref ${REF} --log-filepath 5mC_5hmC_CHG.log --threads 16
```

Transform methylBed file to ViewBS-compatible file: 
```
sed -i 's/ /\t/g' ${BEDMETHYL} | cut -f 1,2,4,6,11,12,13 | \
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $4, $6, $7}' | \
grep 'm' | sed 's/$/\tCG\tCGN/g' | cut -f 1,2,4,5,6,7,8 > ${BED}
cat 5mC_C*.txt | sort -k1,1 -k2,2n | bgzip > 5mC_allContexts.gz #Do for each modification type
tabix -C -p vcf 5mC_allContexts.gz 
```
