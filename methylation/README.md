## Methylation Analyses

Genome assemblies and raw bisulfite sequencing reads were downloaded for the following taxa: 

| Taxon                 | Lineage   | Annotation |References |
| --------------------- | --------- | -------------- |------|
| <i>Arabidopsis thaliana  | Eudicot   | TAIR10.1 | Hsieh et al. 2016, Halter et al. 2021 |
| <i>Marchantia polymorpha | Liverwort | v3.1 | Schmid et al. 2018 |
| <i>Oryza sativa          | Monocot   | v7.0 | Kim et al. 2019, Cui et al. 2021 |
  
Raw reads were trimmed with Trimmomatic v0.39. Most reads were PE, but some were SE and the command was adjusted as needed. 
```
trimmomatic PE [accession number]_1.fastq [accession number]_2.fastq [accession number]_1P.fastq [accession number]_1U.fastq [accession number]_2P.fastq [accession number]_2U.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 MINLEN:35 HEADCROP:5
```
  
Trimmed reads were mapped to the reference genome using Bismark v0.22.3. The genome was first prepared: 
```
bismark_genome_preparation  --verbose /path/to/genome/
```
Reads were then mapped to the prepared genome: 
```
bismark -p 4 -N 1 -L 0,-0.4 --genome /path/to/genome/ -1 [accession number]_1P.fastq -2 [accession number]_2P.fastq 
```
PCR duplicates were removed from the bismark output: 
```
deduplicate_bismark [accession number]_1P_bismark_bt2_pe.bam
```
And methylation calls were performed and reports were generated with samtools v1.15: 
```
bismark_methylation_extractor --bedGraph --CX [deduplicated_bismark].bam
coverage2cytosine -CX -o [accession].bis_rep.cov --genome_folder path/to/genome/ [accession].bismark.cov
bgzip [accession].bis_rep.cov
tabix -C -p vcf [accesion].bis_rep.cov.gz
```
Chloroplast, mitochondrial, and sex-chromosomes were removed prior to downstream analyses as necessary. 

Deduplicated results were read into methylKit v1.20 in R v4.1.3 (see `DiffMethylation.R`) to call differentially methylated bases between gametophyte and sporophytes. 
Methylation reports were used to visualize methylation across contexts in viewBS v0.1.10 (see `VisualizeMethylation.R`). 
