## Organellar Assemblies 

<b>Plastome </b>

Five million 150 PE BGI reads were randomly subsampled for the plastome assembly using seqkit v2.0.0.
```
seqkit sample -s100 DP800002273BR_L01_565_1.fq -n5000000 > sub1.fq
seqkit sample -s100 DP800002273BR_L01_565_2.fq -n5000000 > sub2.fq
```
Novoplasty v4.3.1 was used to assemble the plastome using this subset of reads. Program defaults were used with the exception of: 
- Genome size range of 100-200kb
- Kmer of 33
- Read length of 150, insert size of 300 
```
NOVOplasty4.3.1.pl -c config.txt
```
The plastome was annotated using PGA v1.0 using the default parameters. The annotation was then manually curated and the resulting annotation was visualized with OGDraw v1.3.1 webserver.
```
perl PGA.pl -r plastomes/ -t Lygodium_microphyllum_plastome.fasta 
```

<b>Mitochondria</b>

The mitochondrial genome was assembled using a modified script from (Syme et al. 2021) -- see `assembler.sh`. This pipeline uses a baited iterative approach using both long and short read data and first extracts putative mitochondrial reads from the Nanopore data using minimap2 v2.24 and filters for the longest reads up to 250x coverage with filtlong. As the starting seed, we used several mitochondrial sequences from Lygodium spp. downloaded from GenBank. An initial assembly with flye v2.9.2 is generated and polished with the extracted Nanopore reads using racon v1.4.13. Illumina reads that map to the assembly are extracted from the total dataset, downsampled to at most 250x coverage with rasusa v0.7.0 (Hall 2022), mapped with bwa v0.7.17, and used to polish the assembly twice with pilon v1.24.
