# Genome Size Estimation 

Short read data were used to estimate the genome size of <i>Lygodium microphyllum</i>. We first generated frequency distributions of kmers ranging from 17mers to 61mers with Jellyfish v2.3.0. 
```
jellyfish count -C -o 17mer_reads_10202021.jf -m 17 -s 8G -t 200 *.fq 
jellyfish histo -t 10 17mer_reads_10202021.jf > 17mer_reads.histo
```
A model-based approach was implemented in GenomeScope v2.0 to estimate genome size. 
```
genomescope.R -i 17mer_reads.histo -o genomescope_out -k 17
```
A secondary approach was implemented using a custom R script (see `estimate_genome_size_kmers.R`). 

# Assembly 

Oxford Nanopore data was basecalled with super-acuracy mode either in MinKnow during data collection with guppy or with dorado v0.3.3.

```
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.1.0 *.pod5 --emit-fastq > file.fastq
```

The SUP-called reads were used as input into flye v2.9.2, polished once with medaka v1.8, and pilon v1.24. 
```
flye --nano-hq all_SUP_calls.fastq.gz --out-dir Lygodium_microphyllum_draft -g 5.1g -t 48 -i 2
mini_align -i all_SUP_calls.fastq.gz -r assembly.fasta -P -m -p calls_to_draft -t 60
medaka consensus calls_to_draft.bam all.hdf --model r1041_e82_400bps_sup_v4.2.0 --batch 200 --threads 4
medaka stitch all.hdf assembly.fasta assembly.medaka1.fasta
bwa mem assembly.medaka1.fasta ${forward} ${reverse} -t 40 | samtools view -@ 40 -S -q 15 -b | samtools sort -@ 40 -o assembly.medaka1.fasta.mapped.bam
export _JAVA_OPTIONS="-Xmx400G"
pilon --genome assembly.medaka1.fasta --frags assembly.medaka1.fasta.mapped.bam --output assembly.medaka1.pilon1 --outdir ./ --fix all
```

The draft assembly was then scaffolded with HiC data. 

# Annotation

## Repeat Annotation

## Protein Evidence 

## Transcript Evidence 
Trinity v 2.12 was used to generate <i>de novo</i> and genome-guided transcriptome assemblies: 
```
Trinity --CPU 16 --SS_lib_type RF --output ${out} --seqType fq --left ${leftReads} --right ${rightReads}
Trinity --CPU 16 --SS_lib_type RF --output ${out} --genome_guided_bam ${bam} --genome_guided_max_intron 100000 --max_memory 50G
```

Stringtie 

Mikado 
