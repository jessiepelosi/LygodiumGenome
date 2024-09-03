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
A secondary approach was implemented using a custom R script (see [EstimateGenomeSize.Rmd](https://github.com/jessiepelosi/LygodiumGenome/blob/main/assembly_annotation/EstimateGenomeSize.Rmd)). 

# Assembly 

Oxford Nanopore data was basecalled with super-acuracy mode either in MinKnow during data collection with guppy or with dorado v0.3.3.

```
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.1.0 *.pod5 --emit-fastq > file.fastq
```
### For Lygodium_microphyllum_v0.1.8: 
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

### For Lygodium_microphyllum_v0.1.9:
We generated another assembly by first correcting the raw SUP-called reads with HERRO, which uses a machine-learning method and model. We modified the processing step to keep reads >5kb, rather than >10kb. 

```
./scripts/preprocess_JAP.sh ${reads} ${prefix} 24 24
seqkit seq -ni <reads> > <read_ids>
./scripts/create_batched_alignments.sh all_reads_filt.fastq.gz all_reads_filt.fastq.ids 16 ./aligns
# NOTE THAT WE HAD TO RUN MINIMAP2 OUTSIDE OF HERRO SCRIPTS:
minimap2 -k8g -cx ava-ont -k25 -w17 -e200 -r150 -m4000 -z200 -t64 --dual=yes all_reads_filt.fastq.gz all_reads_filt.fastq.gz > ONT-minimap.out
herro inference --read-alns ./batch -m ./model_v0.1.pt -b 128 all_reads_filt.fastq.gz all_reads_filt_HERRO.fasta
```

We used the corrected reads as input to hifiasm vXX. 
```
hifiasm -o LM_v0.1.9 -t 32 all_reads_filt_HERRO.fasta
```

The draft assemblies were then scaffolded with HiC data. 
```
juicer.sh -z ./references/LMv0.1.8b.medaka1.pilon1.fasta -s none -p references/chrom_sizes.txt --assembly
./yahs references/LMv0.1.8b.medaka1.pilon1.fasta all_merged_dedup_LYMIv0.1.8c.bam
./yahs/juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp references/LMv0.1.8b.medaka1.pilon1.fasta.fai >out_JBAT.log 2>&1
sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part
java -jar -Xmx96G ./scripts/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic <(echo "assembly 1173944728")
```
Import HiC matrix into Juicebox Assembly Tools and manually edit as necessary (e.g., misjoins) and then finalize the assembly with Juicer. 
```
juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp references/LMv0.1.8b.medaka1.pilon1.fasta
```
This yields our scaffolded assembly! To generate subsets of this assembly based on minimum scaffold length, etc. we used SeqKit. 
```
seqkit seq -m 10000 out_JBAT.assembly > LMv1.1.8_m10000.fasta
seqkit grep -f chrms.txt out_JBAT.assembly > LMv1.1.8_chrms.fasta
seqkit stats *.fasta -a
```

# Annotation

### Repeat Annotation

RepeatModeler v2.0 was used to generate a species-specific repeat library for the assembly. RepeatMasker v4.0.5 was then used to mask repetitive elements using the RepeatModeler library. 

```
BuildDatabase -name <db-name> <assembly fasta>
RepeatModeler -database <db-name> -pa 96 -LTRStruct
RepeatMasker -pa 96 -lib <RepeatModeler library> <assembly fasta> -no_is -norna -gff -a --xsmall
```

### Protein Evidence 
Sixteen plant proteomes were downloaded from Phytozome, NCBI, or other repositories to use for protein evidence. 

| Lineage                 | Species                            | Genome Verion |
| ------------------------| ---------------------------------- |---------------
| Green Algae             |<i>Chlamydomonas reinhardtii </i>   | 5.5           |
| Moss / Bryophytes       |<i> Physcomitrella patens</i>       | 3.3           |
| Liverworts / Bryophytes |<i>Marchantia polymorpha</i>        | 3.1           |
| Hornworts / Bryophytes  |<i>Anthoceros agrestis</i> (Oxford) |               |
| Lycophytes              |<i>Selaginella moellendorfii</i>    |               |
| Ferns                   |<i>Marsilea vestita</i>             | 3             |
| Ferns                   |<i>Azolla filiculaoides</i>         | 1.1           |
| Ferns                   |<i>Salvinia cucculata</i>           | 1.2           |
| Ferns                   |<i>Alsophila spinulosa</i>          | 3.1           |
| Ferns                   |<i>Ceratopteris richardii           | 2.1           |
| Ferns                   |<i>Adiantum capillus-verneris       |               |
| ANA / Angiosperms       |<i>Amborella trichopoda             | 1.0           |
| Monocots / Angiosperms  |<i>Oryza sativa</i>                 | 7.0           |
| Monoctos / Angiosperms  |<i>Zea mays</i>                     | 5.0.55        |
| Eudicots / Angiopserms  |<i>Arabidopsis thaliana</i>         | Araport 11    |
| Eudicots / Angiosperms  |<i>Populus trichocarpa</i>          | 4.1           |

### Transcript Evidence 

We used hisat2 to map RNASeq reads to the repeat-masked genome assembly. 
```
hisat2 -x ${INDEX} -1 "$file"_1P.fq.gz -2 "$file"_2P.fq.gz --threads 24 --rna-strandness RF | samtools view -b | samtools sort -o "$file"_${genome}.bam
```

Bad splice junctions were filtered from these alignments with Portcullis v1.2.3. 
```
portcullis full -b --orientation FR --strandedness firststrand ${genome} ${bam} -t 12
```

Trinity v 2.12 was used to generate <i>de novo</i> and genome-guided transcriptome assemblies: 
```
Trinity --CPU 16 --SS_lib_type RF --output ${out} --seqType fq --left ${leftReads} --right ${rightReads}
Trinity --CPU 16 --SS_lib_type RF --output ${out} --genome_guided_bam ${bam} --genome_guided_max_intron 100000 --max_memory 50G
```

Stringtie v2.1.3b was used to stitch transcripts together given a filtered BAM file. 
```
stringtie --rf portcullis.filtered.bam -o stringtie.gtf -p 12 
```
Mikado v2.3.3 was used to pick high-quality transcripts from multiple transcriptome assemblies (<i>de novo</i> and genome-guided Trinity and Stringtie) to be used as transcript evidence for annotation. 
```
mikado configure --list config.txt --reference ${genome} --mode permissive --scoring plants.yaml --copy-scoring
mikado prepare --json-conf configuration.yaml
TransDecoder.LongOrfs -t mikado_prepared.fasta
TransDecoder.Predict -t mikado_prepared.fasta
mikado serialize --json-conf configuration.yaml --orfs mikado.orfs.gff3
mikado pick --configuration configuration.yaml --subloci_out mikado.subloci.gff3
gffread -w transcripts.fa -g ${genome} mikado.loci.gff3
```

## Gene Prediction and Annotation 

Transcript and protein evidence were fed to MAKER3 v3.01.33, which were used to directly predict genes during the first round using the hard-masked genome assembly. The relevant fields in the `maker_opts.ctl` file that were modified are: 
```
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
```
