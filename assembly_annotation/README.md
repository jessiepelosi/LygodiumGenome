# Assembly 

Oxford Nanopore data was basecalled with super-acuracy mode either in MinKnow during data collection with guppy or with dorado v0.3.3.

```
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.1.0 *.pod5 --emit-fastq > file.fastq
```
