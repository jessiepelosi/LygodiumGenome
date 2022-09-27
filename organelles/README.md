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
