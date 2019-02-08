# PopulationReferenceAlignment
Effects of reference genome on alignment across the genome

Align 10 males and 10 females to GRCh38
Align 10 females to Yoruban reference genome

## Data
* 10 females sample names: A2, A3, A11, A13, A16, A17, A27, A29, A30, A31
* 10 males sample names: A10, A100, A18, A21, A22, A23, A24, A25, A32, A34

## Download the Yoruban assembly:

* ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/709/685/GCA_003709685.1_NA19240_EEE_SV-Pop.1

* We want    *_genomic.fna.gz file
       FASTA format of the genomic sequence(s) in the assembly. Repetitive 
       sequences in eukaryotes are masked to lower-case (see below).
       The FASTA title is formatted as sequence accession.version plus 
       description. The genomic.fna.gz file includes all top-level sequences in
       the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds,
       unplaced scaffolds, and any alternate loci or patch scaffolds). Scaffolds
       that are part of the chromosomes are not included because they are
       redundant with the chromosome sequences; sequences for these placed 
       scaffolds are provided under the assembly_structure directory.
       
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/709/685/GCA_003709685.1_NA19240_EEE_SV-Pop.1/GCA_003709685.1_NA19240_EEE_SV-Pop.1_genomic.fna.gz

gunzip GCA_003709685.1_NA19240_EEE_SV-Pop.1_genomic.fna.gz
```

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/022/975/GCA_002022975.1_10x.supernova.msNA19240/GCA_002022975.1_10x.supernova.msNA19240_genomic.fna.gz

gunzip GCA_002022975.1_10x.supernova.msNA19240_genomic.fna.gz
```

