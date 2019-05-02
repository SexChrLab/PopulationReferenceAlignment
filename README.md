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

## Variant calling with GATK

* Download GATK version 4.1.0.0
```
wget https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip

unzip gatk-4.1.0.0.zip
```

## Temp analysis (This is for internal purposes)
### 1. Examine statistics from VCF file for chr21, chrX, chrY, and mtDNA

```
 python ~/softwares/tanya_repos/vcfhelper/extract_stats_from_vcf.py QD FS SOR MQ MQRankSum ReadPosRankSum --vcf ../genotyped_vcfs/chr21.gatk.called.raw.vcf.gz --outfile chr21_prefiltering_annotations.txt
 
 python ~/softwares/tanya_repos/vcfhelper/extract_stats_from_vcf.py QD FS SOR MQ MQRankSum ReadPosRankSum --vcf ../genotyped_vcfs/chrX.gatk.called.raw.vcf.gz --outfile chrX_prefiltering_annotations.txt
 
 python ~/softwares/tanya_repos/vcfhelper/extract_stats_from_vcf.py QD FS SOR MQ MQRankSum ReadPosRankSum --vcf ../genotyped_vcfs/chrM.gatk.called.raw.vcf.gz --outfile chrM_prefiltering_annotations.txt
 
 python ~/softwares/tanya_repos/vcfhelper/extract_stats_from_vcf.py QD FS SOR MQ MQRankSum ReadPosRankSum --vcf ../genotyped_vcfs/chrY.gatk.called.raw.vcf.gz --outfile chrY_prefiltering_annotations.txt
```

### 2. Obtain QD for variants that are heterozygous and variants that are homozygous for the alternate allele

```
python find_het_homalt.py genotyped_vcfs/chrY.gatk.called.raw.vcf.gz temp_analyze_vcf_stats/chrY_hets_pos_QD_MQ_DP.txt temp_analyze_vcf_stats/chrY_homoalt_pos_QD_MQ_DP.txt
```

### 3. Obtain annotations for all of the sites

```
python extract_annotation.py genotyped_vcfs/chrY.gatk.called.raw.vcf.gz QD temp_analyze_vcf_stats/chrY_all_QD.txt
python extract_annotation.py genotyped_vcfs/chrY.gatk.called.raw.vcf.gz DP temp_analyze_vcf_stats/chrY_all_DP.txt
python extract_annotation.py genotyped_vcfs/chrY.gatk.called.raw.vcf.gz MQ temp_analyze_vcf_stats/chrY_all_MQ.txt
```

### 4. Examine annotations for chrX for XY individuals

1. Subset the VCF file to include XY individuals
```
bcftools view -s A10_whole_genome,A100_whole_genome,A18_whole_genome,A21_whole_genome,A22_whole_genome,A23_whole_genome,A24_whole_genome,A25_whole_genome,A32_whole_genome,A34_whole_genome genotyped_vcfs/chrX.gatk.called.raw.vcf.gz > temp_analyze_vcf_stats/chrX.gatk.called.raw.males.vcf.gz
```

2. Extract and plot QD, DP, MQ for all of the variants

```
python extract_annotation.py temp_analyze_vcf_stats/chrX.gatk.called.raw.males.vcf.gz QD temp_analyze_vcf_stats/chrX_males_all_QD.txt
python extract_annotation.py temp_analyze_vcf_stats/chrX.gatk.called.raw.males.vcf.gz DP temp_analyze_vcf_stats/chrX_males_all_DP.txt
python extract_annotation.py temp_analyze_vcf_stats/chrX.gatk.called.raw.males.vcf.gz MQ temp_analyze_vcf_stats/chrX_males_all_MQ.txt
```

3. Extract and plot QD, DP, MQ for variants that are heterozygous and variants that are homozygous for the alternate allele

```
python find_het_homoalt_extract_annotation.py temp_analyze_vcf_stats/chrX.gatk.called.raw.males.vcf.gz QD temp_analyze_vcf_stats/chrX_males_hets_QD.txt temp_analyze_vcf_stats/chrX_males_homo_QD.txt

python find_het_homoalt_extract_annotation.py temp_analyze_vcf_stats/chrX.gatk.called.raw.males.vcf.gz DP temp_analyze_vcf_stats/chrX_males_hets_DP.txt temp_analyze_vcf_stats/chrX_males_homo_DP.txt

python find_het_homoalt_extract_annotation.py temp_analyze_vcf_stats/chrX.gatk.called.raw.males.vcf.gz MQ temp_analyze_vcf_stats/chrX_males_hets_MQ.txt temp_analyze_vcf_stats/chrX_males_homo_MQ.txt
```

### 4. Examine annotations for chrY for XY individuals

```
for anno in QD DP MQ; do
python extract_annotation.py genotyped_vcfs/chrY.gatk.called.raw.vcf.gz ${anno} temp_analyze_vcf_stats/chrY_all_${anno}.txt
done;
```

```
for anno in QD DP MQ; do
python find_het_homoalt_extract_annotation.py genotyped_vcfs/chrY.gatk.called.raw.vcf.gz ${anno} temp_analyze_vcf_stats/chrY_hets_${anno}.txt temp_analyze_vcf_stats/chrY_homo_${anno}.txt
done;
```

### 5. Compare between hard-filter and VQSR (for chr22)
#### a. Concordance analyses

* Use the array VCF generated by Emma Howell

#### b. VQSR for chr22 (pilot)
* Download data for VQSR from https://software.broadinstitute.org/gatk/download/bundle

#### c. 
* For the array VCF file, subset for 20 individuals and remove sites that are homozygous reference

```
bcftools view -s A2,A3,A11,A13,A16,A17,A27,A29,A30,A31,A10,A100,A18,A21,A22,A23,A24,A25,A32,A34 chr22_array.vcf > chr22_array_20.individuals.vcf
python remove_homozygous_reference.py ../data/chr22_array_20.individuals.vcf ../data/chr22_array_20.individuals_rmhomoref.vcf
```

* Obtain sites from VCF for Venn diagram:

```
python extract_snp_from_vcf.py ../data/chr22_array_20.individuals_rmhomoref.vcf ../results/chr22_array_snps.txt
python extract_snp_from_vcf.py ../data/chr22.gatk.called.raw_fixheader_gatk.rec.hardfilter_selectvariant.vcf ../results/chr22_gatk.hard.filter_snps.txt
python extract_snp_from_vcf.py ../results/vqsr/chr22.gatk.called.raw_vqsr_sv.vcf ../results/chr22_vqsr_snps.txt
python extract_snp_from_vcf.py ../data/chr22.gatk.called.raw.vcf ../results/chr22_before.filtering_snps.txt
```

- Then, plot Venn diagram (PopulationReferenceAlignment/compare_hardfilter_vqsr/scripts/plot_venn.R)

* Obtain the number of variants where the genotypes are the same between array and wholegenome:

```
bcftools view -s A2,A3,A11,A13,A16,A17,A27,A29,A30,A31,A10,A100,A18,A21,A22,A23,A24,A25,A32,A34 chr22.gatk.called.raw_fixheader.vcf > chr22.gatk.called.raw_fixheader_fixorder.vcf
python check_genotype_concordance.py
```

* SFS

```
python ~/softwares/tanya_repos/popgen_tools/popgen_tools.py --vcf_file ../../data/chr22_array_20.individuals.vcf --sfs_all --sfs_all_out chr22_array_20.individuals_sfs.out
python ~/softwares/tanya_repos/popgen_tools/popgen_tools.py --vcf_file ../../data/chr22.gatk.called.raw.vcf --sfs_all --sfs_all_out chr22_before.filtering_sfs.out
python ~/softwares/tanya_repos/popgen_tools/popgen_tools.py --vcf_file ../../data/chr22.gatk.called.raw_fixheader_gatk.rec.hardfilter_selectvariant.vcf --sfs_all --sfs_all_out chr22_hard.filtering_sfs.out
python ~/softwares/tanya_repos/popgen_tools/popgen_tools.py --vcf_file ../vqsr/chr22.gatk.called.raw_vqsr_sv.vcf --sfs_all --sfs_all_out chr22_vqsr_sfs.out
```
 - Plot (PopulationReferenceAlignment/compare_hardfilter_vqsr/scripts/compare_sfs.R)

### 6. Restricting whole genome VCF to contain sites that are used in the array

#### Identify sites in the array
```
sed 's/chr22_/chr22/g' chr22_array_positions_GRCh38_fmtchr.bed > chr22_array_positions_GRCh38_fmtchrf_ix.bed
sed -i 's/chr22KI270879v1_alt/chr22/g' chr22_array_positions_GRCh38_fmtchr_fix.bed
```

#### Subset whole genome VCF to contain only sites that are in the array
```
./select_variants_in_array.sh
```

#### Subset array-VCF and wholegenome-VCF per individual
```
./subset_vcf_per_individual.sh
```

#### Concordance per individual

* Use python script `check_genotype_concordance_per_individual.py`. See wrapper script `wrapper_check.genotype.concordance.per.individual.sh`.

#### SFS
* Goal is to compare singleton when restricting to sites only in array vs not

```
python ~/softwares/tanya_repos/popgen_tools/popgen_tools.py --vcf_file ../../data/chr22.gatk.called.raw_fixheader_array.vcf --sfs_all --sfs_all_out chr22_before.filtering_array_sfs.out 

python ~/softwares/tanya_repos/popgen_tools/popgen_tools.py --vcf_file ../../data/chr22.gatk.called.raw_fixheader_gatk.rec.hardfilter_array.vcf --sfs_all --sfs_all_out chr22_hard.filtering_array_sfs.out

 python ~/softwares/tanya_repos/popgen_tools/popgen_tools.py --vcf_file ../../data/chr22.gatk.called.raw_vqsr_sv_fixheader_array.vcf --sfs_all --sfs_all_out chr22_vqsr_array_sfs.out
```

#### Concordance analysis
* Date: 05/02 (Use the VCF file for array that has been fixed by Emma)
* Since the VCF file for array by Emma only contain 10 female individual, I will subset my whole genome VCF file for these individuals. 

```
bcftools view -s A2,A3,A11,A13,A16,A17,A27,A29,A30,A31 chr22.gatk.called.raw_fixheader_array.vcf > chr22.gatk.called.raw_fixheader_array_females.vcf
```

* Run the script `snparray_wholegenome_concordance.sh`.

### 7. For chrX and chrY in XY individuals that were genotyped using the haploid option, calculate and plot read balance
* We want to do this per individual
* Step 1: Use bcftools to subset per individual
