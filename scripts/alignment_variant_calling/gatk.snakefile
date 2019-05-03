# This is a snakefile for setting up GATK
# Data: 10 females and 10 males from Northern Kenya
# Right now the names of these individuals are hard-coded here but it would be better to have a json file for this.

FEMALES = ["A11", "A13", "A16", "A17", "A27", "A29", "A2", "A30", "A31", "A3"]
MALES = ["A10", "A100", "A18", "A21", "A22", "A23", "A24", "A25", "A32", "A34"]
chr_females = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20", "chr21", "chr22", "chrX", "chrM"]
chr_males = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
Ref_GRCh38_Y_HardMasked = "/mnt/storage/SAYRES/XY_Trim_Ref/references/gencode.GRCh38.p7_Ymasked/GRCh38_Ymasked_reference.fa"
Ref_GRCh38_Y_PARsMasked = "/mnt/storage/SAYRES/XY_Trim_Ref/references/gencode.GRCh38.p7_minusYPARs/GRCh38_minusYPARs_reference.fa"
FEMALES_MALES = ["A11", "A13", "A16", "A17", "A27", "A29", "A2", "A30", "A31", "A3", "A10", "A100", "A18", "A21", "A22", "A23", "A24", "A25", "A32","A34"]
chr_females_males = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20", "chr21", "chr22", "chrX", "chrM"]
haploid_chr = ["chrX", 'chrY']
# FEMALES_MALES = ["A11", "A13"]
# chr_females_males = ["chr21", "chr22"]


# Functions
def get_variants_str(samples, path, chr):
    finalString = ''
    for i in samples:
        finalString += '--variant ' + path + '/' + chr + '/' + chr + '.' + i + '.merged.g.vcf.gz' + ' '
    return (finalString)

# Set path to programs
gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"

rule all:
    input:
        expand("genotyped_vcfs/{chrm}.gatk.called.raw.haploid.vcf.gz", chrm=haploid_chr)
    input:
        expand("combined_gvcfs/{chrm}.gatk.combinegvcf_haploid.g.vcf.gz", chrm=haploid_chr)
    input:
        expand("gvcfs/{chrm}/{chrm}.{male_sample}.haploid.merged.g.vcf.gz", chrm=haploid_chr, male_sample=MALES)
    input:
        "genotyped_vcfs_allsites/chrY.gatk.called.raw.allsites.vcf.gz"
    input:
        "genotyped_vcfs/chrY.gatk.called.raw.vcf.gz"
    input:
        expand("genotyped_vcfs/{chrm}.gatk.called.raw.vcf.gz", chrm=chr_females_males)
    input:
        "combined_gvcfs/chrY.gatk.combinegvcf.g.vcf.gz"
    input:
        expand("combined_gvcfs/{chrm}.gatk.combinegvcf.g.vcf.gz", chrm=chr_females_males)
    input:
        expand("gvcfs/{chrm}/{chrm}.{male_sample}.merged.g.vcf.gz", chrm=chr_males, male_sample=MALES)
    input:
        expand("gvcfs/{chrm}/{chrm}.{female_sample}.merged.g.vcf.gz", chrm=chr_females, female_sample=FEMALES)
    input:
        expand("bams/{male_sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam", male_sample=MALES),
        expand("bams/{male_sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam.bai", male_sample=MALES)
    input:
        expand("bams/{female_sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam", female_sample=FEMALES),
        expand("bams/{female_sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam.bai", female_sample=FEMALES)

rule mk_symlink_for_bam_females:
    input:
        female_bams = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{female_sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam",
        female_bais = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{female_sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam.bai"
    output:
        female_bams_symlink = "bams/{female_sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam",
        female_bais_symlink = "bams/{female_sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam.bai"
    shell:
        """
        ln -s {input.female_bams} {output.female_bams_symlink};
        ln -s {input.female_bais} {output.female_bais_symlink}
        """

rule mk_symlink_for_bam_males:
    input:
        male_bams = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam",
        male_bais = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam.bai"
    output:
        male_bams_symlink = "bams/{male_sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam",
        male_bais_symlink = "bams/{male_sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam.bai"
    shell:
        """
        ln -s {input.male_bams} {output.male_bams_symlink};
        ln -s {input.male_bais} {output.male_bais_symlink}
        """

rule gatk_gvcf_females:
	input:
		ref = Ref_GRCh38_Y_HardMasked,
		bam = "bams/{female_sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam",
		bai = "bams/{female_sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam.bai"
	output:
		"gvcfs/{chrm}/{chrm}.{female_sample}.merged.g.vcf.gz"
	params:
		gatk = gatk_path,
		chrm_n = "{chrm}"
	shell:
		"{params.gatk} "
		"HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chrm_n} "
		"--emit-ref-confidence BP_RESOLUTION --output {output}"

rule gatk_gvcf_males:
	input:
		ref = Ref_GRCh38_Y_PARsMasked,
		bam = "bams/{male_sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam",
		bai = "bams/{male_sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam.bai"
	output:
		"gvcfs/{chrm}/{chrm}.{male_sample}.merged.g.vcf.gz"
	params:
		gatk = gatk_path,
		chrm_n = "{chrm}"
	shell:
		"{params.gatk} "
		"HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chrm_n} "
		"--emit-ref-confidence BP_RESOLUTION --output {output}"

rule gatk_gvcf_males_haploid:
	input:
		ref = Ref_GRCh38_Y_PARsMasked,
		bam = "bams/{male_sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam",
		bai = "bams/{male_sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam.bai"
	output:
		"gvcfs/{chrm}/{chrm}.{male_sample}.haploid.merged.g.vcf.gz"
	params:
		gatk = gatk_path,
		chrm_n = "{chrm}"
	shell:
		"{params.gatk} "
		"HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chrm_n} "
		"--emit-ref-confidence BP_RESOLUTION -ploidy 1 --output {output}"

# This rule combine gvcf for A, X, and mtDNA
rule gatk_combinegvcfs_both:
    input:
        ref = Ref_GRCh38_Y_HardMasked, #Is this reference ok to use for A and mtDNA?
        gvcfs = lambda wildcards: expand(
			"gvcfs/{chrm}/{chrm}.{sample}.g.vcf.gz",
			sample=FEMALES_MALES,
			chrm=chr_females_males)
    params:
        gatk = gatk_path,
        chrm_n = "{chrm}"

    output:
        "combined_gvcfs/{chrm}.gatk.combinegvcf.g.vcf.gz"

    run:
        variant_files = []
        for i in input.gvcfs:
	           variant_files.append("--variant " + i)
        variant_files = " ".join(variant_files)
        print(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}"""
        )
        shell(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}""")

rule gatk_combinegvcfs_male: #This is for the Y chromosome
    input:
        ref = Ref_GRCh38_Y_PARsMasked,
        gvcfs = lambda wildcards: expand(
			"gvcfs/chrY/chrY.{sample}.merged.g.vcf.gz",
			sample=MALES)
    params:
        gatk = gatk_path,
        chrm_n = "chrY"

    output:
        "combined_gvcfs/chrY.gatk.combinegvcf.g.vcf.gz"

    run:
        variant_files = []
        for i in input.gvcfs:
	           variant_files.append("--variant " + i)
        variant_files = " ".join(variant_files)
        print(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}"""
        )
        shell(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}""")

rule gatk_combinegvcfs_male_haploid: #This is for the X and Y chromosomes in XY individuals
    input:
        ref = Ref_GRCh38_Y_PARsMasked,
        gvcfs = lambda wildcards: expand(
			"gvcfs/{chrm}/{chrm}.{sample}.haploid.merged.g.vcf.gz",
			sample=MALES,
			chrm=haploid_chr)
    params:
        gatk = gatk_path,
        chrm_n = "{chrm}"

    output:
        "combined_gvcfs/{chrm}.gatk.combinegvcf_haploid.g.vcf.gz"

    run:
        variant_files = []
        for i in input.gvcfs:
	           variant_files.append("--variant " + i)
        variant_files = " ".join(variant_files)
        print(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}"""
        )
        shell(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}""")

rule gatk_genotypegvcf_both:
    input:
        ref = Ref_GRCh38_Y_HardMasked,
        gvcf = "combined_gvcfs/{chrm}.gatk.combinegvcf.g.vcf.gz"
    output:
        "genotyped_vcfs/{chrm}.gatk.called.raw.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx10g" """
        """GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output} """

rule gatk_genotypegvcf_male:
    input:
        ref = Ref_GRCh38_Y_PARsMasked,
        gvcf = "combined_gvcfs/chrY.gatk.combinegvcf.g.vcf.gz"
    output:
        "genotyped_vcfs/chrY.gatk.called.raw.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx10g" """
        """GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output} """

rule gatk_genotypegvcf_male_haploid:
    input:
        ref = Ref_GRCh38_Y_PARsMasked,
        gvcf = "combined_gvcfs/{chrm}.gatk.combinegvcf_haploid.g.vcf.gz"
    output:
        "genotyped_vcfs/{chrm}.gatk.called.raw.haploid.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx10g" """
        """GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output} """

rule gatk_genotypegvcf_allsites_male:
    input:
        ref = Ref_GRCh38_Y_PARsMasked,
        gvcf = "combined_gvcfs/chrY.gatk.combinegvcf.g.vcf.gz"
    output:
        "genotyped_vcfs_allsites/chrY.gatk.called.raw.allsites.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx10g" """
        """GenotypeGVCFs -all-sites -R {input.ref} -V {input.gvcf} -O {output} """


# combine gVCF
# /home/tphung3/softwares/gatk-4.1.0.0/gatk --java-options '-Xmx10g' CombineGVCFs \
# -R /mnt/storage/SAYRES/XY_Trim_Ref/references/gencode.GRCh38.p7_Ymasked/GRCh38_Ymasked_reference.fa \
# --variant /scratch/tphung3/PopulationReferenceAlignment/gvcfs/chr21/chr21.A11.merged.g.vcf.gz \
# --variant /scratch/tphung3/PopulationReferenceAlignment/gvcfs/chr21/chr21.A13.merged.g.vcf.gz \
# --output A11_A13.g.vcf.gz

# rule gatk_GenomicsDBImport:
#     input:
#         cohort_sample_map = "cohort.sample_map"
#     output:
#         "{chrm}_database"
#     params:
#         gatk = gatk_path,
#         chrm_n = "{chrm}"
#     shell:
#         "{params.gatk} --java-options -Xmx8g "
#         "GenomicsDBImport "
#         "--genomicsdb-workspace-path {output} "
#         "-L {params.chrm_n} "
#         "--sample-name-map {input.cohort_sample_map} "
#         "--tmp-dir=/scratch/tphung3/PopulationReferenceAlignment "
#         "--reader-threads 5"
#
# rule genotype_gvcfs:
#     input:
#         ref = Ref_GRCh38_Y_HardMasked,
#         database = "{chrm}_database"
#     output:
#         "{chrm}_joint_allsites.vcf"
#     params:
#         gatk = gatk_path
#     shell:
#         "{params.gatk} --java-options -Xmx8g "
#         "GenotypeGVCFs "
#         "--include-non-variant-sites true "
#         "-R {input.ref} "
#         "-V gendb://{input.database} "
#         "-new-qual "
#         "-O {output} "
