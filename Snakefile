#!/usr/bin/env python3
import peppy

###########
## NOTES ##
###########

##enrichr for enrichment of snps in supertranscripts with GO terms
    #https://www.nature.com/articles/s41598-020-70527-8#Sec2
        ##PCA from ^^ paper
            ##https://academic.oup.com/bioinformatics/article/28/24/3326/245844?login=true
##plink to convert to .ped/.map --> pca
    ##and https://www.york.ac.uk/res/dasmahapatra/teaching/MBiol_sequence_analysis/workshop4_2019.html#using_plink_for_pca_of_genotypes 
    #https://speciationgenomics.github.io/pca/
    #https://qcb.ucla.edu/wp-content/uploads/sites/14/2016/03/
##other info
    #https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
    #https://grunwaldlab.github.io/Population_Genetics_in_R/analysis_of_genome.html#genetic-differentiation

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

#containers
trinity_container = 'docker://trinityrnaseq/trinityrnaseq:2.11.0'
gatk_container = 'docker://broadinstitute/gatk:4.1.4.0' # same gatk version as in trinity container
bcftools_container = 'docker://biocontainers/bcftools:v1.9-1-deb_cv1'
plink_container = 'docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'
vcftools_container = 'docker://biocontainers/vcftools:v0.1.16-1-deb_cv1'

#########
# RULES #
#########

rule target:
    input:
        'output/02_filtering/filtered_snps.vcf.gz',
        ##variant stats
        'output/02_filtering/variant_stats.txt',
        'output/02_filtering/filtered_variant_stats.txt',
        ##missing stats on filtered variants
        "output/02_filtering/missing_filtered_snps.out",
        ##Plink --> PCA
        expand('output/04_plink/{vcf_group}/no_ldpruning/filtered_snps_plink_pca.eigenvec', vcf_group=["all_samples", "ruakura", "dunedin"]),
        expand('output/04_plink/{vcf_group}/ld_pruned/filtered_snps_plink_pca.eigenvec', vcf_group=["all_samples", "ruakura", "dunedin"]),
        'output/03_final_vcfs/all_samples_pruned_snps.vcf',
        #stats
        #'output/05_stats/nucl_diversity_1kb.windowed.pi',
        #'output/05_stats/TajD_1kb.Tajima.D',
        #'output/05_stats/neg_vs_pos_fst_1kb.windowed.weir.fst'

#further analysis beyond PCA in:
    #https://www.york.ac.uk/res/dasmahapatra/teaching/MBiol_sequence_analysis/workshop4_2019.html#inspecting_the_vcf_file

##TiTv ratio - should be around 2 transitions/transversions

###############################
## 05 - diversity statistics ##
###############################

rule Fst:
    input:
        vcf = 'output/02_filtering/filtered_snps.vcf.gz',
        neg_pop = 'output/04_plink/pc1_negative_samples.txt',
        pos_pop = 'output/04_plink/pc1_positive_samples.txt'
    output:
        'output/05_stats/neg_vs_pos_fst_1kb.windowed.weir.fst'
    params:
        'output/05_stats/neg_vs_pos_fst_1kb'
    log:
        'output/logs/Fst.log'
    singularity:
        vcftools_container
    shell:
        'vcftools '
        '--gzvcf {input.vcf} '
        '--weir-fst-pop {input.neg_pop} '
        '--weir-fst-pop {input.pos_pop} '
        '--fst-window-size 1000 '
        '--out {params} '
        '2> {log}'

##some of these steps may be worth doing on each pop individually - TajD and Pi?

rule TajimasD:
    input:
        'output/02_filtering/filtered_snps.vcf.gz'
    output:
        'output/05_stats/TajD_1kb.Tajima.D'
    params:
        'output/05_stats/TajD_1kb'
    log:
        'output/logs/TajimasD.log'
    singularity:
        vcftools_container
    shell:
        'vcftools '
        '--gzvcf {input} '
        '--TajimaD 1000 '
        '--out {params} '
        '2> {log}'

rule nucl_div:
    input:
        'output/02_filtering/filtered_snps.vcf.gz'
    output:
        'output/05_stats/nucl_diversity_1kb.windowed.pi'
    params:
        'output/05_stats/nucl_diversity_1kb'
    log:
        'output/logs/nucl_div.log'
    singularity:
        vcftools_container
    shell:
        'vcftools '
        '--gzvcf {input} '
        '--window-pi 1000 '
        '--out {params} '
        '2> {log}'

##############################
## 04 - PCA of ASW variants ##
##############################

# how many of these SNPs are in genes DE in location or other analyses? are variants in DE genes or perhaps not influencing expression?

# prune the vcf down to only variants kept in plink analysis
rule prune_vcf:
    input:
        vcf = 'output/03_final_vcfs/renamed/all_samples.vcf',
        prune = 'output/04_plink/all_samples/ld_pruned/filtered_snps_plink.prune.in'
    output:
        vcf = 'output/03_final_vcfs/all_samples_pruned_snps.vcf'
    log:
        'output/logs/prune_vcf.log'
    singularity:
        bcftools_container
    shell:
        'bcftools view '
        '-i \'ID=@{input.prune}\' '
        '{input.vcf} '
        '> {output.vcf} '
        '2> {log}'

rule bcftools_setID:
    input:
        vcf = 'output/03_final_vcfs/all_samples.vcf.gz'
    output:
        'output/03_final_vcfs/renamed/all_samples.vcf'
    log:
        'output/logs/bcftools_setID.log'
    singularity:
        bcftools_container
    shell:
        'bcftools annotate '
        '--set-id +"%CHROM\:%POS\" '
        '{input.vcf} '
        '-O v -o {output} '
        '&> {log}'

# pca on pruned
rule plink_pca_LD:
    input:
        vcf = 'output/03_final_vcfs/{vcf_group}.vcf.gz',
        pruned = 'output/04_plink/{vcf_group}/ld_pruned/filtered_snps_plink.prune.in'
    output:
        'output/04_plink/{vcf_group}/ld_pruned/filtered_snps_plink_pca.eigenvec'
    params:
        'output/04_plink/{vcf_group}/ld_pruned/filtered_snps_plink_pca'
    log:
        'output/logs/plink_pca_LD_{vcf_group}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--vcf {input.vcf} '
        '--double-id '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--extract {input.pruned} '
        '--make-bed '
        '--pca '
        '--out {params} '
        '&> {log}'

# prune dataset of variants in linkage - PCA relies on independent variables
rule plink_prune_linkage:
    input:
        'output/03_final_vcfs/{vcf_group}.vcf.gz'
    output:
        'output/04_plink/{vcf_group}/ld_pruned/filtered_snps_plink.prune.in'
    params:
        indep = '50 10 0.1',     # 50 kb window, 10 SNPs, r2 < 0.1,
        out = 'output/04_plink/{vcf_group}/ld_pruned/filtered_snps_plink'
    log:
        'output/logs/plink_prune_linkage_{vcf_group}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--vcf {input} '
        '--double-id '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--indep-pairwise {params.indep} '
        '--out {params.out} '
        '&> {log}'

##PCA on all snps - some likely in LD though
rule plink_pca_no_ldpruning:
    input:
        vcf = 'output/03_final_vcfs/{vcf_group}.vcf.gz',
    output:
        'output/04_plink/{vcf_group}/no_ldpruning/filtered_snps_plink_pca.eigenvec'
    params:
        'output/04_plink/{vcf_group}/no_ldpruning/filtered_snps_plink_pca'
    log:
        'output/logs/plink_pca_no_ldpruning_{vcf_group}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--vcf {input.vcf} '
        '--double-id '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--make-bed '
        '--pca '
        '--out {params} '
        '&> {log}'

#####################
## 03 - final VCFs ##
#####################

rule bcftools_split_vcf_location:
    input:
        vcf = 'output/02_filtering/filtered_snps.vcf.gz',
        location_samples = 'data/sample_lists/{location}_samples.txt'
    output:
        location_vcf = 'output/03_final_vcfs/{location}.vcf.gz'
    log:
        "output/logs/bcftools_split_vcf_location/{location}.log"
    singularity:
        bcftools_container
    shell:
        'bcftools view -S {input.location_samples} {input.vcf} > {output.location_vcf} 2> {log}'

###################################
## 02 - GATK merging & filtering ##
###################################
# largely following https://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/
# though also consistent with filtering in ASW GBS

# inspect individuals with lots of missing data
    # for filtered SNPs samples are missing 9-17% data - all below 20% threshold for removing samples
rule missing_stats:
    input:
        'output/02_filtering/filtered_snps.vcf.gz'
    output:
        all_samples_vcf = 'output/03_final_vcfs/all_samples.vcf.gz',
        stats = 'output/02_filtering/missing_filtered_snps.out'
    log:
        'output/logs/missing_stats_filtered_snps.log'
    singularity:
        bcftools_container
    shell:
        'cp {input} {output.all_samples_vcf} || exit 1 ; '
        'bcftools stats -s - {input} | grep -E ^PSC | cut -f3,14 > {output.stats} 2> {log}'

# remove missing genotypes on more than 20% of inidividuals and minor allele freq. less than 0.05 (used 0.02 prev but only = 2 ASW)
rule remove_singletons: # returns 6125 SNPs - removes all M.hyp SNPs
    input:
        'output/02_filtering/snp_variants.vcf.gz'
    output:
        'output/02_filtering/filtered_snps.vcf.gz' # quality filtered merged SNP file
    log:
        'output/logs/remove_singletons.log'
    singularity:
        bcftools_container
    shell:
        'bcftools filter '
        '-e "F_MISSING > 0.2 || MAF <= 0.05" '
        '-O z -o {output} {input} 2> {log}'

# Remove multiallelic SNPs and indels, monomorphic SNPs, and SNPs in the close proximity of indels
# AC==0 removes all sites where no alternative alleles called for any samples - not variants
# AC==AN removes all sites where only alternative allele called - not true variant (assembly error)
# SnpGap removes variants close to indels - harder to call with certainty
# -m2 -M2 -v snps keeps only biallelic SNPs - minimum and maximum no. alleles is 2, indels removed
rule only_snps: # returns 1,192,637 SNPs - 1,173,492 from ASW, and 19,145 from Mh
    input:
        'output/02_filtering/merged_filtered_pass.vcf.gz'
    output:
        'output/02_filtering/snp_variants.vcf.gz'
    log:
        'output/logs/only_snps.log'
    singularity:
        bcftools_container
    shell:
        'bcftools filter -e "AC==0 || AC==AN" --SnpGap 10 {input} | bcftools view -m2 -M2 -v snps -O z -o {output} 2> {log}'

rule filtered_variant_stats:
    input:
        'output/02_filtering/merged_filtered_pass.vcf.gz'
    output:
        'output/02_filtering/filtered_variant_stats.txt'
    log:
        'output/logs/filtered_variant_stats.log'
    singularity:
        bcftools_container
    shell:
        'bcftools query {input} -f "%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n" > {output} 2> {log}'

# then select only variants that passed
rule select_pass_merged: # returns 2,011,777 variants
    input:
        'output/02_filtering/merged_filtered.vcf'
    output:
        'output/02_filtering/merged_filtered_pass.vcf.gz'
    singularity:
        bcftools_container
    shell:
        'bcftools view -f "PASS" {input} -o {output} -Oz'

# filters OUT variants that meet expression
# extra filtering DP as suggested by Trinity wiki, and SOR as suggested by GATK
# those variants with missing values fail
rule filter_merged:
    input:
        ref = 'data/asw_mh_transcriptome/supertranscripts.fasta',
        vcf = 'output/01_variants/merged.vcf.gz',
        index = 'output/01_variants/merged.vcf.gz.tbi'
    output:
        fil_vcf = 'output/02_filtering/merged_filtered.vcf'
    log:
        'output/logs/filter_merged.log'
    singularity:
        gatk_container
    threads:
        10
    shell:
        'gatk VariantFiltration '
        '-R {input.ref} '
        '-V {input.vcf} '
        '-O {output.fil_vcf} '
        '-filter-expression "DP < 10.0" --filter-name "DP-10.0" '
        '-filter "SOR > 3.0" --filter-name "SOR3" '
        '--missing-values-evaluate-as-failing true '
        '--invalidate-previous-filters false '
        '2> {log}'

    # analysis in R shows most GATK reccomendations met
    # could filter SOR>3 and depth>10 again (after merging seem to have some below)
rule variant_stats:
    input:
        'output/01_variants/merged.vcf.gz'
    output:
        'output/02_filtering/variant_stats.txt'
    log:
        'output/logs/variant_stats.log'
    singularity:
        bcftools_container
    shell:
        'bcftools query {input} -f "%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n" > {output} 2> {log}'

####################################
## 01 - variant calling & merging ##
####################################

rule index_merged_vcf:
    input:
        gz = 'output/01_variants/merged.vcf.gz'
    output:
        csi = 'output/01_variants/merged.vcf.gz.tbi'
    log:
        'output/logs/index_merged_vcf.log'
    singularity:
        gatk_container
    shell:
        'gatk IndexFeatureFile -F {input.gz} 2> {log}'

# merge individual sample VCFs
# only keep variants with PASS - passed trinity pipeline filtering
rule merge_vcf: # 4,268,761 variants
    input:
        gz = expand('output/01_variants/{sample}/{sample}.vcf.gz', sample=all_samples),
        csi = expand('output/01_variants/{sample}/{sample}.vcf.gz.csi', sample=all_samples)
    output:
        merged = 'output/01_variants/merged.vcf.gz'
    log:
        'output/logs/merge_vcf.log'
    threads:
        20
    shell:
        'bcftools merge {input.gz} -o {output.merged} -Oz -f "PASS" --threads {threads} 2> {log}'

rule index_vcf:
    input:
        gz = 'output/01_variants/{sample}/{sample}.vcf.gz'
    output:
        csi = 'output/01_variants/{sample}/{sample}.vcf.gz.csi'
    singularity:
        bcftools_container
    shell:
        'bcftools index {input.gz}'

rule compress_vcf:
    input:
        renamed = 'output/01_variants/{sample}/renamed_{sample}.vcf'
    output:
        gz = 'output/01_variants/{sample}/{sample}.vcf.gz'
    log:
        'output/logs/compress_vcf/{sample}.log'
    singularity:
        bcftools_container
    shell:
        'bcftools view -I {input.renamed} -Oz -o {output.gz}  &> {log}'

# rename samples to merge
rule rename_samples:
    input:
        fil_vcf = 'output/01_variants/{sample}/filtered_output.vcf'
    output:
        renamed = 'output/01_variants/{sample}/renamed_{sample}.vcf'
    params:
        name_file = temp('output/01_variants/{sample}_sample_name.txt')
    singularity:
        bcftools_container
    shell:
        'echo {wildcards.sample} > {params.name_file} || exit 1 ; '
        'bcftools reheader {input} -s {params.name_file} -o {output} || exit 1 ; '
        'rm output/01_variants/{wildcards.sample}_sample_name.txt'

# uses GATK haplotype caller
# one sample at a time as joint calling not supported for rnaseq
rule GATK_calling_variants:
    input:
        fasta = 'data/asw_mh_transcriptome/supertranscripts.fasta',
        gtf = 'data/asw_mh_transcriptome/supertranscripts.gtf',
        r1 = 'data/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'data/bbduk_trim/{sample}_r2.fq.gz'
    output:
        'output/01_variants/{sample}/filtered_output.vcf'
    params:
        wd = 'output/01_variants/{sample}'
    log:
        'output/logs/01_variants/{sample}.log'
    singularity:
        trinity_container
    threads:
        20
    shell:
        '/usr/local/bin/trinityrnaseq/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py '
        '--st_fa {input.fasta} '
        '--st_gtf {input.gtf} '
        '-p {input.r1} {input.r2} '
        '-o {params.wd} '
        '-t {threads} '
        '-m 167385611648 '
        '&> {log}'

