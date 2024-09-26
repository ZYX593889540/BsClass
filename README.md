# Introduction

BsClass is a workflow tool built on Snakemake for bacterial data analysis, offering a complete quality control module for processing NGS and third-generation sequencing (Pacbio/ONT) fastq data. It also provides assembly modules that integrate third-generation assembly with NGS data correction, as well as an SNP calling module to automate bacterial strain identification workflows. By using BsClass, you can automate tasks such as quality control, assembly, ANI (Average Nucleotide Identity) and SNP calculation, clustering, and strain identification for unknown sample genomes from fastq data.

# Running BsClass

``python BsClass.py --workflow <workflow> --cores <cores> --configfile </path/to/config.yaml>``

The --workflow option in BsClass includes three workflows: NGS, ONT, and Pacbio. The corresponding workflow can be selected based on data from NGS, ONT, and Pacbio.

# Configfile

Before running BsClass, you should first create a configuration file, and an example of the configuration file should be saved in the "test_data" directory.

The explanation of parameters in the configuration file is as follows:

``samples: The sample list is an identifier for samples used to process different sequencing files.``

``name: A dictionary of sample names, where key is a sample name with specific meaning and value is its corresponding database number or sequence file number.``

``ref_fna: The path of the reference genome file.``

``filtered_genus_name: The genus names of the filtered species. It should include ["Phylum", "Class", "Order", "Family", "Genus", "Species"].``

``nums: For paired-end sequencing data, the '1' and '2' correspond to the two reads in the paired-end data: "*_1.fq.gz" and "*_2.fq.gz", respectively. These represent the forward and reverse reads of the paired-end sequencing data. Pacbio and ONT data do not require this parameter.``

``work_dir: The path to the working directory. All intermediate files and result files will be generated in this directory.``

``fastq_dir: The directory where the original Fastq file is stored, used to read sequencing data from that directory.``

``genome_database: The path of the genome database is used for comparison.``

``snp_model: Model for generating comparison pairs. It includes two parameters: full and half. full: Both A vs. B and B vs. A comparisons will be calculated and treated as distinct values. half: Only half of the comparisons will be calculated. For instance, if A vs. B is calculated, then B vs. A will not be, as they are considered equivalent.``

``snp_num: The number of output files into which the input file will be split. When computing the number of SNPs, the input file will be divided into the specified number of output files to guarantee that the number of rows in each output file is approximately the same.``

``calcu_type: Used to specify the type of variant data used for cluster analysis, there are two parameters available: "snp" and "total". When "snp" is selected, the script performs clustering based on single nucleotide polymorphism (SNP) variation data. In contrast, when "total" is chosen, the script clusters based on all types of variation data.``

``cut_num1/cut_num2: SNP threshold for the strain boundary， “cut_num1” is the high precision threshold, and “cut_num2” is the low precision threshold.``

# Download the database

Before running BsClass, it is necessary to download the database of Kraken2 and save it in the "{work_dir}/output/kraken2/db" path.
