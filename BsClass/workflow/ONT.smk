import os
import sys
import argparse

BASEDIR = os.path.dirname(os.path.abspath(__file__))

default_configfile = config["default_configfile"]
parser = argparse.ArgumentParser()
parser.add_argument("--configfile", default=default_configfile, help="Path to the configuration file")
args, unknown = parser.parse_known_args()

configfile: args.configfile

samples = config["samples"]
work_dir = config["work_dir"]
fastq_dir = config["fastq_dir"]
name = config["name"]
ref_fna = config["ref_fna"]
name = {str(k): str(v) for k, v in config["name"].items()}  # 将name中的键和值都转换为字符串
ref_fna = {str(k): str(v) for k, v in config["ref_fna"].items()}  # 将ref_fna中的键和值都转换为字符串
filter_name = config["filtered_genus_name"]
genome_database = config["genome_database"]
snp_model = config["snp_model"]
snpfile_num = config["snp_num"]
calcu_type = config["calcu_type"]
cut_num1 = config["cut_num1"]
cut_num2 = config["cut_num2"]


rule all:
    input:
        expand("{work_path}/output/truth_bed/{sample}/truth_{sample}.bed", sample=samples, work_path=[work_dir]),
        expand("{work_path}/output/multiqc/genome_judge/sum.html",work_path=[work_dir]),
        expand("{work_path}/output/unicycler/{sample}/{sample}.fasta", work_path=work_dir, sample=samples),
        expand("{work_path}/output/skani/sample_db_matrix.txt", work_path=work_dir),
        expand("{work_path}/output/skani/snp_mess.txt", work_path=work_dir),
        expand("{work_path}/output/skani/snp_mess_part_{i}.txt", i=range(1, snpfile_num + 1), work_path=work_dir),
        expand("{work_path}/output/snippy/job_{i}.txt", i=range(1, snpfile_num + 1), work_path=work_dir),
        expand("{work_path}/output/snp_group_low.txt", work_path=work_dir),
        expand("{work_path}/output/non_redundant_high.txt", work_path=work_dir)



rule long_trim:
    input:
        long_raw=lambda wildcards: "{fastq_path}/{sample}.fastq".format(sample=wildcards.sample, fastq_path=fastq_dir)
    output:
        long_trim="{work_path}/output/filtlong/{sample}.fastq.gz"
    threads:1
    conda:
        "envs/ont.yaml"
    wildcard_constraints:
        sample="|".join(samples)
    shell:
        """
        filtlong --min_length 1000 --keep_percent 95 {input.long_raw} | gzip > {output.long_trim}
        """


rule kraken2:
    input:
        long_trim=lambda wildcards: "{work_path}/output/filtlong/{sample}.fastq.gz".format(sample=wildcards.sample, work_path=work_dir)
    output:
        report="{work_path}/output/kraken2_long/{sample}.report",
        out="{work_path}/output/kraken2_long/{sample}.output",
        readlist="{work_path}/output/kraken2_long/{sample}.readlist"
    params:
         kraken_db=lambda wildcards: "{work_path}/output/kraken2/db/minikraken2_v1_8GB".format(work_path=work_dir),
         filter_name=lambda wildcards: filter_name[wildcards.sample]
    threads:6
    wildcard_constraints:
        sample="|".join(samples)
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        kraken2 --threads {threads} --db {params.kraken_db} \
         --use-names --report {output.report} \
         {input.long_trim} > {output.out}
        fgrep {params.filter_name} {output.out} | cut -f2 > {output.readlist}
        """


rule flye:
    input:
        readlist=lambda wildcards: "{work_path}/output/kraken2_long/{sample}.readlist".format(work_path=work_dir, sample=wildcards.sample),
        long_trim=lambda wildcards: "{work_path}/output/filtlong/{sample}.fastq.gz".format(sample=wildcards.sample, work_path=work_dir)
    output:
        long_filter="{work_path}/output/clean_long/{sample}.fastq",
        out_genome="{work_path}/output/flye/{sample}/assembly.fasta"
    params:
        flye_out="{work_path}/output/flye/{sample}"
    threads:8
    wildcard_constraints:
        sample="|".join(samples)
    conda:
        "envs/ont.yaml"
    shell:
        """
        seqtk subseq {input.long_trim} {input.readlist} > {output.long_filter}
        flye -o {params.flye_out} --plasmids --threads {threads} \
         --nano-raw {output.long_filter} --meta
        if test -e {output.out_genome}; then
          echo "Flye assembly completed successfully!"
        else
          echo "Flye assembly failed, please check the logs for issues."
    	fi
        """


rule medaka:
    input:
        long_filter=lambda wildcards: "{work_path}/output/clean_long/{sample}.fastq".format(work_path=work_dir, sample=wildcards.sample),
        genome=lambda wildcards: "{work_path}/output/flye/{sample}/assembly.fasta".format(work_path=work_dir, sample=wildcards.sample)
    output:
        medaka_genome="{work_path}/output/medaka/{sample}/consensus.fasta"
    threads:8
    params: 
        out_path="{work_path}/output/medaka/{sample}"
    conda:
        "envs/medaka.yaml"
    wildcard_constraints:
        sample="|".join(samples)
    shell:
        """
        medaka_consensus -t {threads} -i {input.long_filter} \
         -d {input.genome} \
         -o {params.out_path} -m r941_min_high_g360
        """

rule polish_NGS:
    input:
        medaka_genome=lambda wildcards: "{work_path}/output/medaka/{sample}/consensus.fasta".format(work_path=work_dir,sample=wildcards.sample),
        fq1=lambda wildcards: "{work_path}/output/clean_data/{sample}/{sample}_final_1.fastq".format(work_path=work_dir,sample=wildcards.sample),
        fq2=lambda wildcards: "{work_path}/output/clean_data/{sample}/{sample}_final_2.fastq".format(work_path=work_dir,sample=wildcards.sample)
    output:
        sam1="{work_path}/output/polish/{sample}/{sample}_1.sam",
        sam2="{work_path}/output/polish/{sample}/{sample}_2.sam",
        genome="{work_path}/output/polish/{sample}/{sample}.fasta"
    threads:6
    conda:
        "envs/ont.yaml"
    shell:
        """
        bwa index {input.medaka_genome}
        bwa mem -t {threads} -a {input.medaka_genome} \
         {input.fq1} > {output.sam1}
        bwa mem -t {threads} -a {input.medaka_genome} \
         {input.fq2} > {output.sam2}
        polypolish {input.medaka_genome} {output.sam1} {output.sam2} > {output.genome}
        """

rule checkm:
    input:
        genome=lambda wildcards: "{work_path}/output/polish/{sample}/{sample}.fasta".format(sample=wildcards.sample, work_path=work_dir)
    output:
        tab="{work_path}/output/genome_control/checkm/{sample}.tab"
    threads:8
    conda:
         "envs/genome_control.yaml"
    wildcard_constraints:
        sample="|".join(samples)
    params:
        input_path="{work_path}/output/polish/{sample}",
        output_path="{work_path}/output/genome_control/checkm"
    shell:
        """
        checkm lineage_wf -t {threads} -x fasta \
        {params.input_path}  {params.output_path} \
        --tab_table -f {output.tab}
        """

#这里的long read是ONT数据，如果想实现pacbio的功能，需要将quast等的部分参数写成字典的模式来储存键值对
rule quast:
    input:
        genome=lambda wildcards: "{work_path}/output/polish/{sample}/{sample}.fasta".format(sample=wildcards.sample, work_path=work_dir),
        long_fastq=lambda wildcards: "{work_path}/output/filtlong/{sample}.fastq.gz".format(sample=wildcards.sample, work_path=work_dir)
    output:
        report="{work_path}/output/genome_control/quast/{sample}/report.txt"
    threads:4
    wildcard_constraints:
        sample="|".join(samples)
    conda:
        "envs/genome_control.yaml"
    params:
        ref=lambda wildcards: ref_fna[wildcards.sample],
        output="{work_path}/output/genome_control/quast/{sample}"
    shell:
        """
        quast {input.genome} -t {threads} -r {params.ref} -m 400 \
         --nanopore {input.long_fastq} -o {params.output}
        """

'''
rule multiqc_genome:
    input:
        quast=lambda wildcards: "{work_path}/output/genome_control/quast/{sample}/report.txt".format(sample=wildcards.sample, work_path=work_dir)
        #checkm=lambda wildcards: "{work_path}/output/genome_control/checkm/mess.tab".format(work_path=work_dir)
    output:
        html="{work_path}/output/multiqc/genome_judge/sum.html"
    threads: 4
    wildcard_constraints:
        sample="|".join(samples)
    conda:
        "envs/genome_control.yaml"
    params:
        path="{work_path}/output/genome_control",
        output_path="{work_path}/output/multiqc/genome_judge"
    shell:
        """
        cd {params.path}/
        multiqc . -o {params.output_path}/ --filename {output.html}
        """
'''

checkpoint quast_completed:
    input:
        expand("{work_path}/output/genome_control/quast/{sample}/report.txt", sample=samples, work_path=[work_dir])
    output:
        touch("{work_path}/output/quast_completed.txt")
    shell:
        """
        """

rule multiqc_genome:
    input:
        quast_completed=expand("{work_path}/output/quast_completed.txt", work_path=[work_dir])
    output:
        html="{work_path}/output/multiqc/genome_judge/sum.html"
    threads: 4
    conda:
        "envs/genome_control.yaml"
    params:
        path="{work_path}/output/genome_control",
        output_path="{work_path}/output/multiqc/genome_judge"
    shell:
        """
        cd {params.path}/
        multiqc . -o {params.output_path}/ --filename {output.html}
        """

rule parsnpANDnucmer:
    input:
        genome=lambda wildcards: "{work_path}/output/polish/{sample}/{sample}.fasta".format(work_path=work_dir,sample=wildcards.sample),
    output:
        filter="{work_path}/output/nucmer/{sample}.filter",
        snps="{work_path}/output/nucmer/{sample}.snps",
        nucmer_vcf="{work_path}/output/nucmer/{sample}_nucmer.vcf",
        par_vcf="{work_path}/output/parsnp/{sample}/parsnp.vcf"
    params:
        script_dir=os.path.dirname(BASEDIR),
        ref=lambda wildcards: ref_fna[wildcards.sample],
        parsnp_path="{work_path}/output/parsnp/{sample}",
        nucmer_path="{work_path}/output/nucmer/{sample}",
    wildcard_constraints:
        sample="|".join(samples)
    threads:6
    conda:
        "envs/parsnp.yaml"
    shell:
        """
        parsnp -c -r {params.ref} \
         -d {input.genome} \
         {input.genome} \
         -p {threads} --vcf -o {params.parsnp_path} --skip-phylogeny
        if test -e {output.par_vcf}; then
          echo "Parsnp completed successfully!"
        else
          echo "Parsnp failed, please check the logs for issues."
        fi
        nucmer {params.ref} {input.genome} -p {params.nucmer_path} 
        delta-filter -1 -q -r {params.nucmer_path}.delta \
         > {output.filter}
        show-snps -Clr -x 1 -T \
         {output.filter} \
         > {output.snps}
        python {params.script_dir}/scripts/mummer2vcf.py {output.snps} {output.nucmer_vcf}
        """

rule build_truth:
    input:
        par_vcf=lambda wildcards: "{work_path}/output/parsnp/{sample}/parsnp.vcf".format(work_path=work_dir,sample=wildcards.sample),
        nuc_vcf=lambda wildcards: "{work_path}/output/nucmer/{sample}_nucmer.vcf".format(work_path=work_dir,sample=wildcards.sample)
    output:
        truth_vcf="{work_path}/output/truth_vcf/truth_{sample}.vcf",
        chro_bed="{work_path}/output/truth_bed/{sample}/chromosome.bed",
        inco_bed="{work_path}/output/truth_bed/{sample}/inconsistent_variants.bed"
    params:
        bed_path="{work_path}/output/truth_bed/{sample}",
        ref=lambda wildcards: ref_fna[wildcards.sample]
    threads: 1
    wildcard_constraints:
        sample="|".join(samples)
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        sed -i "1s/^/##fileformat=VCFv4.1\\n/" {input.par_vcf}
        bgzip {input.par_vcf}
        bgzip {input.nuc_vcf}
        bcftools index {input.par_vcf}.gz
        bcftools index {input.nuc_vcf}.gz
        bcftools isec -n=2 -w 1 -O v -o {output.truth_vcf} \
         {input.nuc_vcf}.gz \
         {input.par_vcf}.gz
        bcftools isec -C -w 1 -O v -o {params.bed_path}/inconsistent_variants.vcf \
         {input.nuc_vcf}.gz {input.par_vcf}.gz
        awk -F '\\t' 'BEGIN {{OFS="\\t"}} /^#/ {{next}} {{print $1, $2-1, $2}}' \
         {params.bed_path}/inconsistent_variants.vcf > {params.bed_path}/inconsistent_variants.bed
        awk '{{print $1 "\\t0\\t" $2}}' {params.ref}.fai > {params.bed_path}/chromosome.bed
        if test -e {output.chro_bed} && test -e {output.inco_bed}; then
          echo "Construct bed files completed successfully!"
        else
          echo "Construct bed files failed, please check the logs for issues."
        fi
        """

rule bedANDhap:
    input:
        chromosome_bed=lambda wildcards: "{work_path}/output/truth_bed/{sample}/chromosome.bed".format(work_path=work_dir,sample=wildcards.sample),
        inconsistent_bed=lambda wildcards: "{work_path}/output/truth_bed/{sample}/inconsistent_variants.bed".format(work_path=work_dir,sample=wildcards.sample),
        truth_vcf=lambda wildcards: "{work_path}/output/truth_vcf/truth_{sample}.vcf".format(work_path=work_dir,sample=wildcards.sample)
    output:
        truth_bed="{work_path}/output/truth_bed/{sample}/truth_{sample}.bed"
    params:
        ref=lambda wildcards: ref_fna[wildcards.sample]
    threads: 1
    wildcard_constraints:
        sample="|".join(samples)
    conda:
        "envs/hap.yaml"
    shell:
        """
        bedtools subtract -a {input.chromosome_bed} -b {input.inconsistent_bed} > {output.truth_bed}
        """


rule skani:
    input:
        genome=expand("{work_path}/output/unicycler/{sample}/{sample}.fasta", sample=samples, work_path=work_dir)
    output:
        sample_ani="{work_path}/output/skani/sample_matrix.txt",
        database_ani="{work_path}/output/skani/sample_db_matrix.txt",
        full_matrix="{work_path}/output/skani/full_matrix.txt",
        ani_group="{work_path}/output/skani/ani_group.txt",
        snp_data="{work_path}/output/skani/snp_mess.txt"
    threads: 32
    conda:
        "envs/skani.yaml"
    params:
       genome_db=lambda wildcards: "{genome_db}".format(genome_db=genome_database),
       model=lambda wildcards: "{model}".format(model=snp_model),
       path=lambda wildcards: "{path}".format(path=genome_database),
       script_dir=os.path.dirname(BASEDIR)
    shell:
        """
        skani triangle -t {threads} --full-matrix --no-learned-ani \
          {input.genome} -o {output.sample_ani}
        skani triangle -t {threads} --full-matrix --no-learned-ani \
          {params.genome_db}/* -o {output.database_ani}
        python {params.script_dir}/scripts/skani_matrix.py -i {output.database_ani}.af -o {output.full_matrix}
        python {params.script_dir}/scripts/skani_group.py -i {output.full_matrix} -o {output.ani_group}
        python {params.script_dir}/scripts/snp_cal_matrix.py -i {output.ani_group} -o {output.snp_data} -path {params.path} --model {params.model}
        """


rule cut_snp:
    input:
        expand("{work_path}/output/skani/snp_mess.txt", work_path=work_dir)
    output:
        expand("{work_path}/output/skani/snp_mess_part_{i}.txt", i=range(1, snpfile_num + 1), work_path=work_dir)
    threads: 24
    conda:
        "envs/skani.yaml"
    params:
        num_files=snpfile_num,
        script_dir=os.path.dirname(BASEDIR),
        output_prefix=lambda wildcards: "{work_path}/output/skani/snp_mess_part".format(work_path=work_dir)
    shell:
        """
        python {params.script_dir}/scripts/snp_split.py -n {params.num_files} -i {input} -o {params.output_prefix}
        for i in $(seq 1 {params.num_files}); do
            if [ ! -f "{params.output_prefix}_$i.txt" ]; then
                echo "Error: Output file does not exist: {params.output_prefix}_$i.txt"
                exit 1
            fi
        done
        for i in $(seq 1 {params.num_files}); do
            touch {params.output_prefix}_$i.txt.done
        done
        """        


rule snippy:
    input:
        txt="{work_path}/output/skani/snp_mess_part_{i}.txt"
    output:
        vcf_txt="{work_path}/output/snippy/job_{i}.txt"
    threads: 8
    conda:
       "envs/snippy.yaml"
    params:
       snippy_path=lambda wildcards: "{work_path}/output/snippy".format(work_path=work_dir)
    shell:
       """
        while read -r line; do
            strain1=$(echo "$line" | cut -f1)
            strain2=$(echo "$line" | cut -f2)
            strain1_name=$(basename "$strain1" | sed 's/\\.[^.]*$//')
            strain2_name=$(basename "$strain2" | sed 's/\\.[^.]*$//')

            snippy --cpus {threads} --outdir {params.snippy_path}/${{strain1_name}}_${{strain2_name}} --ref ${{strain1}} --ctgs ${{strain2}}

        done < {input.txt}
        touch {output.vcf_txt}
       """


rule snp_group:
    input:
        txt=expand("{work_path}/output/snippy/job_{i}.txt", work_path=work_dir, i=range(1, snpfile_num + 1))
    output:
        snpmess_txt="{work_path}/output/snp.txt",
        snp_group_low="{work_path}/output/snp_group_low.txt",
        snp_group_high="{work_path}/output/snp_group_high.txt",
        data="{work_path}/output/input_genome.txt"
    threads: 12  
    params:
        script_dir=os.path.dirname(BASEDIR),
        snp_path=lambda wildcards: "{work_path}/output/skani".format(work_path=work_dir),
        snippy_path=lambda wildcards: "{work_path}/output/snippy".format(work_path=work_dir),
        uni_path=lambda wildcards: "{work_path}/output/unicycler".format(work_path=work_dir),
        model=lambda wildcards: "{model}".format(model=snp_model),
        type=lambda wildcards: "{iden_type}".format(iden_type=calcu_type),
        num1=lambda wildcards: "{num1}".format(num1=cut_num1),
        num2=lambda wildcards: "{num2}".format(num2=cut_num2)
    shell:
      """
      rm {input.txt}
      bash {params.script_dir}/scripts/gain_snippy.sh -i {params.snp_path}/snp_mess.txt -p {params.snippy_path} -o {output.snpmess_txt}
      python {params.script_dir}/scripts/snp_group.py -i {output.snpmess_txt} -m {params.model} -t {params.type} \
        -n {params.num1} -o {output.snp_group_low}
      python {params.script_dir}/scripts/snp_group.py -i {output.snpmess_txt} -m {params.model} -t {params.type} \
        -n {params.num2} -o {output.snp_group_high}
      cd {params.uni_path}
      ls | grep fasta >> {output.data}
      """


rule strain_mess:
    input:
        sample=expand("{work_path}/output/input_genome.txt", work_path=work_dir),
        ani_group=expand("{work_path}/output/skani/ani_group.txt", work_path=work_dir),
        snp_group_low=expand("{work_path}/output/snp_group_low.txt", work_path=work_dir),
        snp_group_high=expand("{work_path}/output/snp_group_high.txt", work_path=work_dir)
    output:
        strain_low="{work_path}/output/strain_mess_low.txt",
        strain_high="{work_path}/output/strain_mess_high.txt",
        database_low="{work_path}/output/database_low.txt",
        database_high="{work_path}/output/database_high.txt",
        non_redundant_low="{work_path}/output/non_redundant_low.txt",
        non_redundant_high="{work_path}/output/non_redundant_high.txt"
    threads: 12
    params:
        script_dir=os.path.dirname(BASEDIR),
        genomedb_path=lambda wildcards: "{genome_db}".format(genome_db=genome_database),
        path=lambda wildcards: "{work_path}/output".format(work_path=work_dir)
    shell:
      """
      python {params.script_dir}/scripts/strain_mess.py -i {input.sample} -a {input.ani_group} -s {input.snp_group_low} \
        --newstrain {output.strain_low} --database {output.database_low}
      python {params.script_dir}/scripts/strain_mess.py -i {input.sample} -a {input.ani_group} -s {input.snp_group_high} \
        --newstrain {output.strain_high} --database {output.database_high}
      python {params.script_dir}/scripts/database_dun.py -i {output.database_low} -o {output.non_redundant_low}
      python {params.script_dir}/scripts/database_dun.py -i {output.database_high} -o {output.non_redundant_high}
      mkdir -p {params.path}/non_redundant_low_database
      for i in $(cat {output.non_redundant_low}); do
        cd {params.genomedb_path}
        cp $i {params.path}/non_redundant_low_database
      done
      """



