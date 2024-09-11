import os
import sys
import argparse

BASEDIR = os.path.dirname(os.path.abspath(__file__))


default_configfile = config["default_configfile"]
parser = argparse.ArgumentParser()
parser.add_argument("--configfile", default=default_configfile, help="Path to the configuration file")
args, unknown = parser.parse_known_args()

configfile: args.configfile

samples = list(map(str, config["samples"]))
nums = config["nums"]
work_dir = config["work_dir"]
fastq_dir = config["fastq_dir"]
name = config["name"]
ref_fna = config["ref_fna"]
name = {str(k): str(v) for k, v in config["name"].items()}  # 将name中的键和值都转换为字符串
ref_fna = {str(k): str(v) for k, v in config["ref_fna"].items()}  # 将ref_fna中的键和值都转换为字符串
filtered_genus_name = config["filtered_genus_name"]
genome_database = config["genome_database"]
snp_model = config["snp_model"]
snpfile_num = config["snp_num"]
calcu_type = config["calcu_type"]
cut_num1 = config["cut_num1"]
cut_num2 = config["cut_num2"]


rule all:
    input:
        expand("{work_path}/output/unicycler/{sample}/{sample}.fasta", work_path=work_dir, sample=samples),
        expand("{work_path}/output/skani/sample_db_matrix.txt", work_path=work_dir),
        expand("{work_path}/output/skani/snp_mess.txt", work_path=work_dir),
        expand("{work_path}/output/skani/snp_mess_part_{i}.txt", i=range(1, snpfile_num + 1), work_path=work_dir),
        expand("{work_path}/output/snippy/job_{i}.txt", i=range(1, snpfile_num + 1), work_path=work_dir),
        expand("{work_path}/output/snp_group_low.txt", work_path=work_dir),
        expand("{work_path}/output/non_redundant_high.txt", work_path=work_dir)


rule trim:
    input:
        fastq1=lambda wildcards: "{fastq_path}/{srr}_1.fastq".format(srr=name[wildcards.sample], fastq_path=fastq_dir),
        fastq2=lambda wildcards: "{fastq_path}/{srr}_2.fastq".format(srr=name[wildcards.sample], fastq_path=fastq_dir)
    output:
        fastq1="{work_path}/output/fastp/{sample}/{sample}_trim_1.fastq",
        fastq2="{work_path}/output/fastp/{sample}/{sample}_trim_2.fastq"
    threads: 8
    conda:
        "envs/fastp.yaml"
    wildcard_constraints:
        sample="|".join(samples)
    shell:
        """
        fastp --thread {threads} -D --dup_calc_accuracy 4 --detect_adapter_for_pe \
         -i {input.fastq1} \
         -I {input.fastq2} \
         -o {output.fastq1} \
         -O {output.fastq2}
        """


rule kraken2:
    input:
        fastq1=lambda wildcards: "{work_path}/output/fastp/{sample}/{sample}_trim_1.fastq".format(sample=wildcards.sample, work_path=work_dir),
        fastq2=lambda wildcards: "{work_path}/output/fastp/{sample}/{sample}_trim_2.fastq".format(sample=wildcards.sample, work_path=work_dir)
    output:
        report="{work_path}/output/kraken2/{sample}.report",
        out="{work_path}/output/kraken2/{sample}.output",
        readlist1="{work_path}/output/kraken2/{sample}_p1.readlist",
        readlist2="{work_path}/output/kraken2/{sample}_p2.readlist"
    threads: 8
    params:
        kraken_db=lambda wildcards: "{work_path}/output/kraken2/db".format(work_path=work_dir),
        filter_name=lambda wildcards: " ".join(filtered_genus_name[wildcards.sample]),
        dbname=expand("{work_path}/output/kraken2/db", work_path=[work_dir])
    wildcard_constraints:
        sample="|".join(samples)
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        if [ ! -f {params.dbname}/taxo.k2d ]; then
          mkdir -p {params.dbname}
          cd {params.dbname}
          wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20231009.tar.gz
          tar zxvf k2_standard_16gb_20231009.tar.gz
        fi
        kraken2 --threads {threads} --db {params.kraken_db} \
         --use-names --report {output.report} --output {output.out} \
         --paired {input.fastq1} {input.fastq2}
        for filter_name in {params.filter_name}; do
          fgrep $filter_name {output.out} | cut -f2 >> {output.readlist1}
        done
        # 使用sort命令去除readlist1中的重复行
        sort -u {output.readlist1} -o {output.readlist1}.dup
        rm {output.readlist1}
        mv {output.readlist1}.dup {output.readlist1}
        cp {output.readlist1} {output.readlist2}
        #如果原始fastq有/1等信息，就激活下面两行和seqtk的注释行
        #sed -i 's/$/\/1/' {output.readlist1}
        cp {output.readlist1} {output.readlist2}
        #如果原始fastq有/1等信息，就激活下面两行和seqtk的注释行
        #sed -i 's/$/\/1/' {output.readlist1}
        #sed -e 's/1$/2/g' {output.readlist1} > {output.readlist2}
        """


rule seqtk:
    input:
        fastq1=lambda wildcards: "{work_path}/output/fastp/{sample}/{sample}_trim_1.fastq".format(sample=wildcards.sample, work_path=work_dir),
        fastq2=lambda wildcards: "{work_path}/output/fastp/{sample}/{sample}_trim_2.fastq".format(sample=wildcards.sample, work_path=work_dir),
        readlist1=lambda wildcards: "{work_path}/output/kraken2/{sample}_p1.readlist".format(sample=wildcards.sample, work_path=work_dir),
        readlist2=lambda wildcards: "{work_path}/output/kraken2/{sample}_p2.readlist".format(sample=wildcards.sample, work_path=work_dir)
    output:
        fastq1="{work_path}/output/seqtk/{sample}/{sample}_filter_1.fastq",
        fastq2="{work_path}/output/seqtk/{sample}/{sample}_filter_2.fastq"
    conda:
        "envs/ont.yaml"
    wildcard_constraints:
        sample="|".join(samples)
    shell:
        """
        seqtk subseq {input.fastq1} {input.readlist1} > {output.fastq1}
        seqtk subseq {input.fastq2} {input.readlist2} > {output.fastq2}
        #seqtk subseq {input.fastq1} {input.readlist1} | sed -e 's/\.1 / /g' > {output.fastq1}
        #seqtk subseq {input.fastq2} {input.readlist2} | sed -e 's/\.2 / /g' > {output.fastq2}
        """


rule unicycler:
    input:
        fastq1=lambda wildcards: "{work_path}/output/seqtk/{sample}/{sample}_filter_1.fastq".format(sample=wildcards.sample, work_path=work_dir),
        fastq2=lambda wildcards: "{work_path}/output/seqtk/{sample}/{sample}_filter_2.fastq".format(sample=wildcards.sample, work_path=work_dir)
    output:
        genome="{work_path}/output/unicycler/{sample}/{sample}.fasta"
    threads: 8
    params:
        outpath=lambda wildcards: "{work_path}/output/unicycler/{sample}".format(sample=wildcards.sample, work_path=work_dir),
        unipath=lambda wildcards: "{work_path}/output/unicycler".format(work_path=work_dir),
        db=lambda wildcards: "{genome_db}".format(genome_db=genome_database)
    conda:
        "envs/unicycler.yaml" 
    wildcard_constraints:
        sample="|".join(samples)
    shell:
        """
        unicycler -t {threads} --keep 1 --mode normal \
          -1 {input.fastq1} -2 {input.fastq2} -o {params.outpath}
        mv {params.outpath}/assembly.fasta {output.genome}
        cp {output.genome} {params.unipath}
        cp {output.genome} {params.db}
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



