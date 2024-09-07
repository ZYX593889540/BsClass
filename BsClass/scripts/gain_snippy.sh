#!/bin/bash

usage() {
    echo "Usage: $0 -i <input_file> -p <snippy_working_directory> -o <output_file>"
    echo ""
    echo "This script extracts Variant-SNP and VariantTotal values from snippy results."
    echo ""
    echo "Options:"
    echo "  -i <input_file>              Path to the snp_mess.txt file."
    echo "  -p <snippy_working_directory> Path to the snippy working directory."
    echo "  -o <output_file>             Path to the output file."
    echo "  -h                           Display this help message."
    exit 1
}

# 解析参数
while getopts ":i:p:o:h" opt; do
    case ${opt} in
        i )
            input_file=$OPTARG
            ;;
        p )
            snippy_working_directory=$OPTARG
            ;;
        o )
            output_file=$OPTARG
            ;;
        h )
            usage
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# 检查是否所有必需的参数都已提供
if [ -z "${input_file}" ] || [ -z "${snippy_working_directory}" ] || [ -z "${output_file}" ]; then
    echo "Missing required arguments" 1>&2
    usage
fi

# 初始化输出文件并写入表头
echo -e "strain1\tstrain2\tvariant_snp\tvariant_total" > "$output_file"

# 读取输入文件的每一行
while IFS=$'\t' read -r strain1 strain2; do
    strain1_name=$(basename "$strain1" | sed 's/\.[^.]*$//')
    strain2_name=$(basename "$strain2" | sed 's/\.[^.]*$//')
    directory="${snippy_working_directory}/${strain1_name}_${strain2_name}"
    snps_file="${directory}/snps.txt"

    # 检查snps.txt文件是否存在
    if [[ -f "$snps_file" ]]; then
        # 提取Variant-SNP和VariantTotal的值
        variant_snp=$(grep 'Variant-SNP' "$snps_file" | cut -f 2)
        variant_total=$(grep 'VariantTotal' "$snps_file" | cut -f 2)
        if [[ -z "$variant_snp" ]]; then
            variant_snp=0
        fi

        # 写入结果到输出文件
        echo -e "${strain1}\t${strain2}\t${variant_snp}\t${variant_total}" >> "$output_file"
    else
        echo "Warning: $snps_file not found, skipping ${strain1_name}_${strain2_name}"
    fi
done < <(tail -n +2 "$input_file") # 忽略输入文件的第一行（列名）

echo "Extraction completed. Results saved to $output_file"
