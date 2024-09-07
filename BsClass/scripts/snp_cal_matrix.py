import pandas as pd
import argparse
from itertools import combinations, product

# 解析命令行参数
parser = argparse.ArgumentParser(description="Generate SNP comparison pairs based on ANI group.")
parser.add_argument('-i', '--input', required=True, help='Input file with group content.')
parser.add_argument('-path', '--path', required=True, help='Path to the genome data sites.')
parser.add_argument('-o', '--output', required=True, help='Output file for the SNP comparison list.')
parser.add_argument('-model', '--model', choices=['full', 'half'], required=True, help='Model for generating comparison pairs.')
args = parser.parse_args()

# 读取输入文件
df = pd.read_csv(args.input, sep='\t')

# 根据 ANI_group 分组
groups = df.groupby('ANI_group')

# 准备输出文件
output_pairs = []

# 对每个 ANI_group 进行处理
for group_id, strains in groups:
    # 获取该组的所有 strain 名称
    strains_list = strains['Strain'].tolist()
    if len(strains_list) < 2:
        continue  # 如果组的大小小于2，则无法生成比较对，跳过
    
    # 根据 -model 参数选择生成比较对的方式
    if args.model == 'full':
        # full 模型：生成所有可能的两两比较（不包括自身比较）
        for strain1, strain2 in product(strains_list, strains_list):
            if strain1 != strain2:  # 剔除自身比较
                output_pairs.append(f"{args.path}/{strain1}\t{args.path}/{strain2}")
    elif args.model == 'half':
        # half 模型：只保留左下三角部分的比较对
        for strain1, strain2 in combinations(strains_list, 2):
            output_pairs.append(f"{args.path}/{strain1}\t{args.path}/{strain2}")

# 去重并写入输出文件
unique_pairs = set(output_pairs)  # 将列表转换为集合以去除重复项
with open(args.output, 'w') as f:
    f.write("strain1\tstrain2\n")  # 写入表头
    for pair in sorted(unique_pairs):  # 排序后写入文件
        f.write(pair + '\n')

print(f"SNP comparison pairs are saved to {args.output}")
