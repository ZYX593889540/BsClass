import pandas as pd
import argparse

# 更新群组内所有样本对的ANI平均值的函数
def update_group_ani(group, df):
    """更新群组内所有样本对的ANI平均值"""
    if len(group) <= 1:
        return 100.0  # 单个样本与自己的ANI视为100%
    sum_ani = 0
    count = 0
    # 遍历组内的每对样本
    for i in range(len(group)):
        for j in range(i + 1, len(group)):
            if group[i] != group[j]:  # 排除自比较
                sum_ani += df.loc[group[i], group[j]]
                count += 1
    return sum_ani / count if count > 0 else 0  # 防止除零错误

# 找到最佳群组以添加新样本的函数
def find_best_group_for_strain(strain, groups, df):
    """找到最佳群组以添加新样本，返回群组索引和更新后的平均ANI"""
    for idx, group in enumerate(groups):
        test_group = group + [strain]
        avg_ani = update_group_ani(test_group, df)
        # 新样本至少和该group内某个非自身样本大于等于99
        if avg_ani >= 99 and any(df.loc[strain, member] >= 99 for member in group if member != strain):
            return idx, avg_ani
    return None, None

# 解析命令行参数
parser = argparse.ArgumentParser(description="Process ANI matrix and group strains based on ANI values.")
parser.add_argument('-i', '--input', required=True, help='Input file path for ANI matrix.')
parser.add_argument('-o', '--output', required=True, help='Output file path for grouping results.')
args = parser.parse_args()

# 加载处理后的ANI矩阵
df = pd.read_csv(args.input, sep='\t', index_col=0)

# 初始化群组列表
groups = []

# 为每个菌株找到最佳群组
for strain in df.index:
    best_group_idx, _ = find_best_group_for_strain(strain, groups, df)
    if best_group_idx is not None:
        groups[best_group_idx].append(strain)
    else:
        groups.append([strain])

# 将群组信息写入文件
output_file = args.output
with open(output_file, 'w') as f:
    for idx, group in enumerate(groups, start=1):
        f.write(f"Group {idx}: {', '.join(group)}\n")

# 为每个基因组样本列出其对应的ANI群组编号，并保存到文件
group_list_file = args.output
with open(group_list_file, 'w') as f:
    f.write("Strain\tANI_group\n")
    for idx, group in enumerate(groups, start=1):
        for strain in group:
            f.write(f"{strain}\t{idx}\n")

print(f"Grouping results are saved to {args.output}")
