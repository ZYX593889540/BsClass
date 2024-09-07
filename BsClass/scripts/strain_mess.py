import pandas as pd
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description="Process input genome data, ANI group data and SNP group data.")
parser.add_argument('-i', '--input_genome', required=True, help='Input file path for genome names.')
parser.add_argument('-a', '--ani_group', required=True, help='Input file path for ANI group results.')
parser.add_argument('-s', '--snp_group', required=True, help='Input file path for SNP group results.')
parser.add_argument('-o1', '--newstrain', required=True, help='Output file path for new strain result.')
parser.add_argument('-o2', '--database', required=True, help='Output file path for database result.')
args = parser.parse_args()

# 读取输入文件
input_genomes = pd.read_csv(args.input_genome, header=None, names=['strain'])
ani_group = pd.read_csv(args.ani_group, sep='\t')
snp_group = pd.read_csv(args.snp_group, sep='\t', index_col=0)

# 确保所有基因组名称统一格式
input_genomes['strain'] = input_genomes['strain'].apply(lambda x: x.strip())
ani_group['Strain'] = ani_group['Strain'].apply(lambda x: x.strip())
snp_group.index = snp_group.index.map(lambda x: x.strip() + '.fasta' if not x.endswith('.fasta') else x.strip())

# 结果1处理：输入样本的新菌株判断
result1 = []
for strain in input_genomes['strain']:
    ani_group_strains = ani_group[ani_group['Strain'] == strain]
    if ani_group_strains.empty:
        continue
    ani_group_id = ani_group_strains['ANI_group'].values[0]
    same_ani_group = ani_group[ani_group['ANI_group'] == ani_group_id]

    # 判断是否为新菌株
    if len(same_ani_group) == 1:
        new_strain = "yes"
        notes = "None"
    else:
        if strain in snp_group.index:
            snp_group_id = snp_group.loc[strain].values[0]
            same_snp_group = snp_group[snp_group.iloc[:, 0] == snp_group_id]
            same_snp_group_strains = same_snp_group.index[same_snp_group.index != strain].tolist()
        else:
            same_snp_group_strains = []

        if len(same_snp_group_strains) == 0:
            new_strain = "yes"
            notes = "None"
        else:
            new_strain = "no"
            notes = ' '.join(same_snp_group_strains)

    result1.append([strain, new_strain, notes])

# 结果2处理：数据库样本的重复菌株判断
input_genome_set = set(input_genomes['strain'])
database_genomes = ani_group[~ani_group['Strain'].isin(input_genome_set)]
result2 = []
for strain in database_genomes['Strain']:
    ani_group_strains = ani_group[ani_group['Strain'] == strain]
    ani_group_id = ani_group_strains['ANI_group'].values[0]
    same_ani_group = ani_group[ani_group['ANI_group'] == ani_group_id]

    # 排除输入样本后的同组菌株
    same_ani_group = same_ani_group[~same_ani_group['Strain'].isin(input_genome_set)]
    
    # 如果只有一个数据库样本，不考虑新菌株
    if len(same_ani_group) == 1:
        duplicate_strain = "no"
        notes = "None"
    else:
        if strain in snp_group.index:
            snp_group_id = snp_group.loc[strain].values[0]
            same_snp_group = snp_group[snp_group.iloc[:, 0] == snp_group_id]
            same_snp_group_strains = same_snp_group.index[same_snp_group.index != strain].tolist()
            same_snp_group_strains = [s for s in same_snp_group_strains if s not in input_genome_set]
        else:
            same_snp_group_strains = []

        if len(same_snp_group_strains) == 0:
            duplicate_strain = "no"
            notes = "None"
        else:
            duplicate_strain = "yes"
            notes = ' '.join(same_snp_group_strains)

    result2.append([strain, duplicate_strain, notes])

# 保存结果文件
result1_df = pd.DataFrame(result1, columns=['Strain', 'New_Strain', 'Notes'])
result1_df.to_csv(args.newstrain, sep='\t', index=False)

result2_df = pd.DataFrame(result2, columns=['Strain', 'Duplicate_Strain', 'Notes'])
result2_df.to_csv(args.database, sep='\t', index=False)

print(f"Results are saved to {args.newstrain} and {args.database}")
