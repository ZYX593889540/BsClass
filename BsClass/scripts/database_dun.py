import pandas as pd
import argparse

def extract_non_redundant(database_file, output_file):
    df = pd.read_csv(database_file, sep='\t')
    non_redundant_strains = []

    # 创建一个集合用于存储所有已经处理的样本
    processed_strains = set()

    for index, row in df.iterrows():
        strain = row['Strain']
        duplicate_strain = row['Duplicate_Strain']
        notes = row['Notes']

        if duplicate_strain == 'no':
            non_redundant_strains.append(strain)
        else:
            if strain not in processed_strains:
                non_redundant_strains.append(strain)
                if notes != 'None':
                    notes_strains = notes.split()
                    processed_strains.update(notes_strains)

    # 将结果写入输出文件
    with open(output_file, 'w') as f:
        for strain in non_redundant_strains:
            f.write(f"{strain}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract non-redundant strains from database results.")
    parser.add_argument('-i', '--input', required=True, help='Input database result file path.')
    parser.add_argument('-o', '--output', required=True, help='Output file path for non-redundant strains.')
    args = parser.parse_args()
    
    extract_non_redundant(args.input, args.output)
