import pandas as pd
import numpy as np
import argparse

# 定义一个函数用于处理菌株信息，提取'/'后的部分
def process_strain_name(name):
    # 去除"//"中的内容，只保留文件名
    return name.split('/')[-1]

# 设置命令行参数解析
parser = argparse.ArgumentParser(description='Process a skani ANI matrix file.')
parser.add_argument('--input', '-i', required=True, help='Path to the input matrix file')
parser.add_argument('--output', '-o', required=True, help='Path to the output matrix file')
args = parser.parse_args()

# 读取数据
with open(args.input, 'r') as f:
    lines = f.readlines()

# 假设第一行是表头，我们跳过它
lines = lines[1:]

# 提取并处理菌株信息
names = [process_strain_name(line.strip().split()[0]) for line in lines]

# 检查 names 列表的长度是否与剩余行数一致
if len(names) != len(lines):
    raise ValueError("The number of names does not match the number of rows in the data.")

# 提取矩阵数据
data = []
for line in lines:
    row = line.strip().split()[1:]  # 跳过文件名，只获取相似度数据
    data.append([float(value) for value in row])

# 检查 data 的列数是否一致
if not all(len(row) == len(data[0]) for row in data):
    raise ValueError("The number of columns in the data is not consistent.")

# 转换为DataFrame
df = pd.DataFrame(data, index=names, columns=names)

# 将对角线上的值设置为100
np.fill_diagonal(df.values, 100)

# 处理缺失值，这里假设缺失值用np.nan表示，将它们替换为70.0
df.fillna(70.0, inplace=True)

# 保存到新的文件
df.to_csv(args.output, sep='\t')

print(f"The symmetric matrix is generated and saved to {args.output}")
