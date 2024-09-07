import argparse

def split_file(input_file, output_files_prefix, num_files, lines_per_file):
    with open(input_file, 'r') as f:
        _ = f.readline()
        lines = f.readlines()
    
    # 计算每个文件应有的行数，最后一个文件可能包含多出的行
    total_lines = len(lines)
    lines_per_file = [total_lines // num_files + (1 if i < total_lines % num_files else 0) for i in range(num_files)]

    # 分割文件
    start_index = 0
    for i in range(num_files):
        end_index = start_index + lines_per_file[i]
        with open(f"{output_files_prefix}_{i+1}.txt", 'w') as f:
            f.writelines(lines[start_index:end_index])
        start_index = end_index

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a file into multiple files.")
    parser.add_argument('-i', '--input', required=True, help='Input file to split.')
    parser.add_argument('-n', '--num-files', type=int, required=True, help='Number of output files.')
    parser.add_argument('-o', '--output-files-prefix', required=True, help='Prefix for output files.')
    args = parser.parse_args()

    split_file(args.input, args.output_files_prefix, args.num_files, None)
