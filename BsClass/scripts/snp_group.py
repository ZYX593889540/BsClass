import argparse
import pandas as pd
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description="Cluster SNP variation data.")
    parser.add_argument('-i', '--input_file', required=True, help='Input file containing SNP variation data.')
    parser.add_argument('-m', '--model', required=True, choices=['full', 'half'], help='Clustering model to use.')
    parser.add_argument('-t', '--type', required=True, choices=['snp', 'total'], help='Type of variation to cluster on.')
    parser.add_argument('-n', '--num', type=float, required=True, help='Threshold value for clustering.')
    parser.add_argument('-o', '--output_file', help='Output file to save the clustering results.')
    return parser.parse_args()

def preprocess_data(df, model, variation_type):
    if variation_type == 'snp':
        value_column = 'variant_snp'
    else:
        value_column = 'variant_total'

    df['strain1'] = df['strain1'].apply(lambda x: x.split('/')[-1].split('.')[0])
    df['strain2'] = df['strain2'].apply(lambda x: x.split('/')[-1].split('.')[0])

    if model == 'full':
        df['pair'] = df.apply(lambda x: tuple(sorted([x['strain1'], x['strain2']])), axis=1)
        df = df.groupby('pair')[value_column].mean().reset_index()
        df[['strain1', 'strain2']] = pd.DataFrame(df['pair'].tolist(), index=df.index)
        df.drop(columns=['pair'], inplace=True)
    
    return df, value_column

def cluster_data(df, value_column, threshold):
    strains = list(set(df['strain1']).union(set(df['strain2'])))
    matrix = pd.DataFrame(np.inf, index=strains, columns=strains)

    for _, row in df.iterrows():
        matrix.loc[row['strain1'], row['strain2']] = row[value_column]
        matrix.loc[row['strain2'], row['strain1']] = row[value_column]

    clusters = []
    for strain in strains:
        added = False
        for cluster in clusters:
            avg_dist = np.mean([matrix.loc[strain, s] for s in cluster])
            if avg_dist <= threshold:
                cluster.append(strain)
                added = True
                break
        if not added:
            clusters.append([strain])

    return clusters

def main():
    args = parse_arguments()
    
    df = pd.read_csv(args.input_file, sep='\t')
    df, value_column = preprocess_data(df, args.model, args.type)
    
    clusters = cluster_data(df, value_column, args.num)

    output_df = pd.DataFrame([
        {'strain': strain, f'group_{args.model}_{args.type}': i+1}
        for i, cluster in enumerate(clusters)
        for strain in cluster
    ])

    if args.output_file:
        output_file = args.output_file
    else:
        output_file = args.input_file.replace('.txt', f'_clustered_{args.model}_{args.type}.txt')

    output_df.to_csv(output_file, sep='\t', index=False)
    print(f'Clustering completed. Results saved to {output_file}')

if __name__ == '__main__':
    main()
