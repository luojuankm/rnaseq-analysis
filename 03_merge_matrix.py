import pandas as pd
import numpy as np
import os, re, argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-in_dir', required=True)
    parser.add_argument('-out_dir', required=True)
    parser.add_argument('-prefix', default='rsem_expression_count_')
    parser.add_argument('-suffix', default='.csv')
    return parser.parse_args()

def merge_files(in_dir, prefix, suffix):
    data = pd.DataFrame()
    for root, _, files in os.walk(in_dir):
        for f in files:
            if re.match(f'{prefix}.*?{suffix}', f):
                sample = f.replace(prefix, '').replace(suffix, '')
                df = pd.read_csv(os.path.join(root, f), names=['id', sample], header=0)
                df[sample] = np.floor(pd.to_numeric(df[sample], errors='coerce')).astype('Int64')
                data = df if data.empty else data.merge(df, on='id', how='outer')
    return data

def remove_duplicates(data):
    dup = data[data.duplicated('id', keep=False)]
    if dup.empty:
        return data
    cols = [c for c in data.columns if c != 'id']
    dup['sum'] = dup[cols].sum(axis=1)
    keep = dup.loc[dup.groupby('id')['sum'].idxmax()].drop('sum', axis=1)
    return pd.concat([data[~data['id'].isin(keep['id'])], keep])

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    data = merge_files(args.in_dir, args.prefix, args.suffix)
    data = remove_duplicates(data)
    data.to_csv(os.path.join(args.out_dir, 'expression_matrix.csv'), index=False)

if __name__ == '__main__':
    main()
