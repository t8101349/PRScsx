import pandas as pd
import sys
import numpy as np

if len(sys.argv) != 4:
    print("Usage: python script.py <input_file> <format_type> <output_prefix>")
    sys.exit(1)

input_file = sys.argv[1]
format_type = sys.argv[2]
output_prefix = sys.argv[3]

print(f"Reading file: {input_file}")
df = pd.read_csv(input_file, sep='\t', low_memory=False)

if 'OR' in df.columns:
    effective = 'OR'
elif 'BETA' in df.columns:
    effective = 'BETA'
else:
    print("Error: Missing effect size column (OR or BETA)")
    sys.exit(1)

if format_type == 'plink_v1':
    df.rename({'#CHROM':'CHR'}, axis=1, inplace=True)
    df2 = df[['CHR','SNP', 'A1', 'A2', effective, 'P']].copy()
elif format_type == 'plink_v2':
    df.rename({'#CHROM':'CHR','ID': 'SNP'}, axis=1, inplace=True)
    df['A2'] = np.where(df.A1 == df.ALT, df.REF, df.ALT)
    df2 = df[['CHR','SNP', 'A1', 'A2', effective, 'P']].copy()
elif format_type == 'SAIGE':
    df.rename({'MarkerID': 'SNP', 'p.value': 'P', 'Allele1': 'A2', 'Allele2': 'A1'}, axis=1, inplace=True)
    df[['CHR', 'SNP']] = df['SNP'].str.extract(r'^(.*?):(.*?)_.*$')
    df2 = df[['CHR', 'SNP', 'A1', 'A2', effective, 'P']].copy()

else:
    print('Aberrant format...')
    sys.exit(1)

output_file = f"{output_prefix}.sumstats.prscs.txt"
df2.to_csv(output_file, sep='\t', index=False)
print(f"Saved to: {output_file}")
