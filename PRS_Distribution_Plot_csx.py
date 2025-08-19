# =======================================
# Describe: PRS plot in category trait
# Author: Weber 
# Build_Date: 2025.08.15
# Last_update:2025.08.15
# Update:

# =======================================
#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re
import math
from argparse import ArgumentParser
from statannotations.Annotator import Annotator
from sklearn import preprocessing

# set seaborn background style
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

def argparser():
    parser = ArgumentParser()
    parser.add_argument('--prs', '-r', type=str, required=True, help='PRS score file')
    parser.add_argument('--pheno', '-p', type=str, required=True, help='Phenotype file')
    parser.add_argument('--phename', '-m', type=str, required=True, help='Phenotype column name in pheno file')
    parser.add_argument('--normalize', '-a', type=str, default=False, help='Normalization method [z_std|min_max|arctan]')
    parser.add_argument('--out', '-o', type=str, required=True, help='Output prefix')
    return parser.parse_args()

if __name__ == '__main__':
    args = argparser()

    # 1. 讀 PRS file (FID, IID, PRS)
    prs_df = pd.read_csv(args.prs, sep=r'\s+', header=None)
    
    prs_df = prs_df.iloc[:, [1, -1]]  # 取 IID, PRS
    prs_df.columns = ['IID', 'PRS']

    # 2. 讀 phenotype file
    pheno_df = pd.read_csv(args.pheno, sep="\t", usecols=[0, 1], names=["IID", args.phename])

    # 3. 合併
    merged_df = pd.merge(prs_df, pheno_df, left_on='IID', right_on='IID', how='inner')

    # 4. phenotype 轉成字串
    merged_df[args.phename] = merged_df[args.phename].astype(int).astype(str)

    # 5. PRS 正規化
    if args.normalize == 'z_std':
        merged_df['PRS'] = (merged_df.PRS - merged_df.PRS.mean()) / merged_df.PRS.std()
    elif args.normalize == 'min_max':
        scaler = preprocessing.MinMaxScaler(feature_range=(-1, 1))
        merged_df['PRS'] = scaler.fit_transform(merged_df[['PRS']])
    elif args.normalize == 'arctan':
        merged_df['PRS'] = merged_df.PRS.apply(lambda x: math.atan(x) * (2 / math.pi))

    # 6. 畫圖
    order = ['1', '2']  # 1=Control, 2=Case
    fig, ax = plt.subplots(1, 2, figsize=(10, 6))

    kde1 = sns.kdeplot(x=merged_df[merged_df[args.phename] == '1'].PRS,
						ax=ax[0],color="blue",cmap='Blues',fill=True)
    kde2 = sns.kdeplot(x=merged_df[merged_df[args.phename] == '2'].PRS,
						ax=ax[0],color="red",cmap='Reds',fill=True)
    handles = [mpatches.Patch(facecolor=plt.cm.Blues(100), label="Control"),
				mpatches.Patch(facecolor=plt.cm.Reds(100), label="Case")]
    ax[0].legend(title=args.phename,handles=handles,loc='upper right',shadow=True,edgecolor='black')

    prs_plot = sns.boxplot(data=merged_df, y='PRS', x=args.phename, order=order, ax=ax[1], hue=args.phename)
    prs_plot.set_xticklabels(['Control', 'Case'])

    annotator = Annotator(ax=prs_plot, pairs=[('1', '2')], order=order,
                          data=merged_df, x=args.phename, y='PRS')
    annotator.configure(test='Mann-Whitney', text_format='star',
                         loc='inside', comparisons_correction='fdr_bh')
    annotator.apply_and_annotate()

    fig.savefig(f'{args.out}.prs_dist.tiff', dpi=600)
