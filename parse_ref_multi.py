import pandas as pd

def parse_ref_multi(ref_files: dict, outfile):
    """
    合併多個族群的 reference MAF 檔，根據 SNP, BP, A1, A2 對齊。
    
    Parameters:
    - ref_files: dict，例如 {'EAS': 'eas_1kg.hm3', 'EUR': 'eur_1kg.hm3', ...}
    - outfile: str，輸出檔名
    
    Returns:
    - ref_dict: dict，包含 CHR, SNP, BP, A1, A2, FRQ_EAS, FRQ_EUR, ...
    """

    print(f"... parse reference files ...")

    pops = list(ref_files.keys())
    base_pop = pops[0]  # 第一個族群作為 reference 對齊用
    ref_dfs = {}

    for pop, file in ref_files.items():
        df = pd.read_csv(file, sep='\t')
        df = df[['CHR', 'SNP', 'BP', 'A1', 'A2', 'MAF']].copy()
        df = df.rename(columns={'MAF': f'FRQ_{pop}'})
        ref_dfs[pop] = df

    # 對所有族群 inner join：用 CHR, SNP, BP 做主鍵
    merged = ref_dfs[base_pop]
    for pop in pops[1:]:
        merged = pd.merge(merged, ref_dfs[pop], on=['CHR', 'SNP', 'BP'], suffixes=('', f'_{pop}'))

    # 計算是否需要 flip（A1/A2 不同）
    for pop in pops[1:]:
        flp_col = f'FLP_{pop}'
        merged[flp_col] = (
            (merged['A1'] != merged[f'A1_{pop}']) |
            (merged['A2'] != merged[f'A2_{pop}'])
        ).astype(int)

        # 丟掉不同的 A1/A2，統一用 base_pop 的 alleles
        merged.drop(columns=[f'A1_{pop}', f'A2_{pop}'], inplace=True)

    print(f"... {len(merged)} SNPs found in all reference files ...")

    ref_dict = {col: merged[col].tolist() for col in merged.columns}

    """
    將 ref_dict 輸出為 .hm3 格式
    
    - ref_dict: dict，從 parse_ref_multi 回傳的
    - 輸出欄位: 'CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_', 'FLP_'
    """

    df = pd.DataFrame(ref_dict)
    cols = ['CHR', 'SNP', 'BP', 'A1', 'A2']
    for pop in pops:
        freq_col = f'FRQ_{pop}'
        flp_col = f'FLP_{pop}' if pop != base_pop else None

        if freq_col in df.columns:
            cols.append(freq_col)
        if flp_col and flp_col in df.columns:
            cols.append(flp_col)

    df[cols].to_csv(outfile, sep='\t', index=False)
    print(f"輸出完成：{outfile}")
