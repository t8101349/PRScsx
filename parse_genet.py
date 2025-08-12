#!/usr/bin/env python

"""
Parse the reference panel, summary statistics, and validation set.

"""


import os
import numpy as np
from scipy.stats import norm
from scipy import linalg
import h5py


def parse_ref(ref_file, chrom, ref):
    print('... parse reference file: %s ...' % ref_file)

    # 適用簡化格式: CHR SNP BP A1 A2 MAF
    ref_dict = {'CHR': [], 'SNP': [], 'BP': [], 'A1': [], 'A2': [],
                'FRQ_EAS': [], 'FLP_EAS': []}  #使用前需調整,依所使用的population決定欄位, ref file 要先經parse_ref_multi.py處理(則可直接讀取使用)
    
    with open(ref_file) as ff:
        header = next(ff).strip().split()
        col_idx = {col: idx for idx, col in enumerate(header)}

        for line in ff:
            ll = line.strip().split()
            if int(ll[col_idx['CHR']]) == chrom:
                ref_dict['CHR'].append(ll[col_idx['CHR']])
                ref_dict['SNP'].append(ll[col_idx['SNP']])
                ref_dict['BP'].append(ll[col_idx['BP']])
                ref_dict['A1'].append(ll[col_idx['A1']])
                ref_dict['A2'].append(ll[col_idx['A2']])
                
                # 用 MAF 當作頻率，如果 MAF > 0.5，則轉成次要等位基因頻率
                maf = float(ll[col_idx['MAF']])
                frq = maf if maf <= 0.5 else 1.0 - maf
                ref_dict['FRQ_EAS'].append(frq)
                
                # flip 定義為 1，如果等位基因順序需要調整就設 -1（簡化）
                ref_dict['FLP_EAS'].append(1)

    print('... %d SNPs on chromosome %d read from %s ...' % (len(ref_dict['SNP']), chrom, ref_file))
    return ref_dict


def parse_bim(bim_file, chrom):
    print('... parse bim file: %s ...' % (bim_file + '.bim'))

    vld_dict = {'SNP':[], 'A1':[], 'A2':[]}
    with open(bim_file + '.bim') as ff:
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == chrom:
                vld_dict['SNP'].append(ll[1])
                vld_dict['A1'].append(ll[4])
                vld_dict['A2'].append(ll[5])

    print('... %d SNPs on chromosome %d read from %s ...' % (len(vld_dict['SNP']), chrom, bim_file + '.bim'))
    return vld_dict


def parse_sumstats(ref_dict, vld_dict, sst_file, pop, n_subj):
    print('... parse ' + pop.upper() + ' sumstats file: %s ...' % sst_file)

    ATGC = ['A', 'T', 'G', 'C']
    sst_dict = {'SNP':[], 'A1':[], 'A2':[]}
    with open(sst_file) as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            if ll[2] in ATGC and ll[3] in ATGC:
                sst_dict['SNP'].append(ll[1])
                sst_dict['A1'].append(ll[2])
                sst_dict['A2'].append(ll[3])

    print('... %d SNPs read from %s ...' % (len(sst_dict['SNP']), sst_file))


    idx = [ii for (ii,frq) in enumerate(ref_dict['FRQ_'+pop.upper()]) if frq>0]
    snp_ref = [ref_dict['SNP'][ii] for ii in idx]
    a1_ref = [ref_dict['A1'][ii] for ii in idx]
    a2_ref = [ref_dict['A2'][ii] for ii in idx]


    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    vld_snp = set(zip(vld_dict['SNP'], vld_dict['A1'], vld_dict['A2']))

    ref_snp = set(zip(snp_ref, a1_ref, a2_ref)) | set(zip(snp_ref, a2_ref, a1_ref)) | \
              set(zip(snp_ref, [mapping[aa] for aa in a1_ref], [mapping[aa] for aa in a2_ref])) | \
              set(zip(snp_ref, [mapping[aa] for aa in a2_ref], [mapping[aa] for aa in a1_ref]))

    sst_snp = set(zip(sst_dict['SNP'], sst_dict['A1'], sst_dict['A2'])) | set(zip(sst_dict['SNP'], sst_dict['A2'], sst_dict['A1'])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A1']], [mapping[aa] for aa in sst_dict['A2']])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A2']], [mapping[aa] for aa in sst_dict['A1']]))

    comm_snp = vld_snp & ref_snp & sst_snp

    print('... %d common SNPs in the %s reference, %s sumstats, and validation set ...' % (len(comm_snp), pop.upper(), pop.upper()))


    n_sqrt = np.sqrt(n_subj)
    sst_eff = {}
    with open(sst_file) as ff:
        header = next(ff).strip().split()
        header = [col.upper() for col in header]
        idx_snp = header.index('SNP')
        idx_a1 = header.index('A1')
        idx_a2 = header.index('A2')
        idx_beta = header.index('BETA') if 'BETA' in header else header.index('OR')
        idx_se = header.index('SE') if 'SE' in header else None
        idx_p = header.index('P') if 'P' in header else None

        for line in ff:
            ll = line.strip().split()
            snp = ll[idx_snp]
            a1 = ll[idx_a1].upper()
            a2 = ll[idx_a2].upper()
            if a1 not in ATGC or a2 not in ATGC:
                continue

            matched = False
            flip = False
            if (snp, a1, a2) in comm_snp or (snp, mapping[a1], mapping[a2]) in comm_snp:
                flip = False
                matched = True
            elif (snp, a2, a1) in comm_snp or (snp, mapping[a2], mapping[a1]) in comm_snp:
                flip = True
                matched = True

            if not matched:
                continue

            if 'BETA' in header:
                beta = float(ll[idx_beta])
            elif 'OR' in header:
                beta = np.log(float(ll[idx_beta]))

            if idx_se is not None:
                se = float(ll[idx_se])
                beta_std = beta / se / n_sqrt
            elif idx_p is not None:
                p = max(float(ll[idx_p]), 1e-323)
                beta_std = np.sign(beta) * abs(norm.ppf(p / 2.0)) / n_sqrt

            if flip:
                beta_std *= -1

            sst_eff[snp] = beta_std



    sst_dict = {'SNP': [], 'FRQ': [], 'BETA': [], 'FLP': []}
    for ii, snp in enumerate(ref_dict['SNP']):
        if snp in sst_eff:
            a1 = ref_dict['A1'][ii]
            a2 = ref_dict['A2'][ii]

            # 用 comm_snp 確認是否 allele 順序需要 flip
            if (snp, a1, a2) in comm_snp or (snp, mapping[a1], mapping[a2]) in comm_snp:
                frq = ref_dict['FRQ_' + pop.upper()][ii]
                flp = ref_dict['FLP_' + pop.upper()][ii]
            elif (snp, a2, a1) in comm_snp or (snp, mapping[a2], mapping[a1]) in comm_snp:
                frq = 1 - ref_dict['FRQ_' + pop.upper()][ii]
                flp = -1 * ref_dict['FLP_' + pop.upper()][ii]
            else:
                continue  

            sst_dict['SNP'].append(snp)
            sst_dict['BETA'].append(sst_eff[snp])
            sst_dict['FRQ'].append(frq)
            sst_dict['FLP'].append(flp)

    return sst_dict


def parse_ldblk(ldblk_dir, sst_dict, pop, chrom, ref):
    print('... parse %s reference LD on chromosome %d ...' % (pop.upper(), chrom))

    # 注意路徑組合是否需調整
    if ref == '1kg':
        chr_name = ldblk_dir + '/ldblk_1kg_chr' + str(chrom) + '.hdf5'
    elif ref == 'ukbb':
        chr_name = ldblk_dir + '/ldblk_ukbb_chr' + str(chrom) + '.hdf5'

    hdf_chr = h5py.File(chr_name, 'r')
    n_blk = len(hdf_chr)
    ld_blk = [np.array(hdf_chr['blk_'+str(blk)]['ldblk']) for blk in range(1,n_blk+1)]

    snp_blk = []
    for blk in range(1,n_blk+1):
         snp_blk.append([bb.decode("UTF-8") for bb in list(hdf_chr['blk_'+str(blk)]['snplist'])])

    blk_size = []
    mm = 0
    for blk in range(n_blk):
        idx = [ii for (ii,snp) in enumerate(snp_blk[blk]) if snp in sst_dict['SNP']]
        blk_size.append(len(idx))
        if idx != []:
            idx_blk = range(mm,mm+len(idx))
            flip = [sst_dict['FLP'][jj] for jj in idx_blk]
            ld_blk[blk] = ld_blk[blk][np.ix_(idx,idx)]*np.outer(flip,flip)

            _, s, v = linalg.svd(ld_blk[blk])
            h = np.dot(v.T, np.dot(np.diag(s), v))
            ld_blk[blk] = (ld_blk[blk]+h)/2

            mm += len(idx)
        else:
            ld_blk[blk] = np.array([])

    return ld_blk, blk_size


def align_ldblk(ref_dict, vld_dict, sst_dict, n_pop, chrom):
    print('... align reference LD on chromosome %d across populations ...' % chrom)

    # 建 index lookup 加快比對 & 安全處理
    vld_index = {snp: idx for idx, snp in enumerate(vld_dict['SNP'])}
    snp_dict = {'CHR': [], 'SNP': [], 'BP': [], 'A1': [], 'A2': []}

    for ii, snp in enumerate(ref_dict['SNP']):
        for pp in range(n_pop):
            if snp in sst_dict[pp]['SNP']:
                if snp in vld_index:
                    snp_dict['SNP'].append(snp)
                    snp_dict['CHR'].append(ref_dict['CHR'][ii])
                    snp_dict['BP'].append(ref_dict.get('BP', ['0'] * len(ref_dict['SNP']))[ii])  # 若無 BP 欄位則補 0
                    idx = vld_index[snp]
                    snp_dict['A1'].append(vld_dict['A1'][idx])
                    snp_dict['A2'].append(vld_dict['A2'][idx])
                break  # 該 SNP 在某一族群有就記一次即可

    n_snp = len(snp_dict['SNP'])
    print('... %d valid SNPs across populations ...' % n_snp)

    beta_dict = {}
    frq_dict = {}
    idx_dict = {}

    snp_set = set(snp_dict['SNP'])  # 加快 lookup

    for pp in range(n_pop):
        # 對齊該族群中的 SNP
        valid_idx = [i for i, snp in enumerate(sst_dict[pp]['SNP']) if snp in snp_set]
        beta_dict[pp] = np.array([sst_dict[pp]['BETA'][i] for i in valid_idx], ndmin=2).T
        frq_dict[pp] = np.array([sst_dict[pp]['FRQ'][i] for i in valid_idx], ndmin=2).T
        idx_dict[pp] = [snp_dict['SNP'].index(sst_dict[pp]['SNP'][i]) for i in valid_idx]

    return snp_dict, beta_dict, frq_dict, idx_dict

