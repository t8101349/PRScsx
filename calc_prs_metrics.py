#!/usr/bin/env python
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve

# 參數
parser = argparse.ArgumentParser(description="計算 PRS AUC 與最佳 F1 score，並繪圖")
parser.add_argument("--score", required=True, help="plink --score 輸出檔案")
parser.add_argument("--pheno", required=True, help="phenotype 檔 (IID, PHENO)")
parser.add_argument("--pheno_name", required=True, help="phenotype name (e.g., ISCOLONCANCER)")
parser.add_argument("--out_prefix", default="prs_eval", help="輸出檔名前綴")
args = parser.parse_args()

# 讀取資料
pheno_name = args.pheno_name
score_df = pd.read_csv(args.score, delim_whitespace=True)
score_df = score_df.iloc[:, [1, -1]]  # 取 IID, PRS
score_df.columns = ["IID", "PRS"]
pheno_df = pd.read_csv(args.pheno, sep="\t", usecols=[0, 1], names=["IID", pheno_name], header=0)


# 對齊樣本（以 IID 為 key）
df = pd.merge(pheno_df, score_df, on="IID", how="inner")

# y = phenotype (假設 ISCOLONCANCER 為 0/1)
y_true = (df[pheno_name]-1).astype(int)
y_score = df["PRS"]

print(np.unique(y_true, return_counts=True))
print(np.min(y_score), np.max(y_score))


# 計算 AUC
auc = roc_auc_score(y_true, y_score)

# Precision-Recall & F1
precision, recall, thresholds = precision_recall_curve(y_true, y_score)
f1_scores = np.divide(2 * precision * recall, precision + recall, out=np.zeros_like(precision), where=(precision + recall) != 0)
best_idx = f1_scores.argmax()
best_threshold = thresholds[best_idx] if best_idx < len(thresholds) else thresholds[-1]
best_f1 = f1_scores[best_idx]

# 儲存數值結果
with open(f"{args.pheno_name}{args.out_prefix}.txt", "w") as f:
    f.write(f"AUC\t{auc:.4f}\n")
    f.write(f"Best_threshold\t{best_threshold:.6f}\n")
    f.write(f"Best_F1\t{best_f1:.4f}\n")

# 繪製 ROC 曲線
sns.set(style="whitegrid")
fpr, tpr, _ = roc_curve(y_true, y_score)

plt.figure(figsize=(6, 6))
sns.lineplot(x=fpr, y=tpr, label=f"AUC = {auc:.3f}")
sns.lineplot(x=[0, 1], y=[0, 1], linestyle="--", color="gray")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend(loc="lower right")
plt.savefig(f"{args.pheno_name}{args.out_prefix}_roc.png", dpi=300)
plt.close()

# 繪製 Precision-Recall 曲線
plt.figure(figsize=(6, 6))
sns.lineplot(x=recall, y=precision, label=f"Best F1 = {best_f1:.3f}")
plt.scatter(recall[best_idx], precision[best_idx], color="red", zorder=10, label="Best Threshold")
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall Curve")
plt.legend(loc="upper right")
plt.savefig(f"{args.pheno_name}{args.out_prefix}_pr.png", dpi=300)
plt.close()

print(f"[Result] AUC = {auc:.4f}")
print(f"[Result] Best threshold = {best_threshold:.6f}")
print(f"[Result] Best F1 score = {best_f1:.4f}")
print(f"已輸出：\n  - 指標檔案 {args.pheno_name}{args.out_prefix}.txt\n  - 圖檔 {args.pheno_name}{args.out_prefix}_roc.png, {args.pheno_name}{args.out_prefix}_pr.png")

"""
demo
python /home/Weber/Pipeline/PRS/calc_prs_metrics.py \
				--score PRS_EAS/COLONCANCER_PRS.validate.score\
				--pheno ISCOLONCANCERpheno.matchit.txt \
				--pheno_name coloncancer_valid
"""