# ===============================================================
# Describe: calculated PRS by target dataset using PRSCS
# Author: Rogen
# Date: 2022.12.02
# Parameters: 
# Note: It must use python3 and PLINL 1.X in environment, and
# 		PRSce recepted PLINK 1.X summary statistics 
#		format and input bfile only.
# Update:
# 2023.07.19: Completed plink --score function. 
# 2024.02.02: Adding --target_list / --val_list arguments.
# ===============================================================
#!/bin/bash

# === 初始設定 ===
prscsx="/home/Weber/PRScsx"
refld="/home/Weber/PRScsx/ldblk_1kg_eas"
BETA_STD=false
pop="EAS"
threads=15

# === 參數解析 ===
re='^(--help|-h)$'
if [[ $1 =~ $re ]]; then
    Help
else
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -r|--ref_ld_dir) ref_ld_dir="$2"; shift 2;;
            -b|--bfile) bfile="$2"; shift 2;;
            -g|--gwas) gwas="$2"; shift 2;;
            -v|--plink_version) plink_version="$2"; shift 2;;
            -n|--n_gwas) n_gwas="$2"; shift 2;;
            -c|--chr) chr="$2"; shift 2;;
            -p|--phi) phi="$2"; shift 2;;
            -m|--MCMC_ITER) MCMC_ITER="$2"; shift 2;;
            -u|--MCMC_BURNIN) MCMC_BURNIN="$2"; shift 2;;
            -z|--BETA_STD) BETA_STD=true; shift 1;;
            -t|--threads) threads="$2"; shift 2;;
            -o|--out) out="$2"; shift 2;;
            --target_list) target_list="$2"; shift 2;;
            --val_list) val_list="$2"; shift 2;;
            --pop) pop="$2"; shift 2;;
            *) echo "unknown option: $1" >&2; exit 1;;
        esac
    done

    # === 檢查必要參數 ===
    if [[ -z "$bfile" || -z "$gwas" || -z "$n_gwas" || -z "$out" ]]; then
        echo "❌ 必要參數缺失，請確認 --bfile, --gwas, --n_gwas, --out 是否有指定。" >&2
        exit 1
    fi
    if [[ -z "$pop" ]]; then
        echo "Alarm: --pop cannot be empty"
        pop="EAS"
    fi

    echo "ref_dir:${ref_ld_dir}"
    echo "${pop}"

    # === plink version 判斷 ===
    if [ -z "$plink_version" ] || [ "$plink_version" = "v2" ]; then
        plink_version="plink_v2"
    elif [ "$plink_version" = "v1" ]; then
        plink_version="plink_v1"
    elif [ "$plink_version" = "SAIGE" ]; then
        plink_version="SAIGE"
    else
        echo "The \"--plink_version\" input is aberrant."
        exit 1
    fi

    tmpdir="./temp/prscs_temp_$RANDOM"
    mkdir -p "$tmpdir"

    output_dir=$(dirname "$out")
    outname=$(basename "$out")
    mkdir -p "$output_dir"

    eval "python /home/Weber/Pipeline/PRS/PRSCS_sumstat_adjust.py $gwas $plink_version $out"
    awk 'NR==1 || (length($3)==1 && length($4)==1) { $1=""; sub(/^ /, ""); print }' ${out}.sumstats.prscs.txt > ${out}.sumstats.txt

    export MKL_NUM_THREADS=$threads
    export NUMEXPR_NUM_THREADS=$threads
    export OMP_NUM_THREADS=$threads

    # === chr 處理 ===
    if [ -z "$chr" ] || [ "$chr" == "all" ]; then
        list=$(seq 1 22)
    elif [[ "$chr" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
        list=$(echo "$chr" | tr "," " ")
    elif [[ "$chr" =~ ^[0-9]+\.\.[0-9]+$ ]]; then
        start=${chr%%..*}
        end=${chr##*..}
        list=$(seq "$start" "$end")
    else
        echo "The \"--chr\" input aberrant."
        exit 1
    fi

    # === 平行 PRScsx ===
    run_prscsx () {
        chr=$1
        phi=$2
        n_gwas=$3
        bfile=$4
        outdir=$5
        outname=$6
        refld=$7
        iter=$8
        burnin=$9
        std=${10}
        prscsx=${11}
        pop=${12}
    
        mkdir -p "$outdir"
    
        cmd="python $prscsx/PRScsx.py \
            --bim_prefix=$bfile \
            --sst_file=${outdir}/${outname}.sumstats.prscs.txt \
            --n_gwas=$n_gwas \
            --chrom=$chr \
            --ref_dir=$refld \
            --out_dir=$outdir \
            --out_name=${outname}_chr${chr} \
            --pop=$pop \
            --write_pst TRUE \
            --seed=42"
    
        if [ ! -z "$iter" ]; then
            cmd="$cmd --n_iter=$iter"
        fi
        if [ ! -z "$burnin" ]; then
            cmd="$cmd --n_burnin=$burnin"
        fi
        if [ "$std" = true ]; then
            cmd="$cmd --beta_std"
        fi
        if [ ! -z "$phi" ]; then
            cmd="$cmd --phi=$phi"
        fi
    
        echo "Running: $cmd"
        eval $cmd 2>&1 | tee "${outdir}/${outname}_chr${chr}.log"
    }
    export -f run_prscsx

    echo "POP is: $pop"

    parallel -j "3" run_prscsx \
    ::: $list ::: "$phi" ::: "$n_gwas" ::: "$bfile" ::: "$output_dir" ::: "$outname" ::: "$refld" ::: "$MCMC_ITER" ::: "$MCMC_BURNIN" ::: "$BETA_STD" ::: "$prscsx" ::: "$pop"


    ## === 合併結果並製作 scorefile ===
    cat ${output_dir}/${outname}_chr*_${pop}_pst_eff_a1_b0.5_phiauto_chr*.txt > ${output_dir}/${outname}_merged.txt  
    awk 'NR>1 {print $2, $4, $6}' ${output_dir}/${outname}_merged.txt > ${output_dir}/${outname}.scorefile.txt

    # === 準備 keep 檔案 ===
    awk 'NR>1{print $1}' "$target_list" > "$tmpdir/target.iid"
    awk 'NR>1{print $1}' "$val_list" > "$tmpdir/val.iid"
    cat "$tmpdir/target.iid" "$tmpdir/val.iid" | sort | uniq > "$tmpdir/merged.iid"
    awk 'BEGIN{OFS="\t"} {print 0, $1}' "$tmpdir/merged.iid" > "$tmpdir/merged.keep"

    # === 執行 plink --score ===
    plink --bfile "$bfile" \
          --score "${output_dir}/${outname}_merged.txt" 2 4 6 sum \
          --threads "$threads" \
          --keep "$tmpdir/merged.keep" \
          --out "${output_dir}/${outname}"

    # === 拆分 target / validate ===
    awk 'NR==FNR {a[$1]; next} ($2 in a)' "$target_list" "${output_dir}/${outname}.profile" > "${output_dir}/${outname}.target.score"
    awk 'NR==FNR {a[$1]; next} ($2 in a)' "$val_list" "${output_dir}/${outname}.profile" > "${output_dir}/${outname}.validate.score"

    echo "PRS-CSX completed at $(date)"

    # read -p "phenotype name is: " phename
    # echo "phenotype name is $phename"
    upper_phename=$(echo "${outname%%_*}" | tr '[:lower:]' '[:upper:]')
    echo "${upper_phename} plotting"

    lower_phename=$(echo "$upper_phename" | tr '[:upper:]' '[:lower:]')

    # 作圖  PRS_Distribution_Plot_csx.py
    python /home/Weber/Pipeline/PRS/PRS_Distribution_Plot_csx.py \
    --prs "${output_dir}/${outname}.validate.score" \
    --pheno IS${upper_phename}pheno.matchit.txt \
    --phename "${lower_phename}" \
    --normalize z_std \
    --out ${lower_phename}_validate

     echo "PRS Distribution Plot completed"
    
fi


: << Demo
bash /home/Weber/Pipeline/PRS/PRSCSX_for_target.sh \
				-b /SNParray/SourceShare/20240321_50w_Imputation/step12-pgen2bed/Axiom_imputed_r2.MAF \
				-g /home/Weber/Cancer/coloncancer/20250709/GWAS2/ISCOLONCANCER_TPMI_imputed_adjGWAS.glm.logistic \
				-n 59840 \
				-c all \
				-o PRS_EAS/COLONCANCER_PRS \
                   --target_list /home/Weber/Cancer/coloncancer/20250709/coloncancer_split.ind_target_set.tsv \
                   --val_list /home/Weber/Cancer/coloncancer/20250709/coloncancer_split.ind_validate_set.tsv
Demo