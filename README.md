# Mini CRISPR–Cas9 NSC Aging Screens 

This is a re-creation and a compact analysis of the Ruetz et al. neural stem cell (NSC) CRISPR–Cas9 screens using the **processed per-sample count CSVs** from GEO (**GSE189251**).  
To keep the project light and reproducible, I selected **~¼ of the sample CSV files per condition** (¼ of counts), pooled **Young (Y)** vs **Old (O)** for each condition (**RA**, **Kp**, **Q**), and computed CPM → log2FC → permutation p-values → BH-FDR.  

All code lives under `cas9_analysis/workflow/`. 

**Project structure:**
```bash
crispr_nsc_aging/
├── README.md
└── cas9_analysis/
    ├── data/
    │   ├── raw/
    │   └── processed/
    ├── questions.yaml
    ├── answers.yaml
    ├── metadata.yaml      
    └── workflow/
         ├── run_workflow.py
         │  
         │
         └── outputs/
            ├── RA_Y_vs_O_pooled/
            ├── Kp_Y_vs_O_pooled/
            ├── Q_Y_vs_O_pooled/
            └── submission_artifacts_ALL_conditions_quarterFILES.zip
```

## 1) Data sources 

- **GEO accession:** `GSE189251`  
- **Inputs:** processed per-sample counts (`*_counts.csv.gz`) bundled in  
  `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189251/suppl/GSE189251_RAW.tar`
- **Conditions analyzed:** in-vitro SC: **RA**, **Kp**, **Q**  
- **Groups:** **Y** (Young) vs **O** (Old)

> Started from processed **CSV** counts (no re-alignment).  
> Reduced **the number of files** (¼ per condition)

## 2) Downloading the processed data

File location: `cas9_analysis/data/`.

```bash
# from repo root: crispr_nsc_aging/
mkdir -p cas9_analysis/data/raw cas9_analysis/data/processed

# Download processed counts tarball (skip if you already have it)
wget -O cas9_analysis/data/raw/GSE189251_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189251/suppl/GSE189251_RAW.tar"
```

## 3) Pre-processing / subsampling
Miniaturize by selecting ~¼ of the sample CSV files per condition (Y and O, minimum 1 each), rather than thinning counts:
1) Extracted all *_counts.csv.gz files to cas9_analysis/data/processed/ ;
2) Kept SC (in-vitro) samples only; dropped lib_counts and in-vivo/WT ;
3) For each condition (RA/Kp/Q), kept ceil(0.25 × N) files for Y and for O ;
4) Pooled the counts within Y and within O (sum per sgRNA element) ;
5) Computed CPM, element log2FC (O/Y), element p via label shuffling, and gene-level BH-FDR.

## 4) Workflow
Outputs go to cas9_analysis/workflow/outputs/.

### Run the analysis
File location: cas9_analysis/workflow/run_workflow.py

**What it does, per condition (RA/Kp/Q):**
- Discovers SC files (ignores lib_counts and in-vivo) ;
- Keeps ¼ of files for Y and ¼ of files for O (min 1 each) ;
- Pools counts within Y and within O ;
- Computes CPM; element-level log2FΩ ;
- Permutation testing → empirical p ;
- Gene-level BH-FDR ;
- Saves per-condition artifacts + a master ZIP

```bash
python3 cas9_analysis/workflow/run_workflow.py \
  --project cas9_analysis \
  --extract \
  --fraction 0.25 \
  --permutations 2000
```
**Outputs (per condition):** cas9_analysis/workflow/outputs/<COND>_Y_vs_O_pooled/

**Files:** 1) elements_log2fc.csv (2) genes_aggregated.csv (3) volcano_gene.png (4) FILES_USED.txt (5) submission_artifacts_<COND>_pooled_quarterFILES.zip
**Plus a master archive:** cas9_analysis/workflow/outputs/submission_artifacts_ALL_conditions_quarterFILES.zip

## 5) Reproduce locally

```bash
python3 -m venv venv && source venv/bin/activate
pip install numpy pandas matplotlib pyyaml

python cas9_analysis/workflow/run_workflow.py \
  --project cas9_analysis \
  --extract \
  --fraction 0.25
```

**Notes:**

1) Uses ¼ of files per condition (Y/O) — not count thinning.
2) Outputs are small (MBs), safe for sharing.
3) RNG seed fixed for reproducibility.
