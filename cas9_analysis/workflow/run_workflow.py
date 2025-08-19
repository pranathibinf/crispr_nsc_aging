#!/usr/bin/env python3
"""
run_workflow.py — Mini CRISPR NSC analysis on processed counts (GSE189251)

Expected layout (in the *cas9_analysis* folder):
  data/raw/         : GSE189251_RAW.tar lives here
  data/processed/   : extracted *_counts.csv.gz
  workflow/outputs/<COND>_Y_vs_O_pooled/ : results

Pipeline per condition (RA, Kp, Q):
  • pick ~1/4 of SC files for Y and for O (min 1 each)
  • pool counts within Y and within O (sum per sgRNA element)
  • CPM, element log2FC (O/Y), permutation p-values
  • gene-level median p, BH-FDR
  • save CSVs, volcano plot, FILES_USED.txt, per-condition zip + master zip
"""
import argparse, os, tarfile, glob, re, math, zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def smart_read_counts(path: str) -> pd.DataFrame:
    """Read csv/csv.gz → two columns: element, count (int)."""
    for sep in [",", "\t", ";"]:
        for header in [None, 0]:
            try:
                df = pd.read_csv(path, sep=sep, header=header, compression="infer", engine="python")
                if header is None:
                    df = df.rename(columns={0: "element", 1: "count"})
                cols = [str(c).lower() for c in df.columns]
                if "element" in cols and "count" in cols:
                    df = df[[df.columns[cols.index("element")], df.columns[cols.index("count")]]]
                    df.columns = ["element", "count"]
                else:
                    df = df.iloc[:, :2]; df.columns = ["element", "count"]
                df["element"] = df["element"].astype(str)
                df["count"]   = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(int)
                return df
            except Exception:
                pass
    raise RuntimeError(f"Could not parse {path}")

def to_cpm(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    total = int(out["count"].sum())
    out["cpm"] = (out["count"] / max(total, 1)) * 1_000_000.0
    return out

def fdr_bh(pvals):
    p = np.asarray(pvals, float); n = len(p)
    order = np.argsort(p); ranked = p[order]
    q = ranked * n / (np.arange(1, n+1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q); out[order] = q
    return out

def pool_group(paths):
    pooled = None
    for p in paths:
        df = smart_read_counts(p)
        if pooled is None:
            pooled = df.copy()
        else:
            pooled = pooled.merge(df, on="element", how="outer", suffixes=("","_x"))
            pooled["count"] = pooled["count"].fillna(0).astype(int) + pooled["count_x"].fillna(0).astype(int)
            pooled = pooled[["element","count"]]
    if pooled is None:
        pooled = pd.DataFrame({"element":[], "count":[]})
    return pooled

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def run(args):
    # --project should be the cas9_analysis folder
    PROJECT  = args.project
    RAW_DIR  = os.path.join(PROJECT, "data", "raw")
    PROC_DIR = os.path.join(PROJECT, "data", "processed")
    OUT_ROOT = os.path.join(PROJECT, "workflow", "outputs")  # <-- fixed

    for d in (RAW_DIR, PROC_DIR, OUT_ROOT): ensure_dir(d)

    tar_path = os.path.join(RAW_DIR, "GSE189251_RAW.tar")
    if args.download and not os.path.exists(tar_path):
        import urllib.request
        print("[download] fetching GSE189251_RAW.tar ...")
        urllib.request.urlretrieve(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189251/suppl/GSE189251_RAW.tar",
            tar_path
        )

    if args.extract:
        print("[extract] extracting tar → data/processed/")
        with tarfile.open(tar_path, "r") as tar:
            tar.extractall(PROC_DIR)

    # SC files only (ignore library/in vivo)
    all_sc = [p for p in glob.glob(os.path.join(PROC_DIR, "GSM*_SC*.csv.gz"))
              if "lib_counts" not in p]

    def parse_labels(path):
        b = os.path.basename(path)
        cond  = re.search(r"_(RA|Kp|Q)_", b)
        group = re.search(r"_(Y|O)_", b)
        return (cond.group(1) if cond else None,
                group.group(1) if group else None)

    cond_groups = {}
    for p in all_sc:
        cond, grp = parse_labels(p)
        if cond in ("RA","Kp","Q") and grp in ("Y","O"):
            cond_groups.setdefault(cond, {"Y":[], "O":[]})[grp].append(p)

    print("[discover] SC files per condition:")
    for c in ("RA","Kp","Q"):
        y = len(cond_groups.get(c,{}).get("Y",[]))
        o = len(cond_groups.get(c,{}).get("O",[]))
        print(f"  {c}: Y={y} O={o}")

    # keep ~1/4 of files (min 1)
    subset = {}
    for c, groups in cond_groups.items():
        Y = sorted(groups["Y"]); O = sorted(groups["O"])
        y_keep = max(1, math.ceil(len(Y) * args.fraction))
        o_keep = max(1, math.ceil(len(O) * args.fraction))
        subset[c] = {"Y": Y[:y_keep], "O": O[:o_keep]}
        print(f"[subset] {c}: keep {y_keep}/{len(Y)} (Y), {o_keep}/{len(O)} (O)")

    rng = np.random.default_rng(args.seed)
    created = []
    for c, groups in subset.items():
        out_dir = os.path.join(OUT_ROOT, f"{c}_Y_vs_O_pooled")
        ensure_dir(out_dir)

        # record inputs used
        with open(os.path.join(out_dir, "FILES_USED.txt"), "w") as fh:
            for grp in ("Y","O"):
                fh.write(f"[{grp}]\n")
                for p in groups[grp]: fh.write(os.path.basename(p)+"\n")
                fh.write("\n")

        # pool
        dfY = pool_group(groups["Y"])
        dfO = pool_group(groups["O"])
        if dfY.empty or dfO.empty:
            print(f"[warn] {c}: missing pooled data; skipping")
            continue

        # CPM + merge
        dfY_cpm = to_cpm(dfY).rename(columns={"count":"count_Y","cpm":"cpm_Y"})
        dfO_cpm = to_cpm(dfO).rename(columns={"count":"count_O","cpm":"cpm_O"})
        m = dfY_cpm.merge(dfO_cpm, on="element", how="inner")

        # element stats
        m["log2FC"] = np.log2((m["cpm_O"] + 1.0) / (m["cpm_Y"] + 1.0))
        m["gene"]   = m["element"].str.split("_", n=1).str[0]

        # permutations (label shuffling at element level)
        cY = m["cpm_Y"].to_numpy(); cO = m["cpm_O"].to_numpy()
        obs = np.log2((cO + 1.0)/(cY + 1.0))
        abs_obs = np.abs(obs)
        null_vals = []
        for _ in range(args.permutations):
            flip = rng.random(len(cY)) < 0.5
            Operm = np.where(flip, cY, cO)
            Yperm = np.where(flip, cO, cY)
            l2 = np.log2((Operm + 1.0)/(Yperm + 1.0))
            null_vals.append(np.abs(l2))
        null = np.vstack(null_vals)
        m["p_element"] = (null >= abs_obs).mean(axis=0)

        # gene aggregate
        agg = (m.groupby("gene", as_index=False)
                 .agg(mean_log2FC=("log2FC","mean"),
                      n_elements=("element","count"),
                      p_gene=("p_element","median")))
        agg["q_gene"] = fdr_bh(agg["p_gene"].values)

        # save artifacts
        m.to_csv(os.path.join(out_dir, "elements_log2fc.csv"), index=False)
        agg.to_csv(os.path.join(out_dir, "genes_aggregated.csv"), index=False)

        plt.figure()
        plt.scatter(agg["mean_log2FC"], -np.log10(np.clip(agg["p_gene"], 1e-12, 1.0)))
        plt.xlabel("Mean log2FC (gene)")
        plt.ylabel("-log10 p-value (gene)")
        plt.title(f"Volcano — {c} (O vs Y, pooled ~1/4 files)")
        plt.savefig(os.path.join(out_dir, "volcano_gene.png"), bbox_inches="tight")
        plt.close()

        zpath = os.path.join(out_dir, f"submission_artifacts_{c}_pooled_quarterFILES.zip")
        if os.path.exists(zpath): os.remove(zpath)
        with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
            for fn in ("elements_log2fc.csv","genes_aggregated.csv","volcano_gene.png","FILES_USED.txt"):
                z.write(os.path.join(out_dir, fn), arcname=fn)
        created.append(zpath)

    # master zip
    master = os.path.join(OUT_ROOT, "submission_artifacts_ALL_conditions_quarterFILES.zip")
    if os.path.exists(master): os.remove(master)
    with zipfile.ZipFile(master, "w", zipfile.ZIP_DEFLATED) as z:
        for zp in created:
            z.write(zp, arcname=os.path.basename(zp))
    print("[done] master zip:", master)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="CRISPR NSC mini-pipeline on processed counts. "
                    "--project should point to the cas9_analysis folder."
    )
    ap.add_argument("--project", default="cas9_analysis",
                    help="path to cas9_analysis (contains data/ and workflow/)")
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--extract", action="store_true")
    ap.add_argument("--fraction", type=float, default=0.25)
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    run(args)
