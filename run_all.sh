#!/usr/bin/env bash
# End-to-end reproduction: downloads public data, runs the pipeline, writes
# figures, tables and the compiled manuscript.
#
#   bash run_all.sh                 full run (downloads ~3.2 GB)
#   bash run_all.sh --skip-download reuse data already in data/raw/
set -euo pipefail

RAW="data/raw"; RES="results"; FIG="figures"; TAB="tables"; MAN="manuscript"
mkdir -p "$RAW" "$RES" "$FIG" "$TAB" "$MAN"
SKIP=0; [[ "${1:-}" == "--skip-download" ]] && SKIP=1

GSE58885="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58885/matrix/GSE58885_series_matrix.txt.gz"
GSE58885_SOFT="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58885/soft/GSE58885_family.soft.gz"
GSE74193_MX="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74193/matrix/GSE74193_series_matrix.txt.gz"
GSE74193_PROC="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74193/suppl/GSE74193_GEO_procData.csv.gz"
MANIFEST="https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"
BRAINSPAN="https://www.brainspan.org/api/v2/well_known_file_download/267666525"

get () { [[ -s "$2" ]] || { echo "[get] $2"; curl -fsSL --retry 3 -o "$2" "$1"; }; }

if [[ $SKIP -eq 0 ]]; then
  echo "== downloading public data (about 3.2 GB) =="
  get "$GSE58885"      "$RAW/GSE58885_series_matrix.txt.gz"
  get "$MANIFEST"      "$RAW/manifest450k.csv.gz"
  get "$BRAINSPAN"     "$RAW/brainspan_genes_matrix_csv.zip"
  get "$GSE74193_MX"   "$RAW/GSE74193_series_matrix.txt.gz"
  get "$GSE74193_PROC" "$RAW/GSE74193_procData.csv.gz"
  get "$GSE58885_SOFT" "$RAW/GSE58885_family.soft.gz"
  [[ -d "$RAW/brainspan" ]] || { mkdir -p "$RAW/brainspan"; \
      unzip -oq "$RAW/brainspan_genes_matrix_csv.zip" -d "$RAW/brainspan"; }
fi

echo "== building sample metadata =="
python3 code/00_sample_metadata.py --soft "$RAW/GSE58885_family.soft.gz" \
    --matrix "$RAW/GSE58885_series_matrix.txt.gz" \
    --meta74193 "$RAW/GSE74193_series_matrix.txt.gz" --outdir "$RAW"

echo "== 01 per-probe statistics =="
python3 code/01_extract_methylation_stats.py --matrix "$RAW/GSE58885_series_matrix.txt.gz" \
    --meta "$RAW/gse58885_meta.csv" --out "$RES/methylation_probe_stats.npz"

echo "== 01b beta matrix =="
python3 code/01b_dump_beta_matrix.py --matrix "$RAW/GSE58885_series_matrix.txt.gz" \
    --stats "$RES/methylation_probe_stats.npz" --outdir "$RES"

echo "== 02 BrainSpan trajectories =="
python3 code/02_brainspan_machinery.py --dir "$RAW/brainspan" --outdir "$RES"

echo "== 03 probe annotation =="
python3 code/03_parse_annotation.py --manifest "$RAW/manifest450k.csv.gz" \
    --stats "$RES/methylation_probe_stats.npz" --out "$RES/probe_annotation.csv.gz"

echo "== 04 model-free profiles =="
python3 code/04_empirical_profiles.py --beta "$RES/beta_matrix.npy" \
    --stats "$RES/methylation_probe_stats.npz" --annot "$RES/probe_annotation.csv.gz" \
    --meta "$RES/sample_metadata.csv" --outdir "$RES" --nboot 300

echo "== 04b cumulative remodelling (primary measurement) =="
python3 code/04b_cumulative_remodelling.py --beta "$RES/beta_matrix.npy" \
    --stats "$RES/methylation_probe_stats.npz" --annot "$RES/probe_annotation.csv.gz" \
    --profiles "$RES/empirical_profiles.npz" --outdir "$RES" \
    --nboot 300 --nequal 120 --drop-bins 0 --tag _primary

echo "== 05 permanence kernel =="
python3 code/05_permanence_model.py --remodelling "$RES/cumulative_remodelling_primary.csv" \
    --outdir "$RES" --ndraw 20000

echo "== 06 independent cohort =="
python3 code/06_independent_cohort.py --proc "$RAW/GSE74193_procData.csv.gz" \
    --meta "$RAW/gse74193_meta.csv" --probes "$RES/methylation_probe_stats.npz" \
    --annot "$RES/probe_annotation.csv.gz" \
    --reference "$RES/cumulative_remodelling_primary.csv" --outdir "$RES"

echo "== 06b cohort comparability =="
python3 code/06b_cohort_comparability.py --proc "$RAW/GSE74193_procData.csv.gz" \
    --meta "$RAW/gse74193_meta.csv" --probes "$RES/methylation_probe_stats.npz" \
    --annot "$RES/probe_annotation.csv.gz" --beta "$RES/beta_matrix.npy" --outdir "$RES"

echo "== 07 identifiability and power =="
python3 code/07_deconvolution.py --kernel "$RES/permanence_kernel.csv" \
    --outdir "$RES" --n 60000 --nrep 20

echo "== 08 third-trimester scenarios =="
python3 code/08_third_trimester_scenarios.py \
    --remodelling "$RES/cumulative_remodelling_primary.csv" --outdir "$RES" --ndraw 3000

echo "== 09 figures =="; python3 code/09_figures.py --results "$RES" --outdir "$FIG"
echo "== 10 tables ==";  python3 code/10_tables.py  --results "$RES" --outdir "$TAB"
echo "== 11 references =="; python3 code/11_build_references.py --out "$MAN/references.md"
echo "== 12 manuscript ==";  python3 code/12_build_docx.py --md "$MAN/manuscript.md" \
    --refs "$MAN/references.md" --figures "$FIG" --tables "$TAB" --out "$MAN/manuscript.docx"

echo
echo "Done. Manuscript: $MAN/manuscript.docx"
