# A permanence window for prenatal epigenetic programming

Analysis code and derived results for the manuscript *"A permanence window for
prenatal epigenetic programming: why the critical week is biophysically narrow
and epidemiologically invisible."*

This repository reproduces every number, figure and table in the manuscript
from public data. Nothing in it is pre-computed by hand: the pipeline downloads
the source datasets, performs the analysis and regenerates the manuscript.

---

## What the analysis does

Environmental exposures in pregnancy are associated with neurodevelopmental
outcomes, but the gestational week at which an exposure becomes permanently
fixed in the fetal epigenome is unknown. This work decomposes vulnerability
into two factors that behave oppositely in time:

* **inducibility** `S(τ)` — the methylation machinery must be actively writing
  when the perturbation arrives;
* **retention** `R(τ)` — the displacement must survive to term rather than
  being diluted away by progenitor proliferation.

Their product is the **permanence kernel** `P(τ) = S(τ) R(τ)`, whose maximum
satisfies

```
k(τ*) = − d ln S / dτ |_{τ*}
```

where `k` is the local relaxation rate of the epigenetic state. Because `k`
falls as proliferation ceases while `S` falls as the machinery switches off,
the window is interior: it is neither the period of maximal plasticity nor the
end of gestation.

Headline results, all reproducible from this repository:

| Quantity | Value |
|---|---|
| Remodelling intensity peak | 15.4 post-conception weeks (bootstrap 95% CI 13.9–16.6) |
| Permanence window peak | 16.9 pcw (95% CI 15.8–18.0) |
| Half-maximum window | 13.2–19.0 pcw |
| Retention to term, 9 vs 17 pcw | 12.8% vs 45.1% |
| Machinery mRNA vs measured flux | *r* = 0.18 (no relationship) |
| Cross-dissection slope replication | observed/expected = 0.056 (does not replicate) |
| Cohort size for ±4-week localisation | ≈ 10⁶ pregnancies |

---

## Data sources

No data are redistributed here. All are public and are fetched by
`run_all.sh`.

| Dataset | Accession / URL | Size | Used for |
|---|---|---|---|
| Human fetal brain 450K methylation | [GSE58885](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58885) | 567 MB (gz) | discovery cohort, 179 donors, 3.3–26.3 pcw |
| Human DLPFC 450K methylation | [GSE74193](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74193) | 2.5 GB (gz) | independent replication, prenatal controls |
| Illumina 450K manifest | GEO platform [GPL13534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534) | 61 MB (gz) | probe annotation |
| BrainSpan Developmental Transcriptome | <https://www.brainspan.org> | 62 MB (zip) | proliferation and machinery trajectories |

Peak disk use is about 5 GB and peak memory about 2 GB. The full pipeline runs
in roughly three hours on one CPU core; the two largest downloads dominate the
wall-clock time.

---

## Quick start

```bash
git clone <this repository>
cd epigenetic-permanence-window

python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# downloads the data, runs everything, writes figures/, tables/ and the manuscript
bash run_all.sh
```

To re-run only the modelling once the data are downloaded and the methylation
matrix has been built:

```bash
bash run_all.sh --skip-download
```

---

## Pipeline

Each script is standalone and writes to `results/`. They are numbered in
execution order.

| Script | Purpose | Key outputs |
|---|---|---|
| `01_extract_methylation_stats.py` | one streaming pass over the GSE58885 series matrix; per-probe statistics and per-bin means and variances | `methylation_probe_stats.npz` |
| `01b_dump_beta_matrix.py` | materialises the probe × donor beta matrix as float32 | `beta_matrix.npy` |
| `02_brainspan_machinery.py` | neocortical proliferation, DNMT, TET, reader and one-carbon module trajectories | `brainspan_modules.csv`, `brainspan_genes.csv` |
| `03_parse_annotation.py` | compact per-probe annotation from the Illumina manifest | `probe_annotation.csv.gz` |
| `04_empirical_profiles.py` | model-free dispersion and intensity profiles, dominant temporal mode, donor bootstrap | `empirical_profiles.npz`, `dispersion_trajectory.csv` |
| `04b_cumulative_remodelling.py` | **main measurement**: pairwise methylomic distances, monotone embedding, Φ(t), equal-n control, leave-one-bin-out | `cumulative_remodelling_primary.csv/.npz` |
| `05_permanence_model.py` | inducibility–retention decomposition, Monte Carlo over priors, global sensitivity | `permanence_kernel.csv`, `permanence_summary.json` |
| `06_independent_cohort.py` | replication in GSE74193, composition-adjusted variant | `independent_cohort_*.csv/.json` |
| `06b_cohort_comparability.py` | attenuation-corrected test of whether the cohorts disagree or merely lack power | `cohort_comparability.json` |
| `07_deconvolution.py` | simulated epidemiology: trimester and distributed-lag analyses, statistical power | `deconvolution_*.csv/.json` |
| `08_third_trimester_scenarios.py` | sensitivity to the unmeasured third trimester | `third_trimester_*.csv/.json` |
| `09_figures.py` | all six figures as vector PDF and 600-dpi PNG | `figures/` |
| `10_tables.py` | all eight tables as CSV and Markdown | `tables/` |
| `11_build_references.py` | fetches every cited DOI from Crossref and formats a Vancouver list | `manuscript/references.md` |
| `12_build_docx.py` | compiles the formatted Word manuscript with figures and tables | `manuscript/manuscript.docx` |

### Methodological notes

* **Why pairwise distances rather than adjacent-week differences.** The
  sampling error of a bin mean is comparable to the one-week change, so
  differencing adjacent bins is noise dominated. All 66 pairwise distances are
  used jointly, with an unbiased noise correction, and a monotone
  one-dimensional embedding is fitted to them.
* **Why the earliest bin is excluded.** Samples at 3.3–7.0 pcw sit in the
  middle of the trajectory axis, not at its early end; at 23 days
  post-conception the tissue is neural tube, not neocortex. Excluding it raises
  the embedding fit from 0.874 to 0.962. Set `--drop-bins ""` to include it.
* **Why machinery mRNA is not used to predict flux.** It does not (*r* = 0.18),
  which is reported as a result. Inducibility is anchored on measured
  remodelling instead.
* **Why composition cannot be adjusted away.** In bulk fetal tissue the
  developmental trajectory *is* the leading principal component, so projecting
  it out removes the signal of interest. Running `04b` with `--npc 1` drops the
  embedding fit from 0.962 to 0.450; this is reported as a limitation, not used
  as a correction.

---

## Reproducibility

* Every stochastic step uses the fixed seed `20261027`.
* `requirements.txt` pins exact versions; `environment.yml` provides the same
  environment for conda.
* `11_build_references.py` queries Crossref for each DOI, so the bibliography
  cannot drift from the registry. Any DOI that fails to resolve is reported as
  an error rather than silently dropped.
* All derived result files in `results/`, `figures/` and `tables/` are produced
  by the scripts and are committed so that the manuscript can be checked
  without re-running the pipeline.

## Repository layout

```
code/        analysis scripts, numbered in execution order
results/     derived result files (CSV, JSON, NPZ)
figures/     publication figures (PDF vector + 600-dpi PNG)
tables/      manuscript tables (CSV + Markdown)
manuscript/  manuscript source, verified references, compiled DOCX
run_all.sh   end-to-end driver
```

Large intermediates (`beta_matrix.npy`, downloaded archives) are excluded by
`.gitignore` and are regenerated by the pipeline.

## Licence

Code is released under the MIT Licence (`LICENSE`). Derived result tables and
figures are released under CC BY 4.0. The underlying datasets remain subject to
the terms of their original repositories.

## Citing

See `CITATION.cff`. Please also cite the primary data sources: Spiers et al.
(2015) *Genome Research* for GSE58885, Jaffe et al. (2016) *Nature
Neuroscience* for GSE74193, and Kang et al. (2011) *Nature* for BrainSpan.
