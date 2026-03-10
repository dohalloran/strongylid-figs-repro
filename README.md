# strongylid-figs-repro

Reproducible scripts to regenerate **Figures** from the manuscript using the `share.zip` bundle (Kraken/Bracken outputs + SRR manifest).

## Input data expected

Put your `share.zip` at:

```
data/share.zip
```

The zip should contain:

```
share/
  meta/srr_manifest.txt
  tables/genus_bracken.tsv
  kreports/*.kreport
```

## Setup (recommended: conda)

```bash
conda create -n strongylid-py311 python=3.11 -y
conda activate strongylid-py311
pip install -r requirements.txt
```

## Run

```bash
python scripts/make_figs.py --share-zip data/share.zip --outdir outputs/figures
```

Outputs:

- `outputs/figures/FigureX.png`

## Notes

- Figures are written **PNG**.
- Figure 4 uses thresholds from 0 to 2500 reads per genus in steps of 50 (RA cutoff fixed at 0.1%).
- Figure 5 uses a strict filter of RA ≥ 0.1% and ≥100 reads per genus, and computes PRI_RPM as:
  `PRI_RPM = 1e6 * (sum genus new_est_reads passing filters) / (nonhost_reads)`
  where `nonhost_reads = unclassified_reads + root_classified_reads` from each `.kreport`.


### Manifest note

Some `share.zip` bundles may not include `meta/srr_manifest.txt`. If so, place the manifest next to the zip (e.g. `data/srr_manifest.txt`) and pass `--manifest data/srr_manifest.txt`.
