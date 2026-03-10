from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import zipfile
import pandas as pd

FOCAL_GENERA = ["Ancylostoma", "Necator", "Oesophagostomum", "Trichostrongylus"]

@dataclass(frozen=True)
class Paths:
    share_zip: Path
    workdir: Path
    share_dir: Path
    meta_dir: Path
    tables_dir: Path
    kreports_dir: Path

def ensure_unpacked(share_zip: Path, workdir: Path) -> Paths:
    """Unpack share.zip into workdir if needed.

    Accepts either a zip that contains a top-level `share/` directory, or one that contains
    `kreports/`, `tables/`, etc. directly.
    """
    share_zip = Path(share_zip).resolve()
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    # If the user already unpacked, don't re-extract.
    # Otherwise extract.
    marker_candidates = [
        workdir / "share" / "tables" / "genus_bracken.tsv",
        workdir / "tables" / "genus_bracken.tsv",
    ]
    if not any(p.exists() for p in marker_candidates):
        with zipfile.ZipFile(share_zip, "r") as z:
            z.extractall(workdir)

    # Determine root
    root = (workdir / "share") if (workdir / "share").exists() else workdir

    return Paths(
        share_zip=share_zip,
        workdir=workdir,
        share_dir=root,
        meta_dir=root / "meta",
        tables_dir=root / "tables",
        kreports_dir=root / "kreports",
    )

def load_manifest(manifest_path: Path) -> pd.DataFrame:
    mf = Path(manifest_path)
    if not mf.exists():
        raise FileNotFoundError(f"Missing manifest: {mf}")
    df = pd.read_csv(mf, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    # normalize host
    if "host_group" in df.columns:
        host = df["host_group"].astype(str).str.strip().str.lower()
    elif "host" in df.columns:
        host = df["host"].astype(str).str.strip().str.lower()
    else:
        raise ValueError("Manifest missing host_group/host column")

    host = host.replace({"chimp": "chimpanzee", "chimps": "chimpanzee", "chimpanzees": "chimpanzee"})
    host = host.replace({"gorillas": "gorilla"})
    df["host"] = host

    df["SRR"] = df["SRR"].astype(str).str.strip()

    # order chimp then gorilla, SRR ascending within host
    host_order = pd.Categorical(df["host"], categories=["chimpanzee", "gorilla"], ordered=True)
    df = df.assign(_order=host_order).sort_values(["_order", "SRR"]).drop(columns=["_order"])
    return df

def default_manifest_path(p: Paths) -> Path:
    # preferred location inside share bundle
    candidate = p.meta_dir / "srr_manifest.txt"
    return candidate

def load_genus_table(p: Paths) -> pd.DataFrame:
    fp = p.tables_dir / "genus_bracken.tsv"
    if not fp.exists():
        raise FileNotFoundError(f"Missing genus table: {fp}")
    df = pd.read_csv(fp, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df["SRR"] = df["SRR"].astype(str).str.strip()
    df["name"] = df["name"].astype(str).str.strip()
    for c in ["kraken_assigned_reads", "added_reads", "new_est_reads", "fraction_total_reads"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def nonhost_reads_from_kreport(kreport_path: Path) -> int:
    """Estimate total non-host reads from a Kraken report.

    Uses: nonhost = unclassified_clade_reads + root_clade_reads
    """
    kreport_path = Path(kreport_path)
    if not kreport_path.exists():
        raise FileNotFoundError(f"Missing kreport: {kreport_path}")
    with kreport_path.open("r") as fh:
        line1 = fh.readline().strip().split("\t")
        line2 = fh.readline().strip().split("\t")
    unclassified = int(float(line1[1]))
    root = int(float(line2[1]))
    return unclassified + root
