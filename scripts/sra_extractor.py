#!/usr/bin/env python3
"""
sra_extractor.py
---------------------------------
Unified extractor class for SRA metadata and sequence downloads.
Now automatically loads config.yaml from the project root.
"""

import os
import subprocess
import pandas as pd
import yaml
from tqdm import tqdm
from pathlib import Path
from Bio import Entrez


class SRAExtractor:
    def __init__(self, config_path="config.yaml"):
        """Initialize using values from config.yaml."""
        config_file = Path(config_path).resolve()
        if not config_file.exists():
            raise FileNotFoundError(f"❌ config.yaml not found at {config_file}")

        with open(config_file, "r") as f:
            cfg = yaml.safe_load(f)

        # Assign configuration
        self.project_name = cfg.get("project_name", "GAST")
        self.email = cfg.get("email")
        self.ena_preferred = cfg.get("ena_preferred", True)
        self.default_retmax = cfg.get("default_retmax", 20)

        # Directories
        paths = cfg.get("paths", {})
        self.data_dir = Path(paths.get("data_dir", "./data")).resolve()
        self.data_meta = Path(paths.get("metadata_dir", self.data_dir / "metadata")).resolve()
        self.data_seq = Path(paths.get("sequence_dir", self.data_dir / "sequences")).resolve()
        self.output_dir = Path(paths.get("output_dir", "./output")).resolve()

        # Create directories if missing
        for d in [self.data_meta, self.data_seq, self.output_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Set NCBI email for Entrez
        Entrez.email = self.email
        print(f"[INFO] Project '{self.project_name}' initialized")
        print(f"       Metadata path: {self.data_meta}")
        print(f"       Sequence path: {self.data_seq}")

    # ---------------------------
    # 1️⃣ Fetch metadata
    # ---------------------------
    def fetch_runinfo(self, organism, retmax=None):
        retmax = retmax or self.default_retmax
        print(f"[INFO] Fetching metadata for {organism} (max {retmax})")

        handle = Entrez.esearch(db="sra", term=organism, retmax=retmax)
        record = Entrez.read(handle)
        ids = record["IdList"]
        handle.close()

        if not ids:
            print("[WARNING] No records found.")
            return pd.DataFrame()

        summaries = []
        for sra_id in tqdm(ids, desc="Summaries"):
            handle = Entrez.esummary(db="sra", id=sra_id)
            summary = Entrez.read(handle)
            summaries.append(summary)
            handle.close()

        data = []
        for summary in summaries:
            for docsum in summary:
                data.append({
                    "SRA_ID": docsum.get("Id"),
                    "Title": docsum.get("Title"),
                    "Runs": docsum.get("Runs"),
                })

        df = pd.DataFrame(data)
        print(f"[INFO] Retrieved {len(df)} records.")
        return df

    # ---------------------------
    # 2️⃣ Save metadata
    # ---------------------------
    def save_metadata(self, df, filename="SraRunInfo.csv"):
        path = self.data_meta / filename
        df.to_csv(path, index=False)
        print(f"[INFO] Metadata saved to {path}")
        return path

    # ---------------------------
    # 3️⃣ Preview metadata
    # ---------------------------
    def preview_metadata(self, filename="SraRunInfo.csv", n=5):
        path = self.data_meta / filename
        if not path.exists():
            print(f"[ERROR] Metadata file not found: {path}")
            return
        df = pd.read_csv(path)
        print(df.head(n))
        return df.head(n)

    # ---------------------------
    # 4️⃣ Download SRA runs
    # ---------------------------
    def download_runs(self, run_ids):
        print("[INFO] Starting downloads...")

        for run_id in run_ids:
            out_dir = self.data_seq / run_id
            out_dir.mkdir(parents=True, exist_ok=True)

            fq1 = out_dir / f"{run_id}_1.fastq.gz"
            fq2 = out_dir / f"{run_id}_2.fastq.gz"

            if fq1.exists() and fq2.exists():
                print(f"[SKIP] {run_id} already exists.")
                continue

            print(f"🔽 Downloading {run_id}...")

            try:
                subprocess.run(["prefetch", run_id], check=True)
                subprocess.run(["fasterq-dump", "--split-files", run_id, "-O", str(out_dir)], check=True)
                subprocess.run(["gzip", str(out_dir / f"{run_id}_1.fastq")], check=True)
                subprocess.run(["gzip", str(out_dir / f"{run_id}_2.fastq")], check=True)
                print(f"[DONE] {run_id} successfully downloaded.")
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] Failed {run_id}: {e}")

        print("[INFO] All downloads completed.")
