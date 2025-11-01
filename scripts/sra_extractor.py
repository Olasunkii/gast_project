#!/usr/bin/env python3
"""
sra_extractor.py
Interactive SRA metadata fetcher & sequence downloader.
Works in both Jupyter Notebook and CLI.
"""

import os
import time
import csv
import subprocess
from pathlib import Path
from typing import Optional, List
import pandas as pd
from Bio import Entrez
from tqdm import tqdm
import requests
import io


class SRAExtractor:
    def __init__(self,
                 project_root: Optional[str] = None,
                 email: Optional[str] = None,
                 default_retmax: int = 20,
                 sleep_between_requests: float = 0.34):
        """
        project_root: root folder for the project (creates data/metadata, data/sequences)
        email: NCBI Entrez email (required)
        """
        # ✅ Auto-detect project root two levels up (scripts/ -> project root)
        if project_root is None:
            # ensures correct root even when called from notebooks/
            self.project_root = Path(__file__).resolve().parent.parent
        else:
            self.project_root = Path(project_root).resolve()

        # Create necessary directories
        self.data_meta = self.project_root / "data" / "metadata"
        self.data_seq = self.project_root / "data" / "sequences"
        self.data_meta.mkdir(parents=True, exist_ok=True)
        self.data_seq.mkdir(parents=True, exist_ok=True)

        self.default_retmax = default_retmax
        self.sleep_between_requests = sleep_between_requests

        # ✅ Ensure email setup
        if email:
            self.email = email
        elif getattr(Entrez, "email", None):
            self.email = Entrez.email
        else:
            self.email = input("Enter your email (required for NCBI Entrez): ").strip()

        Entrez.email = self.email

        print(f"[INFO] Project initialized at {self.project_root}")
        print(f"Metadata folder: {self.data_meta}")
        print(f"Sequences folder: {self.data_seq}")

    #######################################################
    # Metadata fetching
    #######################################################

    def search_sra(self, organism: str, retmax: Optional[int] = None) -> List[str]:
        retmax = retmax or self.default_retmax
        term = f'{organism}[Organism]'
        handle = Entrez.esearch(db="sra", term=term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        return record.get("IdList", [])

    def fetch_runinfo(self, organism_or_idlist: str, retmax: Optional[int] = None) -> pd.DataFrame:
        """
        Fetch SRA run metadata for a given organism name or list of IDs.
        Saves CSV metadata file to data/metadata.
        """
        retmax = retmax or self.default_retmax
        print(f"[INFO] Fetching metadata for '{organism_or_idlist}' (retmax={retmax})")

        # Step 1: Search for organism in SRA
        try:
            handle = Entrez.esearch(db="sra", term=organism_or_idlist, retmax=retmax)
            record = Entrez.read(handle)
            handle.close()
            ids = record["IdList"]
            print(f"[INFO] Found {len(ids)} results for '{organism_or_idlist}'")
        except Exception as e:
            print(f"[ERROR] Failed to fetch SRA IDs: {e}")
            return pd.DataFrame()

        if not ids:
            print("[WARNING] No results found.")
            return pd.DataFrame()

        # Step 2: Fetch the runinfo table
        try:
            with Entrez.efetch(db="sra", id=",".join(ids), rettype="runinfo", retmode="text") as handle:
                raw_data = handle.read()

            if isinstance(raw_data, bytes):
                raw_data = raw_data.decode("utf-8", errors="replace")

            df = pd.read_csv(io.StringIO(str(raw_data)))
            self.last_metadata = df.copy()

            print(f"[INFO] Retrieved {len(df)} metadata records.")
            return df

        except Exception as e:
            print(f"[ERROR] Failed to fetch runinfo: {e}")
            return pd.DataFrame()

    def save_metadata(self, df: pd.DataFrame, filename: Optional[str] = None) -> Path:
        fname = filename or f"SraRunInfo_{int(time.time())}.csv"
        out = self.data_meta / fname
        df.to_csv(out, index=False)
        print(f"[INFO] Metadata saved to {out}")
        return out

    def preview_metadata(self, csvpath: Optional[str] = None, n: int = 10):
        if csvpath:
            p = Path(csvpath)
            if not p.is_absolute():
                p = self.data_meta / p
            df = pd.read_csv(p, dtype=str)
        else:
            df = getattr(self, "last_metadata", None)
            if df is None:
                raise ValueError("No metadata loaded and no csvpath provided.")

        display_cols = [c for c in ["Run", "ScientificName", "LibraryLayout", "Platform", "size_MB", "bases"] if c in df.columns]
        if not display_cols:
            display_cols = df.columns.tolist()[:6]

        print(df[display_cols].head(n).to_string(index=False))

    #######################################################
    # Sequence download
    #######################################################
    def _get_ena_fastq_links(self, run_id: str) -> list:
        """Fetch FASTQ download URLs for a given SRR run from ENA."""
        url = (
            "https://www.ebi.ac.uk/ena/portal/api/filereport"
            f"?accession={run_id}&result=read_run&fields=fastq_ftp&format=tsv"
        )
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            lines = r.text.strip().split("\n")
            if len(lines) < 2:
                print(f"[WARN] No FASTQ entries found for {run_id}.")
                return []
            fastq_field = lines[1].split("\t")[-1].strip()
            if not fastq_field:
                return []
            fastq_links = fastq_field.split(";")
            fastq_links = [f"https://{link.strip()}" if not link.startswith("http") else link
                           for link in fastq_links if link.strip()]
            return fastq_links
        except Exception as e:
            print(f"[ERROR] Could not retrieve ENA links for {run_id}: {e}")
            return []

    def _download_file(self, url, outpath):
        outpath = Path(outpath)
        tmp_path = outpath.with_suffix(outpath.suffix + ".part")
        headers = {}
        if tmp_path.exists():
            existing_size = tmp_path.stat().st_size
            headers["Range"] = f"bytes={existing_size}-"
        else:
            existing_size = 0

        try:
            response = requests.get(url, stream=True, headers=headers, timeout=60)
            response.raise_for_status()
        except requests.RequestException as e:
            print(f"[ERROR] Failed to start download from {url}: {e}")
            return None

        total_size = int(response.headers.get("content-length", 0)) + existing_size
        mode = "ab" if existing_size > 0 else "wb"

        with open(tmp_path, mode) as f, tqdm(
            total=total_size, initial=existing_size, unit="B", unit_scale=True,
            desc=outpath.name, leave=True
        ) as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))

        tmp_path.rename(outpath)
        return outpath

    def download_runs(self, run_ids, prefer_ena=True):
        """Download FASTQ files for given SRR IDs (from ENA or NCBI)."""
        import subprocess
        from datetime import datetime

        if isinstance(run_ids, str):
            run_ids = [run_ids]

        log_file = self.data_seq / "download_log.csv"
        log_exists = log_file.exists()

        sra_temp_dir = self.data_seq / "sra_temp"
        sra_temp_dir.mkdir(parents=True, exist_ok=True)

        with open(log_file, "a", newline="") as log:
            writer = csv.writer(log)
            if not log_exists:
                writer.writerow(["Run", "Status", "Files", "Timestamp"])

            for run in run_ids:
                run_dir = self.data_seq / run
                run_dir.mkdir(exist_ok=True)
                downloaded_files = []

                # Skip if already downloaded
                if any(run_dir.glob("*.fastq*")):
                    print(f"✅ {run} already downloaded — skipping.")
                    writer.writerow([run, "exists", "existing", datetime.now().strftime("%Y-%m-%d %H:%M")])
                    continue

                print(f"\n🔽 Downloading {run}...")

                urls = self._get_ena_fastq_links(run) if prefer_ena else []

                # ==============================
                # 🔹 ENA FASTQ DOWNLOAD PATHWAY
                # ==============================
                if urls:
                    print(f"Found {len(urls)} FASTQ file(s) on ENA.")
                    for url in urls:
                        fname = Path(url).name
                        outpath = run_dir / fname
                        if not outpath.exists():
                            self._download_file(url, outpath)
                        downloaded_files.append(str(outpath))
                    writer.writerow([run, "downloaded_ENA", ";".join(downloaded_files),
                                    datetime.now().strftime("%Y-%m-%d %H:%M")])

                # ==============================
                # 🔹 NCBI FASTQ DOWNLOAD PATHWAY
                # ==============================
                else:
                    print(f"⚠️ No ENA FASTQ found for {run}. Attempting NCBI prefetch + fasterq-dump...")
                    try:
                        # Step 1: Prefetch the SRA data
                        subprocess.run(["prefetch", run, "--output-directory", str(sra_temp_dir)], check=True)

                        # Step 2: Detect correct SRA file path
                        sra_subdir = sra_temp_dir / run
                        sra_file = sra_subdir / f"{run}.sra"

                        # Bonus fallback — check recursively if not found
                        if not sra_file.exists():
                            print(f"⚠️ Expected .sra not found at {sra_file}, checking subdirectories...")
                            for path in sra_temp_dir.glob(f"**/{run}.sra"):
                                sra_file = path
                                print(f"✅ Found alternative SRA file at: {sra_file}")
                                break

                        # Step 3: Convert .sra → .fastq
                        if sra_file.exists():
                            subprocess.run(["fasterq-dump", str(sra_file), "-O", str(run_dir)], check=True)
                            # Step 4: Clean up (delete .sra + its folder)
                            try:
                                sra_file.unlink()
                                if sra_subdir.exists():
                                    sra_subdir.rmdir()
                            except Exception as cleanup_err:
                                print(f"[WARN] Could not remove temp files: {cleanup_err}")

                            writer.writerow([run, "downloaded_NCBI", "prefetch/fasterq-dump",
                                            datetime.now().strftime("%Y-%m-%d %H:%M")])
                        else:
                            print(f"❌ Could not locate .sra file for {run} even after search.")
                            writer.writerow([run, "failed", "sra_not_found",
                                            datetime.now().strftime("%Y-%m-%d %H:%M")])

                    except subprocess.CalledProcessError as e:
                        print(f"❌ NCBI prefetch/fasterq-dump failed for {run}: {e}")
                        writer.writerow([run, "failed", str(e),
                                        datetime.now().strftime("%Y-%m-%d %H:%M")])

        print("[INFO] ✅ All downloads completed.")


#######################################################
# CLI entry point
#######################################################

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="SRA metadata fetcher & sequence downloader")
    parser.add_argument("--email", type=str, help="Your email (required by NCBI)")
    parser.add_argument("--organism", type=str, help="Organism name")
    parser.add_argument("--retmax", type=int, default=20, help="Number of samples to fetch")
    parser.add_argument("--download", action="store_true", help="Download sequences after fetching metadata")
    parser.add_argument("--runs", type=str, help="Comma-separated SRR IDs to download")
    args = parser.parse_args()

    ex = SRAExtractor(email=args.email)

    organism = args.organism or input("Enter organism name: ")
    df = ex.fetch_runinfo(organism, retmax=args.retmax)
    fname = f"SraRunInfo_{organism.replace(' ', '_')}.csv"
    ex.save_metadata(df, fname)

    # ✅ dynamic preview length
    try:
        n = int(input("Enter number of sequences to preview (default 5): ") or 5)
    except ValueError:
        n = 5
    ex.preview_metadata(fname, n=n)

    if args.download or input("Download sequences? (y/n): ").lower() in ("y", "yes"):
        if args.runs:
            run_ids = [r.strip() for r in args.runs.split(",")]
        else:
            run_ids = df["Run"].dropna().tolist() if "Run" in df.columns else []
        ex.download_runs(run_ids)
