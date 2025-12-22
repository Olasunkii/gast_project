"""
sra_extractor.py
----------------
Lightweight, reusable SRA metadata and FASTQ downloader.

Usage example (in Jupyter or any script):
-----------------------------------------
from scripts.sra_extractor import SRAExtractor

ex = SRAExtractor(project_root="../", email="your@email.com")
df = ex.fetch_runinfo("Klebsiella pneumoniae", retmax=20)
ex.save_metadata(df, "SraRunInfo_kp.csv")
ex.preview_metadata("SraRunInfo_kp.csv", n=10)
filt = ex.filter_paired_illumina(df)
ex.preview_df(filt, n=10)
ex.download_runs(["SRR35784663"])
"""

from pathlib import Path
from typing import Optional, List
import pandas as pd
from Bio import Entrez
from io import StringIO
import time
import subprocess
import requests
from tqdm import tqdm
import csv


class SRAExtractor:
    def __init__(
        self,
        project_root: Optional[str] = ".",
        email: Optional[str] = None,
        default_retmax: int = 100,
        sleep_between_requests: float = 0.34,
    ):
        """
        project_root: Root folder of the project.
                      Always pass "../" when using inside notebooks.
        email: Email for NCBI Entrez (required).
        default_retmax: Default number of SRA records to fetch.
        sleep_between_requests: Delay between Entrez requests.
        """
        self.project_root = Path(project_root).resolve()
        print(f"[INFO] Using project root: {self.project_root}")

        # Define main folders
        self.data_meta = self.project_root / "data" / "metadata"
        self.data_seq = self.project_root / "data" / "sequences"

        # Ensure folders exist
        for p in (self.data_meta, self.data_seq):
            p.mkdir(parents=True, exist_ok=True)

        # Set Entrez email
        if email:
            Entrez.email = email
        if not getattr(Entrez, "email", None):
            raise ValueError("Entrez.email must be set. Pass email='you@example.com'")

        self.default_retmax = default_retmax
        self.sleep_between_requests = sleep_between_requests

    # ------------------------------------------------------------
    # 1️⃣  FETCH METADATA
    # ------------------------------------------------------------
    def search_sra(self, organism: str, retmax: Optional[int] = None) -> List[str]:
        """Return list of SRA record IDs for a given organism."""
        retmax = retmax or self.default_retmax
        term = f"{organism}[Organism]"
        handle = Entrez.esearch(db="sra", term=term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        return id_list

    def fetch_runinfo(
        self, organism_or_idlist, retmax: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Fetch run info for an organism or list of SRA IDs.
        Returns DataFrame of run info.
        """
        retmax = retmax or self.default_retmax

        if isinstance(organism_or_idlist, (list, tuple)):
            id_list = list(organism_or_idlist)
        else:
            id_list = self.search_sra(organism_or_idlist, retmax=retmax)
            if not id_list:
                raise ValueError(f"No SRA records found for '{organism_or_idlist}'")

        id_str = ",".join(id_list)
        attempts = 3
        for attempt in range(attempts):
            try:
                handle = Entrez.efetch(
                    db="sra", id=id_str, rettype="runinfo", retmode="text"
                )
                text = handle.read()
                handle.close()
                break
            except Exception as e:
                if attempt < attempts - 1:
                    print(f"[WARN] Retry {attempt+1}/{attempts} due to {e}")
                    time.sleep(self.sleep_between_requests * 3)
                else:
                    raise

        df = pd.read_csv(StringIO(text))
        self.last_metadata = df.copy()
        print(f"[INFO] Retrieved {len(df)} records.")
        return df

    def save_metadata(self, df: pd.DataFrame, filename: Optional[str] = None) -> Path:
        """Save metadata CSV to data/metadata."""
        if filename is None:
            filename = f"SraRunInfo_{int(time.time())}.csv"
        out_path = self.data_meta / filename
        df.to_csv(out_path, index=False)
        print(f"[INFO] Metadata saved → {out_path}")
        return out_path

    def preview_metadata(self, csvpath: Optional[str] = None, n: int = 10):
        """Preview metadata from file or last loaded DataFrame."""
        if csvpath:
            p = (
                self.data_meta / csvpath
                if not Path(csvpath).is_absolute()
                else Path(csvpath)
            )
            df = pd.read_csv(p, dtype=str)
        else:
            df = getattr(self, "last_metadata", None)
            if df is None:
                raise ValueError("No metadata in memory. Provide csvpath.")
        cols = [
            c
            for c in [
                "Run",
                "ScientificName",
                "LibraryLayout",
                "Platform",
                "size_MB",
                "bases",
            ]
            if c in df.columns
        ]
        print(df[cols or df.columns[:6]].head(n).to_string(index=False))

    def filter_paired_illumina(self, df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """Filter for paired-end Illumina runs."""
        df = df or getattr(self, "last_metadata", None)
        if df is None:
            raise ValueError("No metadata provided or loaded.")
        mask = df["LibraryLayout"].str.lower().str.contains("paired", na=False)
        mask &= df["Platform"].str.lower().str.contains("illumina", na=False)
        filt = df[mask].reset_index(drop=True)
        print(f"[INFO] Filtered to {len(filt)} paired-end Illumina runs.")
        self.last_filtered = filt
        return filt

    def preview_df(self, df: pd.DataFrame, n: int = 10):
        """Quick preview of DataFrame."""
        cols = [
            c
            for c in [
                "Run",
                "ScientificName",
                "LibraryLayout",
                "Platform",
                "size_MB",
                "bases",
            ]
            if c in df.columns
        ]
        print(df[cols or df.columns[:6]].head(n).to_string(index=False))

    # ------------------------------------------------------------
    # 2️⃣  DOWNLOAD FASTQ FILES
    # ------------------------------------------------------------
    def _get_ena_fastq_links(self, run_id: str):
        """Get FASTQ URLs for a run from ENA."""
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run_id}&result=read_run&fields=fastq_ftp&format=tsv"
        r = requests.get(url, timeout=30)
        if r.status_code != 200 or "fastq_ftp" not in r.text:
            return []
        lines = r.text.strip().split("\n")
        if len(lines) < 2:
            return []
        fastq_links = lines[1].split("\t")[-1].split(";")
        return [
            f"https://{link.strip()}" if not link.startswith("http") else link
            for link in fastq_links
        ]

    def _download_file(self, url, outpath):
        """Download a file with progress bar."""
        outpath = Path(outpath)
        tmp = outpath.with_suffix(outpath.suffix + ".part")
        headers = {}
        if tmp.exists():
            existing = tmp.stat().st_size
            headers["Range"] = f"bytes={existing}-"
        else:
            existing = 0
        r = requests.get(url, stream=True, headers=headers, timeout=60)
        total = int(r.headers.get("content-length", 0)) + existing
        mode = "ab" if existing > 0 else "wb"
        with open(tmp, mode) as f, tqdm(
            total=total, initial=existing, unit="B", unit_scale=True, desc=outpath.name
        ) as pbar:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
        tmp.rename(outpath)
        return outpath

    def download_runs(self, run_ids, prefer_ena=True):
        """Download FASTQ files for SRR IDs."""
        if isinstance(run_ids, str):
            run_ids = [run_ids]
        log_file = self.data_seq / "download_log.csv"
        log_exists = log_file.exists()
        with open(log_file, "a", newline="") as log:
            writer = csv.writer(log)
            if not log_exists:
                writer.writerow(["Run", "Status", "Files", "Timestamp"])
            for run in run_ids:
                run_dir = self.data_seq / run
                run_dir.mkdir(exist_ok=True)
                if any(run_dir.glob("*.fastq*")):
                    print(f"✅ {run} already downloaded — skipping.")
                    writer.writerow(
                        [run, "exists", "existing", time.strftime("%Y-%m-%d %H:%M")]
                    )
                    continue
                print(f"\n🔽 Downloading {run}...")
                urls = []
                if prefer_ena:
                    try:
                        urls = self._get_ena_fastq_links(run)
                    except Exception as e:
                        print(f"⚠️ ENA lookup failed for {run}: {e}")
                if urls:
                    print(f"Found {len(urls)} FASTQ file(s) on ENA.")
                    files = []
                    for url in urls:
                        out = run_dir / Path(url).name
                        self._download_file(url, out)
                        files.append(str(out))
                    writer.writerow(
                        [
                            run,
                            "downloaded",
                            ";".join(files),
                            time.strftime("%Y-%m-%d %H:%M"),
                        ]
                    )
                else:
                    print(f"⚠️ No ENA FASTQ found for {run}. Trying NCBI prefetch...")
                    try:
                        subprocess.run(["prefetch", run], check=True)
                        subprocess.run(
                            ["fasterq-dump", run, "-O", str(run_dir)], check=True
                        )
                        writer.writerow(
                            [
                                run,
                                "downloaded_ncbi",
                                "prefetch/fasterq-dump",
                                time.strftime("%Y-%m-%d %H:%M"),
                            ]
                        )
                    except Exception as e:
                        print(f"❌ NCBI prefetch failed for {run}: {e}")
                        writer.writerow(
                            [run, "failed", str(e), time.strftime("%Y-%m-%d %H:%M")]
                        )
