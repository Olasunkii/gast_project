"""
sra_extractor.py
Lightweight SRA metadata fetcher & helper functions for a Jupyter workflow.

Usage:
    from scripts.sra_extractor import SRAExtractor
    ex = SRAExtractor(project_root=".", email="youremail@example.com")
    df = ex.fetch_runinfo("Klebsiella pneumoniae", retmax=50)
    ex.save_metadata(df, "data/metadata/SraRunInfo_kp.csv")
    ex.preview_metadata("data/metadata/SraRunInfo_kp.csv", n=10)
    filt = ex.filter_paired_illumina(df)
    ex.preview_df(filt, n=10)
"""

from pathlib import Path
from typing import Optional, List
import pandas as pd
from Bio import Entrez
from io import StringIO
import time
###########
import subprocess
import requests
from tqdm import tqdm
import shutil
import csv

class SRAExtractor:
    def __init__(self,
                 project_root: Optional[str] = ".",
                 email: Optional[str] = None,
                 default_retmax: int = 100,
                 sleep_between_requests: float = 0.34):
        """
        project_root: root folder for the project (module will create data/metadata, data/sequences)
        email: email for NCBI Entrez (required by NCBI policy). Set Entrez.email if provided.
        default_retmax: default number of run records to request
        sleep_between_requests: delay between Entrez requests (avoid rate limits)
        """
        self.project_root = Path(project_root).resolve()
        self.data_meta = self.project_root / "data" / "metadata"
        self.data_seq = self.project_root / "data" / "sequences"
        self.sleep_between_requests = sleep_between_requests
        self.default_retmax = default_retmax

        # Create directories
        for p in (self.data_meta, self.data_seq):
            p.mkdir(parents=True, exist_ok=True)

        # Entrez email
        if email:
            Entrez.email = email
        if not getattr(Entrez, "email", None):
            raise ValueError("Entrez.email must be set (pass email argument when creating SRAExtractor).")

    def search_sra(self, organism: str, retmax: Optional[int] = None) -> List[str]:
        """
        Use Entrez.esearch to get SRA IDs (SRA record IDs)
        Returns list of SRA ids (these are SRA internal ids, not SRR run accessions).
        """
        retmax = retmax or self.default_retmax
        term = f'{organism}[Organism]'
        handle = Entrez.esearch(db="sra", term=term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        return id_list

    def fetch_runinfo(self, organism_or_idlist, retmax: Optional[int] = None) -> pd.DataFrame:
        """
        Two modes:
        - if organism_or_idlist is a string and looks like organism name: do esearch -> efetch(runinfo)
        - if it's a list of ids (SRA ids) pass them to efetch directly
        Returns a DataFrame of runinfo (columns as in SRA runinfo).
        """
        retmax = retmax or self.default_retmax

        # if a list is given, treat as id list; if string with comma-separated numeric ids,
        # also allow that; otherwise treat as organism name
        if isinstance(organism_or_idlist, (list, tuple)):
            id_list = list(organism_or_idlist)
        elif isinstance(organism_or_idlist, str) and organism_or_idlist.strip().isdigit():
            id_list = [organism_or_idlist.strip()]
        else:
            # treat as organism name
            id_list = self.search_sra(organism_or_idlist, retmax=retmax)
            if not id_list:
                raise ValueError(f"No SRA records found for '{organism_or_idlist}'")

        # efetch runinfo
        id_str = ",".join(id_list)
        # use retries in case NCBI hiccups
        attempts = 3
        for attempt in range(attempts):
            try:
                handle = Entrez.efetch(db="sra", id=id_str, rettype="runinfo", retmode="text")
                text = handle.read()
                handle.close()
                break
            except Exception as e:
                if attempt < attempts - 1:
                    time.sleep(self.sleep_between_requests * 3)
                    continue
                else:
                    raise

        # ensure string
        if isinstance(text, bytes):
            text = text.decode("utf-8")

        df = pd.read_csv(StringIO(text))
        # keep a copy in the object if desired
        self.last_metadata = df.copy()
        return df

    def save_metadata(self, df: pd.DataFrame, filename: Optional[str] = None) -> Path:
        """
        Save metadata DataFrame to CSV in data/metadata.
        If filename is None, uses a sensible default with timestamp.
        Returns Path to saved CSV.
        """
        if filename is None:
            # simple timestamp filename
            fname = f"SraRunInfo_{int(time.time())}.csv"
            out = self.data_meta / fname
        else:
            out = Path(filename)
            if not out.is_absolute():
                out = self.data_meta / out

        df.to_csv(out, index=False)
        return out

    def preview_metadata(self, csvpath: Optional[str] = None, n: int = 10):
        """
        Print a preview of metadata CSV or of in-memory last_metadata.
        If csvpath provided, load it from data/metadata or absolute path.
        """
        if csvpath:
            p = Path(csvpath)
            if not p.is_absolute():
                p = self.data_meta / p
            if not p.exists():
                raise FileNotFoundError(f"{p} not found")
            df = pd.read_csv(p, dtype=str)
        else:
            if getattr(self, "last_metadata", None) is None:
                raise ValueError("No metadata loaded in memory and no csvpath provided.")
            df = self.last_metadata

        # print a small, relevant subset for the user
        display_cols = []
        for c in ["Run", "ScientificName", "LibraryLayout", "Platform", "size_MB", "bases"]:
            if c in df.columns:
                display_cols.append(c)

        if not display_cols:
            # fallback to first 6 columns
            display_cols = df.columns.tolist()[:6]

        print(df[display_cols].head(n).to_string(index=False))

    def filter_paired_illumina(self, df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        Filter for paired-end Illumina runs with non-empty bases/size columns.
        Returns filtered DataFrame (copy).
        """
        if df is None:
            df = getattr(self, "last_metadata", None)
            if df is None:
                raise ValueError("No metadata provided or loaded.")
        # make a copy to avoid modifying original
        dfc = df.copy()

        # normalize column names if they differ
        col_map = {c.lower(): c for c in dfc.columns}
        # safe checks
        libcol = None
        platcol = None
        for k, v in col_map.items():
            if k in ("librarylayout", "library_layout"):
                libcol = v
            if k in ("platform",):
                platcol = v

        if libcol is None:
            raise ValueError("Library layout column not found in metadata (expected 'LibraryLayout').")

        # filter: paired & platform contains ILLUMINA & bases or size_MB not null
        mask = dfc[libcol].str.lower().str.contains("paired", na=False)
        if platcol:
            mask = mask & dfc[platcol].str.lower().str.contains("illumina", na=False)
        # prefer 'bases' column if present else size_MB
        if "bases" in dfc.columns:
            mask = mask & dfc["bases"].notna()
        elif "size_MB" in dfc.columns:
            mask = mask & dfc["size_MB"].notna()

        filtered = dfc[mask].reset_index(drop=True)
        # store filtered
        self.last_filtered = filtered
        return filtered

    def preview_df(self, df: pd.DataFrame, n: int = 10):
        """Show a friendly preview of a provided DataFrame (subset of useful columns)."""
        display_cols = []
        for c in ["Run", "ScientificName", "LibraryLayout", "Platform", "size_MB", "bases"]:
            if c in df.columns:
                display_cols.append(c)
        if not display_cols:
            display_cols = df.columns.tolist()[:6]
        print(df[display_cols].head(n).to_string(index=False))

    #######################################################

    def _get_ena_fastq_links(self, run_id: str):
        """
        Query ENA for FASTQ download URLs for a given run_id.
        Returns list of URLs (usually 2 for paired-end, 1 for single-end).
        """
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run_id}&result=read_run&fields=fastq_ftp&format=tsv"
        r = requests.get(url, timeout=30)
        if r.status_code != 200 or "fastq_ftp" not in r.text:
            return []
        lines = r.text.strip().split("\n")
        if len(lines) < 2:
            return []
        # second line has FTP links (comma-separated)
        fastq_links = lines[1].split("\t")[-1].split(";")
        # Convert ftp:// to https://
        fastq_links = [f"https://{link.strip()}" if not link.startswith("http") else link for link in fastq_links]
        return fastq_links

    def _download_file(self, url, outpath):
        """
        Download a single file with progress bar using requests + tqdm.
        Supports resume (if server supports Range).
        """
        outpath = Path(outpath)
        tmp_path = outpath.with_suffix(outpath.suffix + ".part")

        headers = {}
        if tmp_path.exists():
            existing_size = tmp_path.stat().st_size
            headers["Range"] = f"bytes={existing_size}-"
        else:
            existing_size = 0

        response = requests.get(url, stream=True, headers=headers, timeout=60)
        total_size = int(response.headers.get("content-length", 0)) + existing_size
        mode = "ab" if existing_size > 0 else "wb"

        with open(tmp_path, mode) as f, tqdm(
            total=total_size, initial=existing_size, unit="B", unit_scale=True, desc=outpath.name
        ) as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))

        tmp_path.rename(outpath)
        return outpath

    def download_runs(self, run_ids, prefer_ena=True):
        """
        Download FASTQ files for given SRR run IDs.
        Checks ENA first, then falls back to NCBI prefetch if needed.
        Logs downloads in data/download_log.csv.
        """
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
                downloaded_files = []

                # Skip if already downloaded
                if any(run_dir.glob("*.fastq*")):
                    print(f"✅ {run} already downloaded — skipping.")
                    writer.writerow([run, "exists", "existing", time.strftime("%Y-%m-%d %H:%M")])
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
                    for url in urls:
                        fname = Path(url).name
                        outpath = run_dir / fname
                        if not outpath.exists():
                            self._download_file(url, outpath)
                        downloaded_files.append(str(outpath))
                    writer.writerow([run, "downloaded", ";".join(downloaded_files), time.strftime("%Y-%m-%d %H:%M")])
                else:
                    print(f"⚠️ No ENA FASTQ found for {run}. Attempting NCBI prefetch...")
                    try:
                        subprocess.run(["prefetch", run], check=True)
                        subprocess.run(["fasterq-dump", run, "-O", str(run_dir)], check=True)
                        writer.writerow([run, "downloaded_ncbi", "prefetch/fasterq-dump", time.strftime("%Y-%m-%d %H:%M")])
                    except Exception as e:
                        print(f"❌ NCBI prefetch failed for {run}: {e}")
                        writer.writerow([run, "failed", str(e), time.strftime("%Y-%m-%d %H:%M")])
