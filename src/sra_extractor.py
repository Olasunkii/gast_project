#!/usr/bin/env python3
"""
sra_extractor.py
Interactive SRA metadata fetcher & sequence downloader.
Works in both Jupyter Notebook and CLI.
"""

import csv
from Bio import Entrez
from Bio.Entrez import Parser
from pathlib import Path
import pandas as pd
import io
import time
from typing import Optional, List
from tqdm import tqdm
import subprocess
import shutil
import requests
import xml.etree.ElementTree as ET

# Helper: run command and raise with readable error
def run_cmd(cmd, **kwargs):
    try:
        subprocess.run(cmd, check=True, **kwargs)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{e}")


class SRAExtractor:
    def __init__(
        self,
        project_root: Optional[str] = ".",
        email=None,
        default_retmax: int = 20,
        sleep_between_requests: float = 0.34,
    ):
        """
        project_root: root folder for the project (creates data/metadata, data/sequences)
        email: NCBI Entrez email (required)
        default_retmax: default number of run records to request
        sleep_between_requests: delay between Entrez requests to avoid rate limits
        """
        self.project_root = Path(project_root).resolve()
        self.data_meta = self.project_root / "data" / "metadata"
        self.data_seq = self.project_root / "data" / "sequences"
        self.ast_dir = self.project_root / "data" / "ast"
        self.host_metadata_dir = self.project_root / "data" / "host_metadata"
        self.data_meta.mkdir(parents=True, exist_ok=True)
        self.data_seq.mkdir(parents=True, exist_ok=True)
        self.ast_dir.mkdir(parents=True, exist_ok=True)
        self.host_metadata_dir.mkdir(parents=True, exist_ok=True)
        self.CARBAPENEMS = {
            "imipenem",
            "meropenem",
            "ertapenem",
            "doripenem",
            "biapenem",
        }

        self.default_retmax = default_retmax
        self.sleep_between_requests = sleep_between_requests

        # Email setup
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
        print(f"Antimicrobial susceptibility testing (AST) folder: {self.ast_dir}")
        print(f"Host metadata folder: {self.host_metadata_dir}")

    #######################################################
    # Metadata fetching
    #######################################################

    def search_sra(self, organism: str, retmax: Optional[int] = None) -> List[str]:
        retmax = retmax or self.default_retmax
        term = f"{organism}[Organism]"
        handle = Entrez.esearch(db="sra", term=term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        return record.get("IdList", [])

    def fetch_runinfo(
        self, organism_or_idlist: str, retmax: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Fetch SRA run metadata for a given organism name or list of IDs.
        Saves CSV metadata file to data/metadata.
        """
        import io

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
            with Entrez.efetch(
                db="sra", id=",".join(ids), rettype="runinfo", retmode="text"
            ) as handle:
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
        if filename is None:
            fname = f"SraRunInfo_{int(time.time())}.csv"
        else:
            fname = filename
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

        display_cols = [
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
        if not display_cols:
            display_cols = df.columns.tolist()[:6]

        print(df[display_cols].head(n).to_string(index=False))

    #######################################################
    # Sequence download
    #######################################################
    def _get_ena_fastq_links(self, run_id: str) -> list:
        """
        Fetch FASTQ download URLs for a given SRR run from ENA.
        """
        url = (
            "https://www.ebi.ac.uk/ena/portal/api/filereport"
            f"?accession={run_id}&result=read_run&fields=fastq_ftp&format=tsv"
        )

        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            lines = r.text.strip().split("\n")

            if len(lines) < 2:
                return []

            # last column is fastq_ftp; may contain semicolon-separated ftp paths
            fastq_field = lines[1].split("\t")[-1].strip()
            if not fastq_field:
                return []

            links = []
            for part in fastq_field.split(";"):
                part = part.strip()
                if not part:
                    continue
                # ENA returns ftp paths like ftp.sra.ebi.ac.uk/vol1/fastq/...
                if (
                    part.startswith("ftp://")
                    or part.startswith("http://")
                    or part.startswith("https://")
                ):
                    links.append(part)
                else:
                    links.append("https://" + part)  # convert ftp path to https
            return links

        except Exception as e:
            print(f"[WARN] Could not retrieve ENA links for {run_id}: {e}")
            return []

    def _download_file(self, url, outpath):
        """
        Robust streaming download with resume support.
        Returns Path or None on failure.
        """
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
            total=total_size,
            initial=existing_size,
            unit="B",
            unit_scale=True,
            desc=outpath.name,
            leave=True,
        ) as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))

        tmp_path.rename(outpath)
        return outpath

    def _compress_fastq(self, fq_path: Path, threads: int = 4):
        """
        Compress a .fastq file to .fastq.gz using pigz if available, else gzip.
        Removes original .fastq on success.
        """
        pigz = shutil.which("pigz")
        if pigz:
            run_cmd([pigz, "-p", str(threads), str(fq_path)])
        else:
            run_cmd(["gzip", "-f", str(fq_path)])
        gz_path = fq_path.with_suffix(fq_path.suffix + ".gz")
        return gz_path

    def download_runs(self, run_ids, prefer_ena=True, threads: int = 4):
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
                run_dir.mkdir(parents=True, exist_ok=True)
                downloaded_files = []

                # quick skip if gz already present
                if any(run_dir.glob("*.fastq.gz")):
                    print(f"✅ {run} already downloaded — skipping.")
                    writer.writerow(
                        [run, "exists", "existing", time.strftime("%Y-%m-%d %H:%M")]
                    )
                    continue

                print(f"\n🔽 Downloading {run}...")
                urls = []
                if prefer_ena:
                    urls = self._get_ena_fastq_links(run)

                if urls:
                    print(f"Found {len(urls)} FASTQ file(s) on ENA.")
                    for url in urls:
                        fname = Path(url).name
                        outpath = run_dir / fname
                        if not outpath.exists():
                            res = self._download_file(url, outpath)
                            if res is None:
                                print(f"[WARN] Failed to download {url}")
                                continue

                        # If file is not gz, compress it
                        if outpath.suffix != ".gz" and not outpath.name.endswith(
                            ".fastq.gz"
                        ):
                            # assume it's fastq; compress in place
                            self._compress_fastq(outpath, threads=threads)
                            outpath = outpath.with_suffix(outpath.suffix + ".gz")
                        downloaded_files.append(str(outpath))
                    writer.writerow(
                        [
                            run,
                            "downloaded",
                            ";".join(downloaded_files),
                            time.strftime("%Y-%m-%d %H:%M"),
                        ]
                    )

                else:
                    print(
                        f"⚠️ No ENA FASTQ found for {run}. Attempting NCBI prefetch..."
                    )
                    try:
                        # Force prefetch to put the .sra file inside the run_dir
                        prefetch_cmd = [
                            "prefetch",
                            run,
                            "--output-directory",
                            str(run_dir),
                        ]
                        run_cmd(prefetch_cmd)

                        # find the sra file under run_dir (prefetch may create nested dirs)
                        sra_candidates = list(run_dir.rglob(f"{run}.sra"))
                        if not sra_candidates:
                            raise FileNotFoundError(
                                f".sra file for {run} not found under {run_dir}"
                            )
                        sra_path = sra_candidates[0]

                        # Convert SRA -> FASTQ
                        fasterq_cmd = [
                            "fasterq-dump",
                            str(sra_path),
                            "-O",
                            str(run_dir),
                            "--threads",
                            str(threads),
                        ]
                        run_cmd(fasterq_cmd)

                        # After successful conversion, delete .sra immediately (Option C)
                        if sra_path.exists():
                            sra_path.unlink()
                            print(f"🗑️ Deleted SRA file: {sra_path}")

                        # Compress any produced .fastq files
                        gz_files = []
                        for fq in run_dir.glob("*.fastq"):
                            gz = self._compress_fastq(fq, threads=threads)
                            gz_files.append(str(gz))

                        if not gz_files:
                            raise RuntimeError("fasterq-dump produced no .fastq files")

                        writer.writerow(
                            [
                                run,
                                "downloaded_ncbi",
                                ";".join(gz_files),
                                time.strftime("%Y-%m-%d %H:%M"),
                            ]
                        )
                        downloaded_files = gz_files

                    except Exception as e:
                        print(f"❌ NCBI fallback failed for {run}: {e}")
                        writer.writerow(
                            [run, "failed", str(e), time.strftime("%Y-%m-%d %H:%M")]
                        )

        print("[INFO] All downloads completed.")

    def find_sra_ids_with_antibiogram(
        self, organism: str, retmax: int, start: int = 0
    ) -> list:
        """Searches for BioSamples that are from the corresponding organism
        and contains a antibiogram table. And links the found biosample ids to a SRA run identifiers
        """
        term = f"{organism}[Organism] AND antibiogram[filter]"
        h = Entrez.esearch(db="biosample", term=term, retstart=start, retmax=retmax)
        r = Entrez.read(h)
        biosample_ids = r.get("IdList", [])
        sra_ids = []
        # link biosample to sra run id
        for bid in biosample_ids:
            try:
                link = Entrez.elink(dbfrom="biosample", db="sra", id=bid)
                record = Entrez.read(link)
            except (IOError, RuntimeError, Parser.ValidationError):
                continue
    
            for item in record[0].get("LinkSetDb", []):
                if item.get("DbTo") == "sra":
                    for link_item in item.get("Link", []):
                        sra_ids.append(link_item["Id"])
        return sra_ids

    def fetch_sra_metadata(self, sra_ids: list) -> pd.DataFrame:
        """retrieves sra metadata based on sra_ids list. Returns a dataframe only containing paired sequences."""
        df = pd.DataFrame()
        if not sra_ids:
            return df
        with Entrez.efetch(
            db="sra", id=",".join(sra_ids), rettype="runinfo", retmode="text"
        ) as handle:
            raw = handle.read()
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        df = pd.read_csv(io.StringIO(str(raw)))
        df = df[df["LibraryLayout"] == "PAIRED"]  # only paired sequences
        return df

    def fetch_antibiogram_and_host(
        self, biosample_id: str
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Fetches BioSample XML once, extracts antibiogram table and host metadata, returns two dataframes."""
        if not biosample_id:
            return pd.DataFrame(), pd.DataFrame()

        try:
            h = Entrez.efetch(db="biosample", id=biosample_id, rettype="xml")
            data = h.read()
            root = ET.fromstring(data)
        except:
            return pd.DataFrame(), pd.DataFrame()

        # extract antibiogram
        table = root.find(".//Table")
        if table is None:
            antibiogram_df = pd.DataFrame()
        else:
            header = table.find("Header").findall("Cell")
            columns = [c.text for c in header]
            rows = []
            for row in table.find("Body").findall("Row"):
                cells = [cell.text or "" for cell in row.findall("Cell")]
                rows.append(cells)
            antibiogram_df = pd.DataFrame(rows, columns=columns)

        if not self._contains_carbapenem(antibiogram_df):  # if no carabapenem present
            antibiogram_df = pd.DataFrame()

        # extract host metadata
        attrs = root.findall(".//Attributes/Attribute")
        host = {}
        for attr in attrs:
            name = (
                attr.get("display_name")
                or attr.get("attribute_name")
                or attr.get("harmonized_name")
            )
            if not name:
                continue
            host[name] = attr.text or ""

        host_df = pd.DataFrame([host])

        return antibiogram_df, host_df

    def _contains_carbapenem(self, df: pd.DataFrame) -> bool:
        if df.empty:
            return False
        # normalize column names
        cols = {c.lower(): c for c in df.columns}
        drug_col = None
        for key in ("antibiotic", "drug", "agent"):
            if key in cols:
                drug_col = cols[key]
                break

        if drug_col is None:
            return False

        drugs = df[drug_col].astype(str).str.lower().str.strip()
        return drugs.isin(self.CARBAPENEMS).any()

    def collect_resistant_metadata(
        self, organism: str, retmax: int, batchsize_: int = 20
    ) -> pd.DataFrame:
        """Keep searching SRA id, antibiogram filtering, host-metadata aggregation.
        Until retmax number has been reached and then saves sample and host metadata."""
        collected_metadata = []
        collected_host_metadata = []
        start = 0

        while len(collected_metadata) < retmax:
            sra_ids = self.find_sra_ids_with_antibiogram(organism, retmax, start)
            sample_metadata_df = self.fetch_sra_metadata(sra_ids)
            self._save_batch(
                sample_metadata_df, collected_metadata, collected_host_metadata
            )
            start += batchsize_  # move to next batch

        return self._write_outputs(collected_metadata, collected_host_metadata)

    def _save_batch(self, metadata_df, collected_metadata, collected_host_metadata):
        """Processes one sample metadata batch: saves only samples with antibiograms and corresponding host metadata."""
        for biosample_id in metadata_df["BioSample"]:

            ast, host = self.fetch_antibiogram_and_host(biosample_id)
            if ast.empty:
                continue

            ast.to_csv(self.ast_dir / f"{biosample_id}.csv", index=False)

            collected_host_metadata.append(host.assign(BioSample=biosample_id))

            row = metadata_df[metadata_df["BioSample"] == biosample_id]
            collected_metadata.append(row)

    def _write_outputs(self, collected_metadata, collected_host_metadata):
        """Combines all host metadata into one file and returns the combined sample metadata table."""
        if not collected_metadata:
            return pd.DataFrame()

        if collected_host_metadata:
            host_df = pd.concat(collected_host_metadata, ignore_index=True)
            host_df.to_csv(
                self.host_metadata_dir / "host_metadata_all.csv", index=False
            )

        return pd.concat(collected_metadata, ignore_index=True)


#######################################################
# Optional CLI entry point
#######################################################

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="SRA metadata fetcher & sequence downloader"
    )
    parser.add_argument("--email", type=str, help="Your email (required by NCBI)")
    parser.add_argument("--organism", type=str, help="Organism name")
    parser.add_argument(
        "--retmax", type=int, default=20, help="Number of samples to fetch"
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download sequences after fetching metadata",
    )
    parser.add_argument("--runs", type=str, help="Comma-separated SRR IDs to download")
    parser.add_argument(
        "--threads", type=int, default=4, help="Threads for fasterq-dump / pigz"
    )
    args = parser.parse_args()

    ex = SRAExtractor(email=args.email)
    organism = args.organism or input("Enter organism name: ")
    # any sequence:fetch_runinfo(); only paired:fetch_runinfo_paired();known AST sequences:resistant_paired_metadata
    df = ex.collect_resistant_metadata(organism, args.retmax)
    fname = f"SraRunInfo_{organism.replace(' ', '_')}.csv"
    ex.save_metadata(df, fname)  # save metadata to csv file: data/metadata
    ex.preview_metadata(fname, n=5)

    if args.download or input("Download sequences? (y/n): ").lower() in ("y", "yes"):
        if args.runs:
            run_ids = [r.strip() for r in args.runs.split(",")]
        else:
            run_ids = df["Run"].dropna().tolist()
        ex.download_runs(run_ids, threads=args.threads)
