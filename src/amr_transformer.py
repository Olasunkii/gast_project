import argparse
from pathlib import Path
import pandas as pd

class AMRTransformer:
    def __init__(self, amr_dir, output_file):
        self.amr_dir = Path(amr_dir)
        self.output_file = Path(output_file)

    def run(self):
        amr_long_df = self._load_long_amr()
        if amr_long_df.empty:
            self.output_file.write_text("")#if no amr files return empty
            return

        presence_df = self._build_gene_presence_matrix(amr_long_df)
        presence_df.to_csv(self.output_file, index=True)

    def _load_long_amr(self):
        """Read all amr files and concatenate into one"""
        dfs = []

        for amr_file in self.amr_dir.glob("*.tsv"):
            run_id = amr_file.stem  # bv. "SRR12345"
            df = pd.read_csv(amr_file, sep="\t")
            if "Element symbol" not in df.columns:#if no gene symbol column present skip
                continue
            tmp = df[["Element symbol"]].copy()
            tmp["run_id"] = run_id
            dfs.append(tmp)

        if not dfs:#if no amr files return empty
            return pd.DataFrame(columns=["run_id", "Element symbol"])

        amr_long = pd.concat(dfs, ignore_index=True)

        # Drop nan values
        amr_long = amr_long.dropna(subset=["Element symbol"])
        amr_long["Element symbol"] = amr_long["Element symbol"].astype(str).str.strip()
        amr_long = amr_long[amr_long["Element symbol"] != ""]#check for empty strings

        return amr_long

    def _build_gene_presence_matrix(self, amr_long):
        """Transforms the annotated genes into present/absent matrix. if gene present 1 otherwise 0.
            run_id-column for merging"""
        # Crosstab: counts per (run_id, gene)
        ct = pd.crosstab(amr_long["run_id"], amr_long["Element symbol"])
        # everything > 0 → 1
        presence = (ct > 0).astype(int)
        presence.columns = [
            f"gene_{col}_present"
            for col in presence.columns
        ]
        presence.index.name = "run_id"
        return presence

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--amr_dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    builder = AMRTransformer(args.amr_dir, args.output)
    builder.run()
