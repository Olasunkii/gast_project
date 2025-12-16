import argparse
from pathlib import Path
import pandas as pd

class AMRTransformer:
    def __init__(self, amr_files, output_file):
        self.amr_files = [Path(f) for f in amr_files]
        self.output_file = Path(output_file)

    def run(self):
        amr_long_df = self._load_long_amr()
        if amr_long_df.empty:
            self.output_file.write_text("")
            return

        presence_df = self._build_gene_presence_matrix(amr_long_df)
        presence_df.to_csv(self.output_file, index=True)

    def _load_long_amr(self):
        dfs = []

        for amr_file in self.amr_files:
            run_id = amr_file.stem
            df = pd.read_csv(amr_file, sep="\t")

            if "Element symbol" not in df.columns:
                continue

            tmp = df[["Element symbol"]].copy()
            tmp["run_id"] = run_id
            dfs.append(tmp)

        if not dfs:
            return pd.DataFrame(columns=["run_id", "Element symbol"])

        amr_long = pd.concat(dfs, ignore_index=True)
        amr_long = amr_long.dropna(subset=["Element symbol"])
        amr_long["Element symbol"] = amr_long["Element symbol"].astype(str).str.strip()
        amr_long = amr_long[amr_long["Element symbol"] != ""]

        return amr_long

    def _build_gene_presence_matrix(self, amr_long):
        ct = pd.crosstab(amr_long["run_id"], amr_long["Element symbol"])
        presence = (ct > 0).astype(int)
        presence.columns = [f"gene_{col}_present" for col in presence.columns]
        presence.index.name = "run_id"
        return presence


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--amr_files", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    AMRTransformer(args.amr_files, args.output).run()
