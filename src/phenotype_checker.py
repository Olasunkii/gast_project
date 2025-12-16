import pandas as pd
from pathlib import Path
import yaml
import argparse

class PhenotypeChecker:
    def __init__(self, input_folder, config_file, output_file):
        self.input_folder = Path(input_folder)
        self.output_file = Path(output_file)

        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)

        self.bp = self.config["breakpoints"]

    def run(self):
        out = []
        for file in self.input_folder.glob("*.csv"):
            df = pd.read_csv(file)
            if not self._has_required(df):
                continue
            out.extend(self._check_file(file.stem, df))

        if out:
            pd.DataFrame(out).to_csv(self.output_file, sep="\t", index=False)
        else:
            with open(self.output_file, "w") as f:
                f.write("All carbapenem resistance phenotypes consistent with EUCAST clinical breakpoints.\n")

    def _has_required(self, df):
        c = df.columns
        return (
            "Antibiotic" in c
            and "Measurement" in c
            and "Resistance phenotype" in c
        )

    def _check_file(self, sample_name, df):
        acc = []

        for _, row in df.iterrows():
            m_raw = str(row["Measurement"])
            if "/" in m_raw:
                continue

            ab = self._normalize_ab(str(row["Antibiotic"]))
            if ab is None:
                continue

            mic = row["Measurement"]
            reported = self._normalize_sir(row["Resistance phenotype"])
            expected = self._classify(ab, mic)

            if expected is None or reported is None:
                continue
            if expected != reported:
                acc.append({
                    "sample": sample_name,
                    "antibiotic": ab,
                    "mic_value": mic,
                    "reported_sir": reported,
                    "expected_sir": expected
                })

        return acc

    def _classify(self, ab, mic):
        if pd.isna(mic):
            return None
        try:
            mic = float(mic)
        except ValueError:
            return None

        s_thr, r_thr = self.bp[ab]
        if mic <= s_thr:
            return "S"
        if mic > r_thr:
            return "R"
        return "I"

    def _normalize_ab(self, name):
        n = name.strip().lower()
        for k in self.bp:
            if k in n:
                return k
        return None

    def _normalize_sir(self, x):
        if not isinstance(x, str):
            return None
        v = x.strip().lower()
        if v in ("s", "susceptible", "sensitive"):
            return "S"
        if v in ("i", "intermediate"):
            return "I"
        if v in ("r", "resistant", "ns", "nonsusceptible"):
            return "R"
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    checker = PhenotypeChecker(args.input, args.config, args.output)
    checker.run()
