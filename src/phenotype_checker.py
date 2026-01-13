import argparse
import pandas as pd
from pathlib import Path
import yaml

class PhenotypeChecker:
    """Checks whether the resistance label and dosage are in line with EUCAST defintions stated in the conf.yaml.
        Input: Folder path to AST CSV/ Antibiogram"""
    def __init__(self, ast_input_folder, config_file, output_file):
        self.ast_input_folder = Path(ast_input_folder)
        self.output_file = Path(output_file)

        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)

        self.bp = self.config["breakpoints"]

    def run(self):
        """Execute consistency checks on input CSV files and write the results.
            Input: path to AST folder
            Output: CSV file containing isolates that fall outside the EUCAST range."""
        out = []
        for file in self.ast_input_folder.glob("*.csv"):
            df = pd.read_csv(file)
            if not self._has_required(df):
                continue
            out.extend(self._check_file(file.stem, df))

        if out:
            pd.DataFrame(out).to_csv(self.output_file, sep="\t", index=False)
        else:
            with open(self.output_file, "w") as f:
                f.write(
                    "All carbapenem resistance phenotypes consistent with EUCAST clinical breakpoints.\n"
                )

    def _has_required(self, df):
        """Check whether the input DataFrame contains the required columns."""
        c = df.columns
        return "Antibiotic" in c and "Measurement" in c and "Resistance phenotype" in c

    def _check_file(self, sample_name, df):
        """Validate reported resistance phenotypes against expected classifications."""
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
                acc.append(
                    {
                        "sample": sample_name,
                        "antibiotic": ab,
                        "mic_value": mic,
                        "reported_sir": reported,
                        "expected_sir": expected,
                    }
                )

        return acc

    def _classify(self, ab, mic):
        """Classifies antibiotic susceptibility based on Minimum Inhibitory Concentration (MIC)."""
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
        """Normalizes various string synonyms of SIR phenotypes to standard codes."""
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
