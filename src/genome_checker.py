import pandas as pd
from pathlib import Path
import yaml
import argparse


class GenomeQC:
    def __init__(self, input_folder, config_file, output_file):
        self.input_folder = Path(input_folder)
        self.output_file = Path(output_file)

        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)["genome_completeness_qc"]

    def run(self):
        """Assess the checkm output for each sample folder in input folder dir."""
        failures = []
        # for each sample folder in input folder direcotry
        for sample_dir in self.input_folder.iterdir():
            if not sample_dir.is_dir():
                continue

            report = sample_dir / "quality_report.tsv"
            if not report.exists():
                continue

            df = pd.read_csv(report, sep="\t")
            failures.extend(self._check(sample_dir.name, df))

        if failures:
            pd.DataFrame(failures).to_csv(self.output_file, sep="\t", index=False)
        else:
            with open(self.output_file, "w") as f:
                f.write(
                    "All genomes pass completeness, contamination, and contig thresholds.\n"
                )

    def _check(self, sample_name, df):
        """if isolate fails on genome checks it will be documented."""
        acc = []
        comp_min = self.config["completeness"]["min"]
        cont_max = self.config["contamination"]["max"]
        contigs_max = self.config["total_contigs"]["max"]

        for _, row in df.iterrows():
            complete = row.get("Completeness")
            cont = row.get("Contamination")
            contigs = row.get("Total_Contigs")

            fail = []
            if pd.notna(complete) and complete < comp_min:
                fail.append(f"Completeness<{comp_min}")

            if pd.notna(cont) and cont > cont_max:
                fail.append(f"Contamination>{cont_max}")

            if pd.notna(contigs) and contigs > contigs_max:
                fail.append(f"Contigs>{contigs_max}")

            if fail:
                acc.append(
                    {
                        "sample": sample_name,
                        "Completeness": complete,
                        "Contamination": cont,
                        "Total_Contigs": contigs,
                        "Failures": ";".join(fail),
                    }
                )

        return acc


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    qc = GenomeQC(args.input, args.config, args.output)
    qc.run()
