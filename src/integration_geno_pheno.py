import pandas as pd
from pathlib import Path
import sys

class Integration:
    def __init__(self, seq_dir, metadata_file, ast_dir, assembly_paths, output):
        self.seq_dir = Path(seq_dir)
        self.metadata_file = Path(metadata_file)
        self.ast_dir = Path(ast_dir)

        # sample ID is the parent directory name of assembly.fasta ;map sample based on directory structure
        self.assembly_map = {Path(p).parent.name: Path(p) for p in assembly_paths}
        self.output = Path(output)

    def run(self):
        meta = pd.read_csv(self.metadata_file)
        meta = self.add_ast(meta)
        meta = self.add_sequences(meta)
        meta.to_csv(self.output, index=False)

    def add_ast(self, df):
        """This function will transform the corresponding antibiogram table to a pandas dataframe.
            And then merge it to the metadata dataframe based on biosample id"""
        # collect flattened AST rows
        rows = []
        for csv_path in self.ast_dir.glob("*.csv"):
            sample = csv_path.stem# sample ID from filename stem
            ast = pd.read_csv(csv_path).set_index("Antibiotic")    
            # flatten table: antibiotic_measurement to value/method/etc.
            flat = {}
            for antibiotic, row in ast.iterrows():
                for col, val in row.items():
                    flat[f"{antibiotic}_{col}"] = val

            rows.append({"sample": sample, **flat})
        ast_all = pd.DataFrame(rows)
        merged = df.merge(ast_all, left_on="BioSample", right_on="sample", how="left")
        return merged

    def add_sequences(self, df):
        """This function add its raw assembled draft genome to the metadata table under oclumn Assembled_seq"""
        # loop over assembly fasta paths mapped by sample
        for sample, fasta_path in self.assembly_map.items():
            if sample not in df["Run"].values:
                continue
            with open(fasta_path) as f:
                seq = self._strip_headers(f.read())#strip headers for each contig
            df.loc[df["Run"] == sample, "Assembled_seq"] = seq#seq contains raw DNA sequence

        return df

    def _strip_headers(self, data):
        """this funciton removes all FASTA header lines"""
        return "\n".join([l for l in data.splitlines() if not l.startswith(">")])


seq_dir = sys.argv[1]
metadata_file = sys.argv[2]
ast_dir = sys.argv[3]
assembly_paths = sys.argv[4].split()  # Snakemake passes space-separated list
output = sys.argv[5]

print("ARGS:", sys.argv)
print("seq_dir:", seq_dir)
print("metadata_file:", metadata_file)
print("ast_dir:", ast_dir)
print("assembly_paths:", assembly_paths)
print("output:", output)


Integration(seq_dir, metadata_file, ast_dir, assembly_paths, output).run()
