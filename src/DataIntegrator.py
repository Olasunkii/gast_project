import pandas as pd
from pathlib import Path
import sys

class Integration:
    def __init__(self, seq_dir, sample_metadata_file, host_metadata_file, ast_dir, assembly_paths,amr_dir, output):
        self.seq_dir = Path(seq_dir)
        self.sample_metadata_file = Path(sample_metadata_file)
        self.host_metadata_file = Path(host_metadata_file)
        self.ast_dir = Path(ast_dir)
        self.amr_dir= amr_dir
        # sample ID is the parent directory name of assembly.fasta ;map sample based on directory structure
        self.assembly_map = {Path(p).parent.name: Path(p) for p in assembly_paths}
        self.output = Path(output)

    def run(self):
        base = pd.read_csv(self.sample_metadata_file)        
        self.integrated_df = base[["BioSample", "Run"]].copy() # Load only the columns required for merging
        print(self.integrated_df)
        self.integrate_host_metadata()
        self.integrate_ast()
        self.integrate_amr()
        self.integrate_sequences()
        self.integrated_df.to_csv(self.output, index=False)

    def integrate_host_metadata(self):
        """Reads the hostmetata and merges it based on biosample id."""
        host_metadata_df = pd.read_csv(self.host_metadata_file)
        self.integrated_df = self.integrated_df.merge(
            host_metadata_df, on="BioSample", how="left"      
        )

    def integrate_ast(self):
        """Transforms the corresponding antibiogram table to a pandas dataframe.
            Then merges it to the metadata dataframe based on biosample id"""
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

            rows.append({"sample_id": sample, **flat})
        ast_all = pd.DataFrame(rows)
        self.integrated_df = self.integrated_df.merge(ast_all, left_on="BioSample", right_on="sample_id", how="left")

    def integrate_sequences(self):
        """Adding the raw assembled draft genome to the metadata table under clumn Assembled_seq"""
        # loop over assembly fasta paths mapped by sample
        for sample, fasta_path in self.assembly_map.items():
            if sample not in self.integrated_df["Run"].values:
                continue
            with open(fasta_path) as f:
                seq = self._strip_headers(f.read())#strip headers for each contig
            self.integrated_df.loc[self.integrated_df["Run"] == sample, "Assembled_seq"] = seq#seq contains raw DNA sequence

    def integrate_amr(self):
        amr_df = self._transform_amr()
        self.integrated_df = self.integrated_df.merge(amr_df, left_on="Run", right_on="run_id", how="left")

    def _transform_amr(self):
        """Transform AMR data into a single row structure suitable for merging by biosample_id."""
        folder_path = Path(self.amr_dir)
        rows = []
        for amr_file in folder_path.glob("*.tsv"):
            id_ = amr_file.stem  # retrieve run/sample identifier from file stem
            df = pd.read_csv(amr_file, sep="\t") 
            df.index = df.index + 1  # start row numbering at 1
            record = {} 

            for idx, amr_row in df.iterrows():
                for col, val in amr_row.items():  # iterate over each column in row
                    key = f"{col}_{idx}"  # construct unique column name
                    record[key] = val  # assign AMR value under flattened name

            record["run_id"] = id_  # add identifier
            rows.append(record)  # collect record into result set

        return pd.DataFrame(rows)

    def _strip_headers(self, data):
        """Removing all FASTA header lines"""
        return "\n".join([l for l in data.splitlines() if not l.startswith(">")])


seq_dir = sys.argv[1]
sample_metadata_file = sys.argv[2]
host_metadata_file = sys.argv[3]
ast_dir = sys.argv[4]
assembly_paths = sys.argv[5].split()  # Snakemake passes space-separated list
amr_dir = sys.argv[6]
output = sys.argv[7]

# print("ARGS:", sys.argv)

Integration(seq_dir, sample_metadata_file, host_metadata_file, ast_dir, assembly_paths,amr_dir, output).run()
