import pandas as pd
from pathlib import Path
import sys

class Integration:
    def __init__(self, seq_dir, metadata_file, ast_dir, assembly_paths,amr_dir, output):
        self.seq_dir = Path(seq_dir)
        self.metadata_file = Path(metadata_file)
        self.ast_dir = Path(ast_dir)
        self.amr_dir= amr_dir
        # sample ID is the parent directory name of assembly.fasta ;map sample based on directory structure
        self.assembly_map = {Path(p).parent.name: Path(p) for p in assembly_paths}
        self.output = Path(output)

    def run(self):
        self.integrated_df = pd.read_csv(self.metadata_file).copy()
        self.integrate_ast()
        self.integrate_sequences()
        self.integrate_amr()
        self.integrated_df.to_csv(self.output, index=False)

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

            rows.append({"sample": sample, **flat})
        ast_all = pd.DataFrame(rows)
        merged = self.integrated_df.merge(ast_all, left_on="BioSample", right_on="sample", how="left")
        return merged

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
        """Transform AMR data into a single row structure suitable for merging by biosample_id."""
        amr_df = self._transform_amr()
        self.integrated_df = self.integrated_df.merge(amr_df, left_on="BioSample", right_on="biosample_id", how="left")
    
    def _transform_amr(self):
        """Collect all AMR files and convert each table into a single row by reindexing columns to 'column_rowId'."""
        dict_amr_dfs = self._retrieve_amr()
        rows = []
        for df_key in dict_amr_dfs:
            # Source - https://stackoverflow.com/a; Posted by Scott Boston
            # Retrieved 2025-12-09, License - CC BY-SA 3.0
            t = dict_amr_dfs[df_key].copy()
            t.index = t.index + 1 #start count from *_1 instead of 0
            s = t.stack() #stack columns into series
            s.index = s.index.map('{0[1]}_{0[0]}'.format)#stack new columns by "columname_row_id"
            s["biosample_id"] = df_key #add biosample id for merging
            rows.append(s)

        return pd.DataFrame(rows)

    def _retrieve_amr(self):
        """Reads all amr files and converts it to a dict with biosample id as key"""
        all_amr_dfs = {}
        folder_path = Path(self.amr_dir)
        for amr_file in folder_path.iterdir():
            if amr_file.is_file():
                if amr_file.suffix.lower() in [".tsv"]:
                    all_amr_dfs[amr_file.stem] = pd.read_csv(amr_file, sep="\t")#amr_file.stem contains biosample id
        return all_amr_dfs

    def _strip_headers(self, data):
        """Removing all FASTA header lines"""
        return "\n".join([l for l in data.splitlines() if not l.startswith(">")])


seq_dir = sys.argv[1]
metadata_file = sys.argv[2]
ast_dir = sys.argv[3]
assembly_paths = sys.argv[4].split()  # Snakemake passes space-separated list
amr_dir = sys.argv[5]
output = sys.argv[6]

print("ARGS:", sys.argv)
print("seq_dir:", seq_dir)
print("metadata_file:", metadata_file)
print("ast_dir:", ast_dir)
print("assembly_paths:", assembly_paths)
print("amr_path:", amr_dir)
print("output:", output)


Integration(seq_dir, metadata_file, ast_dir, assembly_paths,amr_dir, output).run()
