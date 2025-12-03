import pandas as pd
import sys
from pathlib import Path

class integration():
    def __init__(self, sequence_folder, metadata_path):
        self.sequence_folder_path = sequence_folder
        self.metadata_folder_path = metadata_path
    
    def integrate(self):
        df_metadata = pd.read_csv(self.metadata_folder_path)
        df_merged_ast = self.integrate_ast(df_metadata)
        df_integrated = self.integrate_sequences(df_merged_ast)
        print("Writing to:", output_file)
        df_integrated.to_csv(output_file, index=False)

    def integrate_ast(self, df):
        rows = []
        for csv_file_name in Path(ast_folder_path).glob("*.csv"):
            ast_df = pd.read_csv(csv_file_name)
            sample = csv_file_name.stem #get only the biosample id based on filename
            ast_df = ast_df.set_index("Antibiotic") # set index on antibiotic column

            # flatten by turning each col row into a key:value
            flat = {}
            for antibiotic, row in ast_df.iterrows():
                for col, val in row.items():
                    new_key = f"{antibiotic}_{col}"
                    flat[new_key] = val
            rows.append({"sample": sample, **flat})

        ast_df_all = pd.DataFrame(rows)#this merges all antibiotic columns; if no existing key then value is set to NAN
        merged = df.merge(ast_df_all,left_on="BioSample",right_on="sample",how="left")
        return merged

    def integrate_sequences(self, df_metadata):
        list_seq_Id= df_metadata["Run"].values
        for sequence_id in list_seq_Id:
            filepath = assembly_folder_path
            with open(filepath) as f:
                data = f.read()
                cleaned_data = self.remove_headers(data)
    
            df_metadata.loc[df_metadata["Run"] == sequence_id, "Assembled_seq"] = cleaned_data
        return df_metadata

    def remove_headers(self,sequence_data):
        clean_lines = []
        for line in sequence_data.splitlines():
            if line.startswith(">"):
                continue
            clean_lines.append(line)
        clean_data = "\n".join(clean_lines)
        return clean_data
#collect all parameters
seq_folder_path=sys.argv[1]
metadata_folder_path=sys.argv[2]
ast_folder_path=sys.argv[3]
assembly_folder_path=sys.argv[4]
output_file= sys.argv[5]
#start integration of phenotypic and genotype information into one
integration_model = integration(seq_folder_path, metadata_folder_path)
integration_model.integrate()