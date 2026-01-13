import argparse
import pandas as pd
from pathlib import Path
import yaml

class Preprocessor:
    """
    - encoding categorical columns; only manual based on config file
    - dropping columns
    - splitting columns (location, mic values and longitude and latitude)
    """

    def __init__(self, input_file, config_file, output_file):
        self.input_file = Path(input_file)
        self.output_file = Path(output_file)
        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)

    def run(self):
        """Preprocesses the data to splitting, dropping and encoding columns."""
        self.df = pd.read_csv(self.input_file)
        self._run_structural()
        self._drop_columns_initial()  # drops columns that do not have to be encoded
        self._encode()
        self._drop_nan_columns()  # drops columns only containing nan
        self.df.to_csv(self.output_file, index=False)

    def _run_structural(self):
        """Splitting of location, combination of antibiotics and latitude and longitude columns.
        Drops the original columns that were replaced by the newly split columns."""
        if self.config["preprocessing"]["structural"]["split_location"]:
            self._split_location()

        if self.config["preprocessing"]["structural"]["split_mic"]:
            self._split_mic()
        if self.config["preprocessing"]["structural"]["split_latlon"]:
            self._split_latlon()

    def _encode(self):
        """Encode selected categorical values to standardized numeric codes based on NCBI conventions,
        then one-hot encode all remaining categorical columns."""
        enc = self.config.get("preprocessing", {}).get("encoding", {})
        if enc.get("encode_comparison_signs"):
            self._apply_encoding(
                column_substrings=enc["sign_column_contains"],
                mapping=enc["comparison_signs"],
                normalize_lower=False,
            )
        if enc.get("encode_phenotypes"):
            ph = enc["phenotype"]
            self._apply_encoding(
                column_substrings=ph["phenotype_column_contains"],
                mapping=ph["mapping"],
                normalize_lower=True,
            )
        if enc.get("encode_units"):
            un = enc["units"]
            self._apply_encoding(
                column_substrings=un["units_column_contains"],
                mapping=un["mapping"],
                normalize_lower=True,
            )
        self._encode_bool_cols()

    def _apply_encoding(self, *, column_substrings, mapping, normalize_lower):
        """Applies encoding by corresponding mapping defined in the config.yaml"""
        substrings = [s.lower() for s in column_substrings]
        cols = [c for c in self.df.columns if any(s in c.lower() for s in substrings)]

        if normalize_lower:
            mapping = {k.lower(): v for k, v in mapping.items()}
        else:
            mapping = {k.strip(): v for k, v in mapping.items()}

        for col in cols:
            s = self.df[col].astype(str).str.strip()
            if normalize_lower:
                s = s.str.lower()
            self.df[col] = s.replace({"nan": None}).map(mapping)

    def _encode_bool_cols(self):
        bool_cols = self.df.select_dtypes(include=["bool"]).columns.tolist()
        self.df[bool_cols] = self.df[bool_cols].astype(int)

    def _split_latlon(self):
        """Splits latitude and longitude to each a column including corresponding sign convention for direction"""
        col = self.config["preprocessing"]["location"]["latlon_column"]

        lat_list = []
        lon_list = []

        for raw in self.df[col].astype(str):
            parts = raw.split()
            lat_v = float(parts[0])
            lat_d = parts[1]
            lon_v = float(parts[2])
            lon_d = parts[3]
            # adding sign convention for direction
            if lat_d == "S":
                lat_v = -lat_v
            if lon_d == "W":
                lon_v = -lon_v

            lat_list.append(lat_v)
            lon_list.append(lon_v)

        self.df["latitude"] = lat_list
        self.df["longitude"] = lon_list
        self.df.drop(columns=[col], inplace=True)

    def _split_mic(self):
        """Split MIC values into _primary and _secondary columns; for single antibiotic, _secondary is set to NaN."""
        targets = self.config["preprocessing"]["mic_split"]["column_contains"]
        exclude = self.config["preprocessing"]["mic_split"]["exclude_contains"]

        meas_cols = [
            c
            for c in self.df.columns
            if any(t in c.lower() for t in targets)
            and not any(e in c.lower() for e in exclude)
        ]

        for col in meas_cols:
            p_list = []
            i_list = []

            for v in self.df[col].astype(str):
                v = v.strip()
                if "/" in v:
                    p, i = (x.strip() for x in v.split("/", 1))
                    p_list.append(float(p) if p not in ["", "nan"] else None)
                    i_list.append(float(i) if i not in ["", "nan"] else None)
                else:
                    p_list.append(float(v) if v not in ["", "nan"] else None)
                    i_list.append(None)

            self.df[f"{col}_primary"] = p_list
            self.df[f"{col}_secondary"] = i_list

        self.df.drop(columns=meas_cols, inplace=True)

    def _split_location(self):
        """Splits location into two columns: country and province/state"""
        col = self.config["preprocessing"]["location"]["source_column"]
        dest = self.config["preprocessing"]["location"]["target_columns"]
        self.df[dest] = self.df[col].str.split(":", expand=True)
        self.df.drop(columns=[col], inplace=True)

    def _drop_columns_initial(self):
        explicit = self.config["preprocessing"]["drop"]["explicit"]
        existing = [c for c in explicit if c in self.df.columns]
        self.df.drop(columns=existing, inplace=True)

        contains = self.config["preprocessing"]["drop"]["contains"]
        for pattern in contains:
            cols = [c for c in self.df.columns if pattern in c.lower()]
            self.df.drop(columns=cols, inplace=True)
        print("Dropped explicit:", len(existing))
        print("Dropped pattern-matched:", contains)

    def _drop_nan_columns(self):
        """Drops columns only containing nan values"""
        nan_only = [c for c in self.df.columns if self.df[c].isna().all()]
        self.df.drop(columns=nan_only, inplace=True)

        print("Dropped NaN-only:", len(nan_only))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    preprocessor = Preprocessor(args.input, args.config, args.output)
    transformed = preprocessor.run()
