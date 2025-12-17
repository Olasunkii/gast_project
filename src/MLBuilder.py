import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, RobustScaler

class MLBuilder:
    """
      - create one target column of carbapenem resistance
      - dataset split
      - scaling on conintinuous columns (fit on train only and not binary or ordinal columns)
    """
    def __init__(self, input_file, config_file, output_dir):
        self.input_file = Path(input_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.scaler_ = None
        self.target = None
        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)

    def run(self):
        self.df = pd.read_csv(self.input_file)
        self._derive_carbapenem_target() # derive target + remove leakage

        X = self.df.drop(columns=[self.target])
        y = self.df[self.target]

        (
            X_train,
            X_validation,
            X_test,
            y_train,
            y_validation,
            y_test,
        ) = self._split_data(X, y)

        X_train, X_validation, X_test = self._scale(
            X_train, X_validation, X_test
        )

        self._save_splits(
            X_train, y_train,
            X_validation, y_validation,
            X_test, y_test,
        )

    def _derive_carbapenem_target(self):
        """Computes a target column based of any resistance to a carbapenems.
            Then removes the corresponding columns"""
        carb_list = [c.lower() for c in self.config.get("carabapenems", [])]
        if not carb_list:
            raise ValueError("No carbapenems defined in config.")

        phenotype_cols = [
            col for col in self.df.columns
            if col.lower().endswith("_resistance phenotype")
            and any(col.lower().startswith(c + "_") for c in carb_list)
        ]

        if not phenotype_cols:
            raise ValueError("No carbapenem resistance phenotype columns found.")

        self.target = self.config["ml"]["column_target_name"]

        self.df[self.target] = (
            self.df[phenotype_cols].max(axis=1) >= 1
        ).astype(int)

        drop_cols = [
            col for col in self.df.columns
            if any(col.lower().startswith(c + "_") for c in carb_list)
        ]
        self.df.drop(columns=drop_cols, inplace=True)

    def _split_data(self, X, y):
        cfg = self.config["ml"]["split"]

        X_tmp, X_test, y_tmp, y_test = train_test_split(
            X,
            y,
            test_size=cfg["test_size"],
            random_state=cfg["random_state"],
            stratify=y if cfg.get("stratify", True) else None,
        )

        val_rel = cfg["validation_size"] / (1 - cfg["test_size"])

        X_train, X_validation, y_train, y_validation = train_test_split(
            X_tmp,
            y_tmp,
            test_size=val_rel,
            random_state=cfg["random_state"],
            stratify=y_tmp if cfg.get("stratify", True) else None,
        )

        return X_train, X_validation, X_test, y_train, y_validation, y_test

    def _scale(self, X_train, X_validation, X_test):
        """
        Scales only continuous numeric columns based on config.
        Binary / low-cardinality numeric columns are untouched.
        """
        scaling_cfg = self.config["ml"].get("scaling", {})
        method = scaling_cfg.get("method", "none").lower()

        if method == "none":
            self.scaler_ = None
            return X_train, X_validation, X_test

        num_cols = X_train.select_dtypes(include=[np.number]).columns
        #if more then 10 unique values then it is determined as column that should be scaled; excluding binary&ordinal columns
        num_cols = [c for c in num_cols if X_train[c].nunique() >= 10]

        if not num_cols:
            self.scaler_ = None
            return X_train, X_validation, X_test

        if method == "standard":
            scaler = StandardScaler()
        elif method == "robust":
            scaler = RobustScaler()
        else:
            raise ValueError(f"Unknown scaling method: {method}")

        self.scaler_ = scaler.fit(X_train[num_cols])

        for df in (X_train, X_validation, X_test):
            df.loc[:, num_cols] = scaler.transform(df[num_cols])

        return X_train, X_validation, X_test


    def _save_splits(
        self,
        X_train, y_train,
        X_validation, y_validation,
        X_test, y_test,
    ):
        X_train.to_csv(self.output_dir / "X_train.csv", index=False)
        y_train.to_csv(self.output_dir / "y_train.csv", index=False)

        X_validation.to_csv(self.output_dir / "X_validation.csv", index=False)
        y_validation.to_csv(self.output_dir / "y_validation.csv", index=False)

        X_test.to_csv(self.output_dir / "X_test.csv", index=False)
        y_test.to_csv(self.output_dir / "y_test.csv", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    checker = MLBuilder(args.input, args.config, args.output)
    checker.run()
