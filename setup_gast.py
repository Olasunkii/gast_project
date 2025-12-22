#!/usr/bin/env python3
"""
setup_gast.py
---------------------------------
Initial setup script for GAST project structure and environment.

Run this once after cloning:
    python setup_gast.py
"""

from pathlib import Path
import yaml

# --- 1️⃣ Define project structure ---
def init_structure(base_dir="project_gast"):
    base = Path(base_dir).resolve()
    folders = [
        base / "data" / "metadata",
        base / "data" / "sequences",
        base / "data" / "ast",
        base / "results"
    ]

    for folder in folders:
        folder.mkdir(parents=True, exist_ok=True)
        print(f"[OK] Created: {folder}")

    # --- 2️⃣ Create config.yaml if not exists ---
    config_path = base / "configs" / "config.yaml" #config without using settings
    if not config_path.exists():
        config = {
            "paths": {
                "metadata_dir": "data/metadata",
                "host_metadata_dir": "data/host_metadata",
                "sequences_dir": "data/sequences",
                "ast_dir": "data/ast",
                "results_dir": "results",
                "reference_db_dir": "ref_db",
            },
            "carbapenems": [
                "imipenem",
                "meropenem",
                "doripenem",
                "ertapenem",
                "biapenem",
            ],
            "breakpoints": { #current EUCAST carabapenem breakingpoints 2025
                "doripenem": [1, 2],
                "ertapenem": [0.5, 0.5],
                "imipenem": [2, 4],
                "imipenem-relebactam": [2, 2],
                "meropenem": [2, 2],
                "meropenem-vaborbactam": [8, 8],
            },
            "preprocessing": {
                "structural": {
                    "split_location": True,
                    "split_mic": True,
                    "split_latlon": True,
                },
                "location": {
                    "source_column": "geographic location",
                    "target_columns": ["Country", "State/Province"],
                    "latlon_column": "latitude and longitude",
                },
                "mic_split": {
                    "column_contains": ["measurement"],
                    "exclude_contains": ["sign", "unit"],
                },
                "encoding": {
                    "encode_units": True,
                    "units": {
                        "units_column_contains": ["units"],
                        "mapping": {
                            "mg/l": 1,
                            "mm": 2,
                        },
                    },
                    "encode_comparison_signs": True,
                    "sign_column_contains": [
                        "measurement sign",
                        "_sign",
                    ],
                    "comparison_signs": {
                        "<": -2,
                        "<=": -1,
                        "==": 0,
                        "=": 0,
                        ">=": 1,
                        ">": 2,
                    },
                    "encode_phenotypes": True,
                    "phenotype": {
                        "phenotype_column_contains": [
                            "phenotype",
                            "resistance phenotype",
                            "_phen",
                        ],
                        "eucast_phenotype_mapping": {
                            "S": 0,
                            "I": 1,
                            "R": 2,
                        },
                        "mapping": {
                            "n": -1,
                            "nd": -1,
                            "not defined": -1,
                            "s": 0,
                            "sensitive": 0,
                            "susceptible": 0,
                            "i": 1,
                            "intermediate": 1,
                            "ns": 1,
                            "nonsusceptible": 1,
                            "ssd": 1,
                            "susceptible-dose dependent": 1,
                            "r": 2,
                            "resistant": 2,
                            "hlar": 2,
                            "high level aminoglycoside resistance": 2,
                        },
                    },
                },
                "drop": {
                    "explicit": [
                        "BioSample",
                        "Run",
                        "isolate",
                        "collected by",
                        "collection date",
                        "host subject id",
                        "sample",
                        "organism_complex",
                    ],
                    "contains": [
                        "vendor",
                        "laboratory typing",
                        "testing standard",
                        "_id",
                    ],
                },
            },
        }

        with open(config_path, "w") as f:
            yaml.dump(config, f)
        print(f"[OK] Config file created: {config_path}")
    else:
        print(f"[SKIP] Config file already exists: {config_path}")

    # --- 3️⃣ .gitignore setup ---
    gitignore_path = base / ".gitignore"
    if not gitignore_path.exists():
        gitignore = """# Ignore data and cache files
data/*
!data/metadata/*.csv
output/*
__pycache__/
*.pyc
.ipynb_checkpoints/
.env
.DS_Store
"""
        gitignore_path.write_text(gitignore)
        print("[OK] .gitignore created.")
    else:
        print("[SKIP] .gitignore already exists.")

    print("\n✅ GAST project setup complete.\n")


if __name__ == "__main__":
    init_structure(".")
