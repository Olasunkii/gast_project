#!/usr/bin/env python3
"""
setup_gast.py
---------------------------------
Initial setup script for GAST project structure and environment.

Run this once after cloning:
    python setup_gast.py
"""

import os
from pathlib import Path
import yaml

# --- 1️⃣ Define project structure ---
def init_structure(base_dir="project_gast"):
    base = Path(base_dir).resolve()
    folders = [
        base / "data" / "metadata",
        base / "data" / "sequences",
        base / "output",
        base / "notebooks",
        base / "scripts"
    ]

    for folder in folders:
        folder.mkdir(parents=True, exist_ok=True)
        print(f"[OK] Created: {folder}")

    # --- 2️⃣ Create config.yaml if not exists ---
    config_path = base / "config.yaml"
    if not config_path.exists():
        config = {
            "project_name": "Genomic Analysis Support Tool (GAST)",
            "paths": {
                "data_dir": str(base / "data"),
                "metadata_dir": str(base / "data" / "metadata"),
                "sequence_dir": str(base / "data" / "sequences"),
                "output_dir": str(base / "output")
            },
            "email": "your_email@example.com",
            "ena_preferred": True,
            "default_retmax": 50
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
        print(f"[OK] .gitignore created.")
    else:
        print(f"[SKIP] .gitignore already exists.")

    print("\n✅ GAST project setup complete.\n")


if __name__ == "__main__":
    init_structure(".")
