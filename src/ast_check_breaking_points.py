import sys
import pandas as pd
from pathlib import Path

# MIC breakpoints (mg/L): (S_threshold, R_threshold)
EUCAST_BP = {
    "doripenem":    (1, 1),
    "ertapenem":    (0.5, 0.5),
    "imipenem":     (2, 2),
    "meropenem":    (2, 2),   # non-meningitis
}

def classify_eucast(ab, mic):
    if pd.isna(mic):
        return None
    mic = float(mic)
    s_thr, r_thr = EUCAST_BP[ab]
    if mic <= s_thr:
        return "S"
    if mic > r_thr:
        return "R"
    return "I"

def normalize_ab(name):
    n = name.strip().lower()
    for k in EUCAST_BP:
        if k in n:
            return k
    return None

def normalize_phenotypic_label(x):
    if not isinstance(x, str):
        return None
    x = x.strip().lower()
    if x in ("s", "susceptible", "sensitive"):
        return "S"
    if x in ("i", "intermediate"):
        return "I"
    if x in ("r", "resistant", "ns", "nonsusceptible"):
        return "R"
    return None

def main():
    input_folder = Path(sys.argv[1])
    out_file = Path(sys.argv[2])

    inconsistencies = []

    for file in input_folder.glob("*.csv"):
        df = pd.read_csv(file)

        # required columns
        if (
            "Antibiotic" not in df.columns
            or "Measurement" not in df.columns
            or "Resistance phenotype" not in df.columns
        ):
            continue

        for _, row in df.iterrows():
            measure = str(row["Measurement"])
            if "/" in measure: # skip combination antibiotics (e.g.,MIC "8.0/4.0")
                continue

            ab_norm = normalize_ab(str(row["Antibiotic"]))#only retrieve antibiotic name
            if ab_norm is None:
                continue
            reported_sir = normalize_phenotypic_label(row["Resistance phenotype"])
            mic = row["Measurement"]

            try:
                expected = classify_eucast(ab_norm, mic)
            except ValueError:
                continue   # skip malformed MIC values

            if expected is None or reported_sir is None:
                continue
            #found inconsistency with eucast breaking point
            if expected != reported_sir:
                inconsistencies.append({
                    "sample": file.stem,
                    "antibiotic": ab_norm,
                    "mic_value": mic,
                    "reported_sir": reported_sir,
                    "expected_sir": expected
                })

    if inconsistencies:
        pd.DataFrame(inconsistencies).to_csv(out_file, sep="\t", index=False)
    else:
        with open(out_file, "w") as f:
            f.write("All carbapenem SIR calls are consistent with EUCAST.\n")

if __name__ == "__main__":
    main()
