from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
TMP_DIR = BASE_DIR / "tmp"

TMP_DIR.mkdir(exist_ok=True)

def make_input_for_deepccs(dataset: str):
    src = DATA_DIR / f"dataset_{dataset}.csv"
    if not src.exists():
        raise FileNotFoundError(f"Dataset not found: {src}")

    df = pd.read_csv(src)

    out_df = pd.DataFrame()
    out_df["SMILES"] = df["smiles"]

    # DeepCCS espera aductos en formato tipo M+H, M+Na, M-H, etc.
    #toca limpiarlos
    def normalize_adduct(a: str) -> str:
        # ej: "[M+H]+" -> "M+H"
        if isinstance(a, str):
            a = a.strip()
            if a.startswith("[") and "]" in a:
                a = a[1:a.index("]")]
        return a

    out_df["Adducts"] = df["adduct"].apply(normalize_adduct)

    out_path = TMP_DIR / f"{dataset}_for_deepccs.csv"
    out_df.to_csv(out_path, index=False)
    print(f"[INFO] Wrote DeepCCS input -> {out_path}")

if __name__ == "__main__":
    make_input_for_deepccs("allccs")
