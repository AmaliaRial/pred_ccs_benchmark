import pandas as pd
import os
from pathlib import Path

#cuidado con la doble S en el nombre la carpeta

# .../repos/pred_css_benchmark/pred_ccs_benchmark
BASE_DIR = Path(__file__).resolve().parent

# .../repos/pred_css_benchmark
ROOT_DIR = BASE_DIR.parent

# .../repos/pred_css_benchmark/datasets
DATA_DIR = ROOT_DIR / "datasets"

# .../repos/pred_css_benchmark/pred_ccs_benchmark/data
OUT_DIR = BASE_DIR / "data"
OUT_DIR.mkdir(exist_ok=True, parents=True)

#TODO ver por que no mete los nombres de dataset
def clean_ccsbase(input_file, output_file):
    """
    Limpia ccsbase_descriptors:
      id,name,adduct,m/z,ccs,smi,type,z,ref,ccs_type,ccs_method,is_3D
    """
    df = pd.read_csv(input_file, low_memory=False)

    df_clean = pd.DataFrame()
    df_clean["dataset"] = "ccsbase" #TODO ver por que no mete los nombres de dataset
    df_clean["id"] = df["id"]
    df_clean["name"] = df["name"]

    # estructura: solo SMILES aquí
    df_clean["smiles"] = df["smi"]
    df_clean["inchi"] = None  # no hay InChI en este

    df_clean["adduct"] = df["adduct"]
    df_clean["ccs"] = df["ccs"]
    df_clean["is_3D"] = df.get("is_3D", None)
    df_clean["gas"] = "N2"

    # los dejo estos por si aca
    df_clean["mz"] = df.get("m/z", None)
    df_clean["type"] = df.get("type", None)
    df_clean["charge"] = df.get("z", None)

    df_clean.to_csv(output_file, index=False)
    print(f"[INFO] Saved CCSbase cleaned dataset -> {output_file}")


def clean_allccs(input_file, output_file):
    """
    Limpia AllCCS2_experimental_with_inchis_descriptors.csv

    Columnas:
    AllCCS ID,Name,Structure,Formula,Type,Adduct,m/z,CCS,
    Confidence level,Update date,InChI,is_3D
    """
    df = pd.read_csv(input_file, low_memory=False)

    df_clean = pd.DataFrame()
    df_clean["dataset"] = "allccs"
    df_clean["id"] = df["AllCCS ID"]
    df_clean["name"] = df["Name"]

    # Structure la usamos como SMILES
    df_clean["smiles"] = df["Structure"]
    df_clean["inchi"] = df["InChI"]

    df_clean["adduct"] = df["Adduct"]
    df_clean["ccs"] = df["CCS"]
    df_clean["is_3D"] = df.get("is_3D", None)
    df_clean["gas"] = "N2"

    # info opcional
    df_clean["formula"] = df.get("Formula", None)
    df_clean["mz"] = df.get("m/z", None)
    df_clean["type"] = df.get("Type", None)

    df_clean.to_csv(output_file, index=False)
    print(f"[INFO] Saved AllCCS cleaned dataset -> {output_file}")


def clean_metlin_lipids(input_file, output_file):
    """
    Limpia METLIN-CCS-Lipids_descriptors:
      Name,Formula,RT,MainPositiveAdduct,...,InChI,...,
      CCS [M+Na]+, CCS [M+H]+, ..., CCS [M-H+FA]-, ... ,is_3D

    Convierte de formato ancho (varios CCS por fila) a formato largo
    (una fila por (molécula, adduct) con CCS != NaN).
    """
    df = pd.read_csv(input_file, low_memory=False)

    # columnas de CCS para distintos aductos
    ccs_cols = [c for c in df.columns if c.startswith("CCS [")]

    id_vars = ["Name", "Formula", "InChI", "is_3D"]
    # nos aseguramos de que existen
    for col in id_vars:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in METLIN lipids file.")

    melted = df.melt(
        id_vars=id_vars,
        value_vars=ccs_cols,
        var_name="ccs_adduct_col",
        value_name="ccs"
    )

    # quitamos filas sin CCS
    melted = melted.dropna(subset=["ccs"])

    # extraemos el adduct del nombre de columna, ej: "CCS [M+Na]+"
    def extract_adduct(colname: str) -> str:
        # asume formato "CCS [M+Na]+"
        # -> devuelve "[M+Na]+"
        try:
            return colname.split("CCS ")[1]
        except IndexError:
            return None

    melted["adduct"] = melted["ccs_adduct_col"].apply(extract_adduct)

    df_clean = pd.DataFrame()
    df_clean["dataset"] = "metlinlipidims"
    df_clean["id"] = None  # este fichero no tiene ID numérico claro
    df_clean["name"] = melted["Name"]
    df_clean["smiles"] = None  # no hay smiles
    df_clean["inchi"] = melted["InChI"]
    df_clean["adduct"] = melted["adduct"]
    df_clean["ccs"] = melted["ccs"]
    df_clean["is_3D"] = melted["is_3D"]
    df_clean["gas"] = "N2"

    # info opcional
    df_clean["formula"] = melted["Formula"]

    df_clean.to_csv(output_file, index=False)
    print(f"[INFO] Saved METLIN lipids cleaned dataset -> {output_file}")


def clean_metlin_ims(input_file, output_file, sep="\t"):
    """
    Limpia METLIN_IMS_descriptors.tsv:

    Molecule Name Molecular Formula METLIN ID Precursor Adduct CCS1 CCS2 CCS3
    CCS_AVG % CV m/z Adduct m/z.1 Dimer Dimer.1 dimer line CCS m/z.2
    pubChem inchi smiles InChIKEY is_3D
    """ #TODO ver bien la separacion de los nombres en el archivo original
    df = pd.read_csv(input_file, sep=sep, low_memory=False)

    df_clean = pd.DataFrame()
    df_clean["dataset"] = "metlinims"
    df_clean["id"] = df["METLIN ID"]
    df_clean["name"] = df["Molecule Name"]

    # estructura: tenemos smiles e inchi
    df_clean["smiles"] = df["smiles"]
    df_clean["inchi"] = df["inchi"]

    df_clean["adduct"] = df["Adduct"]
    df_clean["ccs"] = df["CCS"]
    df_clean["is_3D"] = df.get("is_3D", None)
    df_clean["gas"] = "N2"

    # opcional: m/z
    df_clean["mz"] = df.get("m/z", None)

    df_clean.to_csv(output_file, index=False)
    print(f"[INFO] Saved METLIN IMS cleaned dataset -> {output_file}")


def main():

    clean_ccsbase(
        DATA_DIR / "ccsbase_descriptors.csv",
        OUT_DIR / "dataset_ccsbase.csv",
    )

    clean_allccs(
        DATA_DIR / "AllCCS2_experimental_with_inchis_descriptors.csv",
        OUT_DIR / "dataset_allccs.csv",
    )

    clean_metlin_lipids(
        DATA_DIR / "METLIN-CCS-Lipids_descriptors.csv",
        OUT_DIR / "dataset_metlinlipidims.csv",
    )

    clean_metlin_ims(
        DATA_DIR / "METLIN_IMS_descriptors.tsv",
        OUT_DIR / "dataset_metlinims.csv",
        sep=",",
    )


if __name__ == "__main__":
    main()
