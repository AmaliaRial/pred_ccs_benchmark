#!/usr/bin/env python
import argparse
import os
import subprocess
from pathlib import Path

import pandas as pd

# ---- CONFIG ----

DEFAULT_DATASETS = ["allccs", "metlinims", "metlinlipidims", "ccsbase"]
DEFAULT_TOOLS = ["darkchem", "allccs", "ccsbase", "deepccs", "ccsp2", "hyperccs"]
# cuando el modelo esté listo: DEFAULT_TOOLS.append("ourmodel")

BASE_DIR = Path(__file__).resolve().parent          # .../pred_ccs_benchmark
DATA_DIR = BASE_DIR / "data"                       # cleaned datasets
PRED_DIR = BASE_DIR / "predictions"                # predictions/<dataset>/<tool>.csv
RESULTS_DIR = BASE_DIR / "results"                 # results/<dataset>/<tool>/
BENCHMARK_DIR = BASE_DIR / "benchmark"             # global tables etc.

CREATE_REPORT_SCRIPT = BASE_DIR / "create_report.py"



def run_create_report(dataset, tool, dataset_file, pred_file, out_dir):
    """
    Llama a create_report:
    python create_report.py <tool> <dataset> <file_original> <file_results>
    y deja los outputs (joined, metrics, report) en out_dir.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python",
        str(CREATE_REPORT_SCRIPT),
        tool,
        dataset,
        str(dataset_file),
        str(pred_file),
    ]
    print(f"[INFO] Running create_report for dataset={dataset}, tool={tool}")
    subprocess.run(cmd, check=True, cwd=str(out_dir))


def run_benchmark(datasets, tools):
    """
    Para cada dataset y herramienta:
      - busca el CSV de dataset limpio
      - busca el CSV de predicciones
      - llama a create_report
    """
    for dataset in datasets:
        # tus ficheros son data/dataset_allccs.csv, dataset_metlinims.csv, etc.
        dataset_file = DATA_DIR / f"dataset_{dataset}.csv"

        if not dataset_file.exists():
            print(f"[WARN] Dataset file not found: {dataset_file}, skipping dataset {dataset}")
            continue

        for tool in tools:
            # CSV de predicciones estándar: predictions/<dataset>/<tool>.csv
            pred_file = PRED_DIR / dataset / f"{tool}.csv"
            if not pred_file.exists():
                print(f"[WARN] Predictions file not found: {pred_file}, skipping tool {tool} for dataset {dataset}")
                continue

            out_dir = RESULTS_DIR / dataset / tool
            run_create_report(dataset, tool, dataset_file, pred_file, out_dir)


def join_metrics(datasets, tools, output_file=None):
    """
    Lee todos los metrics<tool>.csv de results/<dataset>/<tool>/
    y los junta en un solo CSV.
    """
    rows = []
    for dataset in datasets:
        for tool in tools:
            metrics_path = RESULTS_DIR / dataset / tool / f"metrics{tool}.csv"
            if not metrics_path.exists():
                continue
            df = pd.read_csv(metrics_path, sep=";")
            rows.append(df)

    if not rows:
        print("[WARN] No metrics files found.")
        return

    data = pd.concat(rows, ignore_index=True)

    BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)
    if output_file is None:
        output_file = BENCHMARK_DIR / "joined_metrics.csv"

    data.to_csv(output_file, index=False)
    print(f"[INFO] Saved joined metrics to {output_file}")


def build_compounds_table(dataset, tools, output_file=None):
    """
    Construye la “tabla 2”:
    - fila = compuesto (SMILES + Adduct)
    - columnas:
        * info básica (SMILES, Adduct, CCS real)
        * CCS_<tool>, err_perc_<tool>, mz_<tool>, etc.
    """
    base_df = None

    for tool in tools:
        joined_path = RESULTS_DIR / dataset / tool / f"joined{tool}.csv"
        if not joined_path.exists():
            print(f"[WARN] joined file not found: {joined_path}, skipping tool {tool}")
            continue

        df = pd.read_csv(joined_path, sep=";")
        df_tool = df.copy()

        id_cols = ["SMILES", "Adduct"]

        if "CCS" not in df_tool.columns:
            raise ValueError(f"Column 'CCS' not found in {joined_path}")

        col_map = {}
        if "Predicted CCS" in df_tool.columns:
            col_map["Predicted CCS"] = f"CCS_{tool}"
        if "percentage_difference" in df_tool.columns:
            col_map["percentage_difference"] = f"err_perc_{tool}"
        if "Predicted_mz" in df_tool.columns:
            col_map["Predicted_mz"] = f"mz_{tool}"

        df_tool = df_tool.rename(columns=col_map)

        keep_cols = list(set(id_cols + ["CCS"] + list(col_map.values())))
        df_tool = df_tool[keep_cols]

        if base_df is None:
            base_df = df_tool
        else:
            base_df = base_df.merge(df_tool, on=["SMILES", "Adduct", "CCS"], how="outer")

    if base_df is None:
        print(f"[WARN] No joined files found for dataset {dataset}.")
        return

    BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)
    if output_file is None:
        output_file = BENCHMARK_DIR / f"compounds_{dataset}.csv"

    base_df.to_csv(output_file, index=False)
    print(f"[INFO] Saved compounds table for dataset {dataset} to {output_file}")


# ---- CLI ENTRYPOINT ----

def main():
    parser = argparse.ArgumentParser(description="CCS benchmark CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # run-benchmark
    p_run = subparsers.add_parser("run-benchmark", help="Run per-dataset, per-tool benchmark (create_report).")
    p_run.add_argument("--datasets", nargs="+", default=DEFAULT_DATASETS, help="List of dataset names")
    p_run.add_argument("--tools", nargs="+", default=DEFAULT_TOOLS, help="List of tool names")

    # join-metrics
    p_join = subparsers.add_parser("join-metrics", help="Join metrics from all dataset/tool combinations.")
    p_join.add_argument("--datasets", nargs="+", default=DEFAULT_DATASETS)
    p_join.add_argument("--tools", nargs="+", default=DEFAULT_TOOLS)
    p_join.add_argument("--output", help="Output CSV for joined metrics")

    # build-compounds-table
    p_comp = subparsers.add_parser("build-compounds-table", help="Build cross-tool per-compound table for one dataset.")
    p_comp.add_argument("dataset", help="Dataset name (e.g. metlinims)")
    p_comp.add_argument("--tools", nargs="+", default=DEFAULT_TOOLS)
    p_comp.add_argument("--output", help="Output CSV file")

    args = parser.parse_args()

    if args.command == "run-benchmark":
        run_benchmark(args.datasets, args.tools)
    elif args.command == "join-metrics":
        join_metrics(args.datasets, args.tools, args.output)
    elif args.command == "build-compounds-table":
        build_compounds_table(args.dataset, args.tools, args.output)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
