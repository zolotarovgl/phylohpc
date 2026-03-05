#!/usr/bin/env python3

import argparse
import os
import json
import numpy as np
import pandas as pd
from pathlib import Path


def get_nseq(fasta):
    with open(fasta) as f:
        return sum(1 for line in f if line.startswith(">"))


def get_mlen(fasta):
    lengths = []
    with open(fasta) as f:
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    lengths.append(len("".join(seq)))
                    seq = []
            else:
                seq.append(line)
        if seq:
            lengths.append(len("".join(seq)))
    return int(np.median(lengths)) if lengths else 0


def predict_model(coefs, nseq, mlen):
    ln_n = np.log(nseq)
    ln_m = np.log(mlen)

    log_y = (
        coefs["(Intercept)"]
        + coefs["log(nseq)"] * ln_n
        + coefs["log(mlen)"] * ln_m
        + coefs["log(nseq):log(mlen)"] * ln_n * ln_m
    )

    return np.exp(log_y)


def round_base(x, base):
    return np.where(x >= base, np.ceil(x / base) * base, x)


def convert_mem(x):
    return (np.ceil(x / 100) * 100).astype(int).astype(str) + ".MB"


def convert_time(x, min_time=5):
    return np.ceil(np.maximum(x, min_time)).astype(int).astype(str) + ".min"


def predict_resources(job_name, models, defaults, nseq, mlen):

    job_defaults = defaults[job_name]

    default_mem = job_defaults["mem"]
    default_time = job_defaults["time"]

    model = models.get(job_name)

    # ---- Memory prediction ----
    if model and "mem" in model:
        mem_val = predict_model(model["mem"], nseq, mlen)
    else:
        mem_val = default_mem

    # ---- Time prediction ----
    if model and "time" in model:
        time_val = predict_model(model["time"], nseq, mlen)
    else:
        time_val = default_time

    # enforce defaults as minimum
    mem_val = max(mem_val, default_mem)
    time_val = max(time_val, default_time)

    # ---- Large family override ----
    threshold = job_defaults.get("large_nseq_threshold")

    if threshold and nseq >= threshold:
        mem_val = job_defaults.get("large_mem", mem_val)
        time_val = job_defaults.get("large_time", time_val)

    return mem_val, time_val


def main():

    parser = argparse.ArgumentParser(description="Predict resources and write TSV")

    parser.add_argument("--ids_fn", required=True)
    parser.add_argument("--cluster_dir", required=True)

    parser.add_argument("--models_json", required=True)
    parser.add_argument("--defaults_json", required=True)

    parser.add_argument("--outfile", required=True)

    parser.add_argument("--max_mem", type=int, required=True)
    parser.add_argument("--max_time", type=int, required=True)

    parser.add_argument("--increase", type=float, required=True)

    args = parser.parse_args()

    # ---- Load IDs ----
    ids = [x.strip() for x in open(args.ids_fn)]

    fasta_files = {
        Path(f).stem: os.path.join(args.cluster_dir, f)
        for f in os.listdir(args.cluster_dir)
        if f.endswith(".fasta")
    }

    fasta_files = {k: v for k, v in fasta_files.items() if k in ids}

    print(f"{len(ids)} ids. {len(fasta_files)} .fastas found")

    # ---- Compute input features ----
    input_data = []
    for k, f in fasta_files.items():
        input_data.append({
            "id": k,
            "nseq": get_nseq(f),
            "mlen": get_mlen(f)
        })

    df_input = pd.DataFrame(input_data).set_index("id")

    # ---- Load models ----
    with open(args.models_json) as f:
        models = json.load(f)

    # ---- Load defaults ----
    with open(args.defaults_json) as f:
        defaults = json.load(f)

    jobs = sorted(defaults.keys())

    results = {}

    for job_name in jobs:

        mem_preds = []
        time_preds = []

        for idx, row in df_input.iterrows():

            nseq = row["nseq"]
            mlen = row["mlen"]

            mem_val, time_val = predict_resources(
                job_name,
                models,
                defaults,
                nseq,
                mlen
            )

            mem_preds.append(mem_val)
            time_preds.append(time_val)

        pred = pd.DataFrame({
            f"{job_name.lower()}_mem": mem_preds,
            f"{job_name.lower()}_time": time_preds
        }, index=df_input.index)

        # Apply safety increase
        mem_col = f"{job_name.lower()}_mem"
        time_col = f"{job_name.lower()}_time"

        pred[mem_col] *= (1 + args.increase)
        pred[time_col] *= (1 + args.increase)

        results[job_name.lower()] = pred

    pred = pd.concat(results.values(), axis=1)

    memcols = [c for c in pred.columns if "_mem" in c]
    timcols = [c for c in pred.columns if "_time" in c]

    # Apply maximum caps
    pred[memcols] = np.minimum(pred[memcols], args.max_mem)
    pred[timcols] = np.minimum(pred[timcols], args.max_time)

    # Round to scheduler-friendly values
    pred[memcols] = round_base(pred[memcols].values, 1024)
    pred[timcols] = round_base(pred[timcols].values, 60)

    # Convert units
    pred[memcols] = convert_mem(pred[memcols])
    pred[timcols] = convert_time(pred[timcols])

    pred.reset_index().to_csv(args.outfile, sep="\t", index=False)

    print(f"Created: {args.outfile}")


if __name__ == "__main__":
    main()