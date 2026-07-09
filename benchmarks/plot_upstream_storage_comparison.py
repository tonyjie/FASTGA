#!/usr/bin/env python3
"""Old vs New upstream FastGA storage comparison (EXAMPLE dataset).

Panel A: temp-dir disk usage over time (T=8), old (5671357) vs new (ddeea32) overlaid.
Panel B: peak total storage vs thread count, old vs new.

Reads the storage_audit_old_upstream / storage_audit_main_upstream outputs.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import os

BENCH = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/benchmarks"
OLD_DIR = os.path.join(BENCH, "storage_audit_old_upstream")
NEW_DIR = os.path.join(BENCH, "storage_audit_main_upstream")
OUT_DIR = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/docs/benchmark_storage_upstream"
os.makedirs(OUT_DIR, exist_ok=True)

OLD_LABEL = "old upstream (5671357)"
NEW_LABEL = "new upstream (ddeea32)"
OLD_COLOR = "#D32F2F"
NEW_COLOR = "#1976D2"


def read_monitor(path):
    t, mb = [], []
    with open(path) as f:
        r = csv.reader(f, delimiter='\t')
        next(r)
        for row in r:
            if len(row) >= 2:
                t.append(float(row[0]))
                mb.append(float(row[1]) / 1048576)
    return np.array(t), np.array(mb)


def read_summary(path):
    threads, peak_total = [], []
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            threads.append(int(row["threads"]))
            peak_total.append(float(row["peak_total_mb"]))
    return np.array(threads), np.array(peak_total)


fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5))

# ---- Panel A: temp timeline T=8 ----
to, mo = read_monitor(os.path.join(OLD_DIR, "T08_tmpdir_monitor.tsv"))
tn, mn = read_monitor(os.path.join(NEW_DIR, "T08_tmpdir_monitor.tsv"))
axA.plot(to, mo, color=OLD_COLOR, lw=2.2, label=OLD_LABEL, alpha=0.9)
axA.plot(tn, mn, color=NEW_COLOR, lw=1.4, ls="--", label=NEW_LABEL, alpha=0.95)
axA.set_xlabel("elapsed time (s)")
axA.set_ylabel("temp-dir disk usage (MB)")
axA.set_title("A. Temp storage timeline (T=8): curves overlap")
axA.legend(loc="upper left", fontsize=9)
axA.grid(alpha=0.3)

# ---- Panel B: peak vs threads ----
tho, po = read_summary(os.path.join(OLD_DIR, "summary.csv"))
thn, pn = read_summary(os.path.join(NEW_DIR, "summary.csv"))
x = np.arange(len(thn))
w = 0.38
axB.bar(x - w/2, po, w, color=OLD_COLOR, alpha=0.85, label=OLD_LABEL)
axB.bar(x + w/2, pn, w, color=NEW_COLOR, alpha=0.85, label=NEW_LABEL)
axB.set_xticks(x)
axB.set_xticklabels([f"T{t}" for t in thn])
axB.set_xlabel("thread count")
axB.set_ylabel("peak total storage (MB)")
axB.set_title("B. Peak storage vs threads: flat & identical")
axB.set_ylim(0, max(pn) * 1.25)
axB.legend(loc="upper right", fontsize=9)
axB.grid(alpha=0.3, axis="y")
for xi, v in zip(x, pn):
    axB.text(xi, v + max(pn)*0.02, f"{v:.0f}", ha="center", fontsize=8)

fig.suptitle("FastGA storage: old vs new upstream (EXAMPLE, HAP1 vs HAP2)",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(OUT_DIR, "storage_old_vs_new.png")
fig.savefig(out, dpi=130)
print("wrote", out)
