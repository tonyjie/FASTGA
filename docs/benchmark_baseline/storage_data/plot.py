#!/usr/bin/env python3
"""Baseline storage figure for current upstream FastGA (EXAMPLE).
Reads audit/ (this dir) and writes ../storage_profiling.png.
Panel A: temp-dir footprint over time (T=8). Panel B: peak total vs threads."""
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.path.join(HERE, "audit")
OUT = os.path.join(os.path.dirname(HERE), "storage_profiling.png")
C = "#1976D2"

def read_monitor(p):
    t, mb = [], []
    with open(p) as f:
        r = csv.reader(f, delimiter='\t'); next(r)
        for row in r:
            if len(row) >= 2:
                t.append(float(row[0])); mb.append(float(row[1]) / 1048576)
    return np.array(t), np.array(mb)

def read_summary(p):
    th, pk, wk, tm = [], [], [], []
    with open(p) as f:
        for row in csv.DictReader(f):
            th.append(int(row["threads"])); pk.append(float(row["peak_total_mb"]))
            wk.append(float(row["peak_workdir_mb"])); tm.append(float(row["peak_tmpdir_mb"]))
    return np.array(th), np.array(pk), np.array(wk), np.array(tm)

fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5))

t, mb = read_monitor(os.path.join(AUDIT, "T08_tmpdir_monitor.tsv"))
a.fill_between(t, mb, color=C, alpha=0.3); a.plot(t, mb, color=C, lw=2)
a.set_xlabel("elapsed time (s)"); a.set_ylabel("temp-dir disk usage (MB)")
a.set_title("A. Temp-storage timeline (T=8)", fontsize=12)
imax = int(np.argmax(mb))
a.annotate(f"seed-merge peak\n{mb[imax]:.0f} MB", xy=(t[imax], mb[imax]),
           xytext=(t[imax]+3, mb[imax]-120), fontsize=9,
           arrowprops=dict(arrowstyle="->", color="#555"))
a.set_ylim(0, max(mb)*1.25); a.grid(alpha=0.3)

th, pk, wk, tm = read_summary(os.path.join(AUDIT, "summary.csv"))
x = np.arange(len(th))
b.bar(x, wk, color="#90A4AE", label="persistent GIX+GDB", width=0.62)
b.bar(x, tm, bottom=wk, color=C, label="temp (seed pairs)", width=0.62)
b.set_xticks(x); b.set_xticklabels([f"T{t}" for t in th])
b.set_xlabel("thread count"); b.set_ylabel("peak scratch storage (MB)")
b.set_title("B. Peak storage vs threads — flat (thread-independent)", fontsize=12)
b.set_ylim(0, max(pk)*1.2)
for xi, v in zip(x, pk):
    b.text(xi, v+30, f"{v:.0f}", ha="center", fontsize=8.5, fontweight="bold")
b.legend(loc="upper right", fontsize=9); b.grid(alpha=0.3, axis="y")

fig.suptitle("Baseline storage profiling — current upstream FastGA (ddeea32), EXAMPLE (HAP1×HAP2)",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT, dpi=130); print("wrote", OUT)
