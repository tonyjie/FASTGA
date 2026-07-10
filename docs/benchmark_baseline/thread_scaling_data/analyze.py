#!/usr/bin/env python3
"""Thread-scaling analysis for current upstream FastGA (EXAMPLE dataset).
Reads results.tsv (this dir), prints the median table, and writes the figure
../thread_scaling.png. Single series: upstream only (no old-build comparison)."""
import csv, os, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results.tsv")
OUT = os.path.join(os.path.dirname(HERE), "thread_scaling.png")

by_t = {}
with open(RESULTS) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["nonredundant_aln"] and int(r["threads"]) <= 32:
            by_t.setdefault(int(r["threads"]), []).append(r)

T = sorted(by_t)
def med(t, k): return st.median(float(r[k]) for r in by_t[t])
wall  = {t: med(t, "wall_s") for t in T}
cpu   = {t: med(t, "cpu_pct") for t in T}
rss   = {t: med(t, "maxrss_kb") / 1024 for t in T}   # MB
speed = {t: wall[1] / wall[t] for t in T}

# ---- median table ----
print("| T | wall(s) | speedup | CPU% | RSS(MB) |")
print("|--:|--:|--:|--:|--:|")
for t in T:
    print(f"| {t} | {wall[t]:.1f} | {speed[t]:.2f}x | {cpu[t]:.0f} | {rss[t]:.0f} |")

# ---- plot (upstream only) ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
C = "#2563eb"; C_IDEAL = "#d1d5db"

ax1.plot(T, [wall[t] for t in T], "o-", color=C, lw=2, label="upstream (ddeea32)")
ax1.set_xscale("log", base=2); ax1.set_yscale("log", base=2)
ax1.set_xticks(T); ax1.set_xticklabels(T)
ax1.set_xlabel("threads (-T)"); ax1.set_ylabel("wall-clock (s, log2)")
ax1.set_title("FastGA end-to-end runtime\nEXAMPLE HAP1 vs HAP2 (~86 Mbp each)")
ax1.grid(True, which="both", alpha=0.25); ax1.legend(frameon=False)

ax2.plot(T, [speed[t] for t in T], "o-", color=C, lw=2, label="measured")
ax2.plot(T, T, ":", color=C_IDEAL, lw=1.6, label="ideal linear")
ax2.set_xscale("log", base=2); ax2.set_yscale("log", base=2)
ax2.set_xticks(T); ax2.set_xticklabels(T); ax2.set_yticks(T); ax2.set_yticklabels(T)
ax2.set_xlabel("threads (-T)"); ax2.set_ylabel("speedup vs T1 (log2)")
ax2.set_title("Parallel speedup (32-thread hard cap)")
ax2.grid(True, which="both", alpha=0.25); ax2.legend(frameon=False)

fig.tight_layout()
fig.savefig(OUT, dpi=140, bbox_inches="tight")
print("\nsaved:", OUT)
