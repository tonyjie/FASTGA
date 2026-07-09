#!/usr/bin/env python3
"""optimize-memory (rebased onto ddeea32) vs new upstream: storage effects.

Panel A: combined footprint (GIX/GDB + temp) over time, no -k. Shows Opt1 early
         GIX deletion collapsing the sort+align tail.
Panel B: GIX (ktab) size, upstream vs optimized. Shows Opt3 mask-byte removal.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv, os

BENCH = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/benchmarks"
TL = os.path.join(BENCH, "opt1_timeline")
OUT_DIR = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/docs/benchmark_storage_upstream"
os.makedirs(OUT_DIR, exist_ok=True)

UP_C, OPT_C = "#D32F2F", "#1976D2"

def load(p):
    t, m = [], []
    with open(p) as f:
        r = csv.reader(f, delimiter='\t'); next(r)
        for row in r:
            if len(row) >= 2:
                t.append(float(row[0])); m.append(float(row[1]))
    return np.array(t), np.array(m)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={'width_ratios':[2,1]})

# Panel A: timeline
tu, mu = load(os.path.join(TL, "upstream_timeline.tsv"))
to, mo = load(os.path.join(TL, "optimized_timeline.tsv"))
axA.fill_between(tu, mu, color=UP_C, alpha=0.25)
axA.plot(tu, mu, color=UP_C, lw=2.0, label="new upstream (keeps GIX to exit)")
axA.fill_between(to, mo, color=OPT_C, alpha=0.25)
axA.plot(to, mo, color=OPT_C, lw=2.0, label="optimize-memory (Opt1: early GIX delete)")
# annotate the seed-merge peak and the sort+align tail
axA.axvline(14.3, color="gray", ls=":", lw=1)
axA.text(14.3, 2400, " seed-merge peak\n (GIX+temp coexist)", fontsize=8, va="top")
axA.annotate("sort+align tail\n2003 → 146 MB (-93%)",
             xy=(24, 300), xytext=(17, 1300), fontsize=9, color=OPT_C,
             arrowprops=dict(arrowstyle="->", color=OPT_C))
axA.set_xlabel("elapsed time (s)")
axA.set_ylabel("combined disk footprint (MB)")
axA.set_title("A. Opt1 — footprint over time (EXAMPLE, T=8, no -k)")
axA.legend(loc="upper right", fontsize=8.5)
axA.grid(alpha=0.3)

# Panel B: Opt3 ktab size (HAP1)
up_ktab = 840595878 / 1048576
opt_ktab = 775934664 / 1048576
bars = axB.bar(["upstream", "optimized"], [up_ktab, opt_ktab],
               color=[UP_C, OPT_C], alpha=0.85, width=0.55)
axB.set_ylabel("GIX ktab size, HAP1 (MB)")
axB.set_title("B. Opt3 — mask-byte removal")
axB.set_ylim(0, up_ktab*1.15)
for b, v in zip(bars, [up_ktab, opt_ktab]):
    axB.text(b.get_x()+b.get_width()/2, v+8, f"{v:.0f}", ha="center", fontsize=9)
axB.annotate("-7.69%", xy=(1, opt_ktab), xytext=(0.5, up_ktab*0.5),
             fontsize=11, color=OPT_C, fontweight="bold", ha="center")
axB.grid(alpha=0.3, axis="y")

fig.suptitle("optimize-memory (rebased @ ddeea32) vs new upstream — storage",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.96])
out = os.path.join(OUT_DIR, "optmem_vs_upstream.png")
fig.savefig(out, dpi=130)
print("wrote", out)
