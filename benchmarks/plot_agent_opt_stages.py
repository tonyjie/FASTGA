#!/usr/bin/env python3
"""agent-optimization incremental stages: peak storage + wall time (EXAMPLE, T=8).
All stages are bit-exact to the ddeea32 baseline (323,569 alignments)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, os

OUT_DIR = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/docs/agent_optimization"
os.makedirs(OUT_DIR, exist_ok=True)

stages = ["baseline\n(ddeea32)", "+Opt1\nearly-del", "+Opt3\nmask byte",
          "+Opt4\nLCP byte", "+OptC\n-C4", "+OptC\n-C8"]
peak = [2344, 2344, 2212, 2096, 1070, 906]
wall = [31.5, 31.4, 30.6, 65.7, 84.4, 95.5]
x = np.arange(len(stages))
C = "#1976D2"; C2 = "#D32F2F"; GREEN = "#2E7D32"

fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.2))

# Panel A: peak storage
bars = a1.bar(x, peak, color=C, alpha=0.85, width=0.6)
a1.set_xticks(x); a1.set_xticklabels(stages, fontsize=8.5)
a1.set_ylabel("peak scratch storage (MB)")
a1.set_title("A. Peak storage — every stage bit-exact ✓", fontsize=12)
a1.set_ylim(0, 2600)
for xi, v, p in zip(x, peak, [None,"0%","-5.6%","-10.6%","-54%","-61%"]):
    a1.text(xi, v+30, f"{v}", ha="center", fontsize=8.5, fontweight="bold")
    if p and p!="0%": a1.text(xi, v/2, p, ha="center", fontsize=9, color="white", fontweight="bold")
a1.grid(alpha=0.3, axis="y")

# Panel B: wall time
bars2 = a2.bar(x, wall, color=C2, alpha=0.85, width=0.6)
a2.set_xticks(x); a2.set_xticklabels(stages, fontsize=8.5)
a2.set_ylabel("wall time (s)")
a2.set_title("B. Wall time — Opt4 LCP recompute & chunking cost time", fontsize=12)
a2.set_ylim(0, 110)
for xi, v in zip(x, wall):
    a2.text(xi, v+1.5, f"{v:.0f}s", ha="center", fontsize=8.5, fontweight="bold")
a2.annotate("LCP on-the-fly\n≈2× wall", xy=(3, 65.7), xytext=(2.1, 92),
            fontsize=8.5, color=C2, arrowprops=dict(arrowstyle="->", color=C2))
a2.grid(alpha=0.3, axis="y")

fig.suptitle("agent-optimization: cumulative optimizations vs upstream (EXAMPLE, HAP1×HAP2, T=8)\n"
             "correctness: all 6 stages produce identical 323,569 alignments (bit-exact md5)",
             fontsize=12.5, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
out = os.path.join(OUT_DIR, "stages.png")
fig.savefig(out, dpi=130); print("wrote", out)
