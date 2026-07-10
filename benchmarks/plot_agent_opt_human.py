#!/usr/bin/env python3
"""agent-optimization on human (GRCh38 x CHM13, T=32): peak storage + wall.
All configs bit-exact (518,037 alignments, md5 8b6c42...) — matches Ashir's report."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, os
OUT = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/docs/agent_optimization"
os.makedirs(OUT, exist_ok=True)

cfg = ["baseline\n(ddeea32)", "agent nc\n(Opt1+3+4)", "agent\n-C16", "agent\n-C32"]
peak = [64.13, 55.52, 5.13, 3.47]
wall = [589, 589, 793, 963]
red  = [None, "-13%", "-92.0%", "-94.6%"]
x = np.arange(len(cfg))
fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5))

bars = a.bar(x, peak, color="#1976D2", alpha=0.88, width=0.6)
a.set_xticks(x); a.set_xticklabels(cfg, fontsize=9)
a.set_ylabel("peak scratch storage (GiB)")
a.set_title("A. Peak storage — human, all bit-exact ✓", fontsize=12)
a.set_ylim(0, 72)
for xi, v, r in zip(x, peak, red):
    a.text(xi, v+1, f"{v:.1f}", ha="center", fontsize=9, fontweight="bold")
    if r and r!="-13%": a.text(xi, v+6, r, ha="center", fontsize=10, color="#2E7D32", fontweight="bold")
a.grid(alpha=0.3, axis="y")

b.bar(x, wall, color="#D32F2F", alpha=0.88, width=0.6)
b.set_xticks(x); b.set_xticklabels(cfg, fontsize=9)
b.set_ylabel("wall time (s), T=32")
b.set_title("B. Wall time — Opt4 free on human; chunking costs time", fontsize=12)
b.set_ylim(0, 1080)
for xi, v in zip(x, wall):
    b.text(xi, v+12, f"{v}s", ha="center", fontsize=9, fontweight="bold")
b.grid(alpha=0.3, axis="y")

fig.suptitle("agent-optimization on human (GRCh38×CHM13, 3.1 Gbp, T=32)\n"
             "518,037 alignments, md5 8b6c42… — bit-exact at every config incl. -C16/-C32",
             fontsize=12.5, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.92])
p = os.path.join(OUT, "human.png"); fig.savefig(p, dpi=130); print("wrote", p)
