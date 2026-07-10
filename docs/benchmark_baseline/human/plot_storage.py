#!/usr/bin/env python3
"""Human-genome storage footprint over time for stock upstream FastGA (GRCh38 x CHM13, T=32).
Reads storage_data/timeline.tsv (elapsed_s, persistent_mb, temp_du_mb, temp_fd_mb) and
storage_data/storage.Llog (phase boundaries). Writes ../human/storage_timeline.png.

Note: temp files are open()-then-unlink()ed; on local disk `du` misses them, so temp is taken
from the process's open file descriptors (temp_fd). persistent = du(work) is exact."""
import os, re, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SD = os.path.join(HERE, "storage_data")
ORANGE, BLUE = "#f5a623", "#4a90d9"

def read_tl(p):
    t, per, tmp = [], [], []
    with open(p) as f:
        r = csv.reader(f, delimiter='\t'); next(r)
        for row in r:
            if len(row) >= 4:
                t.append(float(row[0])); per.append(float(row[1]) / 1024)  # MB->GB
                tmp.append(float(row[3]) / 1024)                            # temp_fd
    return np.array(t), np.array(per), np.array(tmp)

def _sec(tok):
    v = tok[:-1] if tok.endswith("w") else tok
    if ":" in v:
        q = v.split(":")
        return float(q[0]) * 60 + float(q[1]) if len(q) == 2 else float(q[0]) * 3600 + float(q[1]) * 60 + float(q[2])
    return float(v)

def _wall(line):
    m = re.search(r'(\d+(?::\d+)*\.\d+)w', line); return _sec(m.group(1)) if m else None

def parse_stages(path):
    sec = None; phase = []; g = {"gdb1": 0, "gix1": 0, "gdb2": 0, "gix2": 0}; seen = False
    for ln in open(path):
        if "FAtoGDB" in ln and " -" in ln: sec = "gdb1" if not seen else "gdb2"
        elif "GIXmake" in ln and " -" in ln: sec = "gix1" if not seen else "gix2"
        elif re.search(r'FastGA\b.* -T', ln): sec = "fastga"
        elif "Total Resources" in ln and sec in g:
            g[sec] = _wall(ln) or 0
            if sec == "gix1": seen = True
        elif "Resources for phase" in ln and sec == "fastga":
            w = _wall(ln)
            if w is not None: phase.append(w)
    return g["gdb1"], g["gix1"], g["gdb2"], g["gix2"], (phase[0] if phase else 0), (sum(phase[1:]) if len(phase) > 1 else 0)

t, per, tmp = read_tl(os.path.join(SD, "timeline.tsv"))
tot = per + tmp
ymax = max(tot) * 1.15
fig, ax = plt.subplots(figsize=(11, 5))
ax.set_xlim(0, t[-1]); ax.set_ylim(0, ymax)

Llog = os.path.join(SD, "storage.Llog")
if os.path.exists(Llog):
    gdb1, gix1, gdb2, gix2, merge, align = parse_stages(Llog)
    b1, b2, b3 = gdb1 + gix1, gdb1 + gix1 + gdb2 + gix2, gdb1 + gix1 + gdb2 + gix2 + merge
    for x0, x1, lab, col in [(0, b1, "GIX build\nGRCh38", "#fff3e0"), (b1, b2, "GIX build\nCHM13", "#ffe0b2"),
                             (b2, b3, "seed\nmerge", "#e3f2fd"), (b3, t[-1], "Sort + align", "#e8f5e9")]:
        x1c = min(x1, t[-1]); ax.axvspan(x0, x1c, color=col, zorder=0)
        if x1c - x0 > t[-1] * 0.03:   # only label regions wide enough (seed merge is a sliver)
            ax.text((x0 + x1c) / 2, ymax * 0.97, lab, ha="center", va="top", fontsize=8.5, color="#555", style="italic")
    ax.annotate("seed\nmerge", xy=((b2 + b3) / 2, ymax * 0.55), xytext=((b2 + b3) / 2 + t[-1] * 0.05, ymax * 0.68),
                fontsize=8, color="#1976d2", ha="center", arrowprops=dict(arrowstyle="->", color="#1976d2", lw=1))

ax.fill_between(t, 0, per, color=ORANGE, alpha=0.85, label="persistent (GIX + GDB)", zorder=2)
ax.fill_between(t, per, tot, color=BLUE, alpha=0.85, label="temp (seed pairs, alignment)", zorder=2)
ax.plot(t, tot, color="#2c6cb0", lw=1, zorder=3)
ip = int(np.argmax(tot))
ax.annotate(f"peak {tot[ip]:.1f} GB", xy=(t[ip], tot[ip]), xytext=(t[ip] + t[-1] * 0.03, tot[ip] + ymax * 0.02),
            fontsize=10, fontweight="bold", color="#c0392b", arrowprops=dict(arrowstyle="->", color="#c0392b"))
ax.set_xlabel("elapsed time (s)"); ax.set_ylabel("scratch storage (GB)")
ax.set_title("FastGA human-genome storage footprint over time — GRCh38×CHM13, T=32", fontweight="bold")
ax.legend(loc="center right", fontsize=9); ax.grid(alpha=0.25)
fig.tight_layout()
out = os.path.join(HERE, "storage_timeline.png"); fig.savefig(out, dpi=140); print("wrote", out)
