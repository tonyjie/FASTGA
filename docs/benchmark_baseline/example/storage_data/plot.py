#!/usr/bin/env python3
"""Baseline storage figure for current upstream FastGA (EXAMPLE).
Reads audit/ (this dir) and writes ../storage_profiling.png.
Panel A: combined scratch footprint over time (T=8) — persistent GIX/GDB + temp
         seed-pairs, stacked, with FastGA phase regions marked (from T08.Llog).
Panel B: peak total (persistent + temp) vs thread count."""
import csv, os, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.path.join(HERE, "audit")
OUT = os.path.join(os.path.dirname(HERE), "storage_profiling.png")
ORANGE, BLUE = "#f5a623", "#4a90d9"

def read_timeline(p):
    t, per, tmp = [], [], []
    with open(p) as f:
        r = csv.reader(f, delimiter='\t'); next(r)
        for row in r:
            if len(row) >= 3:
                t.append(float(row[0])); per.append(float(row[1])); tmp.append(float(row[2]))
    return np.array(t), np.array(per), np.array(tmp)

def read_summary(p):
    th, per, tmp, tot = [], [], [], []
    with open(p) as f:
        for row in csv.DictReader(f):
            th.append(int(row["threads"])); per.append(float(row["peak_persistent_mb"]))
            tmp.append(float(row["peak_temp_mb"])); tot.append(float(row["peak_total_mb"]))
    return np.array(th), np.array(per), np.array(tmp), np.array(tot)

def _wsec(tok):
    v = tok[:-1] if tok.endswith("w") else tok
    if ":" in v:
        q = v.split(":")
        return float(q[0]) * 60 + float(q[1]) if len(q) == 2 else float(q[0]) * 3600 + float(q[1]) * 60 + float(q[2])
    return float(v)

def _wall(line):
    m = re.search(r'(\d+(?::\d+)*\.\d+)w', line); return _wsec(m.group(1)) if m else None

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

fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.2))

# ---- Panel A: stacked timeline + phase regions ----
t, per, tmp = read_timeline(os.path.join(AUDIT, "T08_timeline.tsv"))
tot = per + tmp
ymax = max(tot) * 1.18
a.set_ylim(0, ymax); a.set_xlim(0, t[-1])

Llog = os.path.join(AUDIT, "T08.Llog")
if os.path.exists(Llog):
    gdb1, gix1, gdb2, gix2, merge, align = parse_stages(Llog)
    b1, b2, b3 = gdb1 + gix1, gdb1 + gix1 + gdb2 + gix2, gdb1 + gix1 + gdb2 + gix2 + merge
    for x0, x1, lab, col in [(0, b1, "GIX build\nHAP1", "#fff3e0"), (b1, b2, "GIX build\nHAP2", "#ffe0b2"),
                             (b2, b3, "Seed merge", "#e3f2fd"), (b3, t[-1], "Sort + align", "#e8f5e9")]:
        x1c = min(x1, t[-1]); a.axvspan(x0, x1c, color=col, zorder=0)
        a.text((x0 + x1c) / 2, ymax * 0.96, lab, ha="center", va="top", fontsize=8.5, color="#555", style="italic")

a.fill_between(t, 0, per, color=ORANGE, alpha=0.85, label="persistent (GIX + GDB)", zorder=2)
a.fill_between(t, per, tot, color=BLUE, alpha=0.85, label="temp (seed pairs, alignment)", zorder=2)
a.plot(t, tot, color="#2c6cb0", lw=1, zorder=3)
ip = int(np.argmax(tot))
a.annotate(f"peak {tot[ip]:.0f} MB", xy=(t[ip], tot[ip]), xytext=(t[ip] + 2.5, tot[ip] + ymax * 0.02),
           fontsize=9, fontweight="bold", color="#c0392b",
           arrowprops=dict(arrowstyle="->", color="#c0392b"))
a.set_xlabel("elapsed time (s)"); a.set_ylabel("scratch storage (MB)")
a.set_title("A. Storage footprint over time (T=8), by phase", fontsize=12)
a.legend(loc="lower right", fontsize=8.5); a.grid(alpha=0.25)

# ---- Panel B: peak vs threads ----
th, sp, stmp, stot = read_summary(os.path.join(AUDIT, "summary.csv"))
x = np.arange(len(th))
a.set_ylim(0, ymax)
b.bar(x, sp, color=ORANGE, alpha=0.9, label="persistent (GIX + GDB)", width=0.62)
b.bar(x, stmp, bottom=sp, color=BLUE, alpha=0.9, label="temp (seed pairs)", width=0.62)
b.set_xticks(x); b.set_xticklabels([f"T{t}" for t in th])
b.set_xlabel("thread count"); b.set_ylabel("peak scratch storage (MB)")
b.set_title("B. Peak storage vs threads — flat (thread-independent)", fontsize=12)
b.set_ylim(0, max(stot) * 1.2)
for xi, v in zip(x, stot):
    b.text(xi, v + max(stot) * 0.015, f"{v:.0f}", ha="center", fontsize=8.5, fontweight="bold")
b.legend(loc="upper right", fontsize=8.5); b.grid(alpha=0.3, axis="y")

fig.suptitle("Baseline storage profiling — current upstream FastGA (ddeea32), EXAMPLE (HAP1×HAP2)",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT, dpi=130); print("wrote", OUT)
