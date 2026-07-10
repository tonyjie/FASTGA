#!/usr/bin/env python3
"""Human-genome performance breakdown for stock upstream FastGA (GRCh38 x CHM13, T=32).
Parses perf_data/logs/rep*.Llog (FastGA -L per-stage resource lines: user/sys/wall/CPU%)
and perf_data/time/rep*.time (/usr/bin/time). Prints a per-stage table (wall, share, CPU%,
median of reps) and writes ../human/perf_breakdown.png.

Stages: GDB (FAtoGDB x2, serial), GIX (GIXmake x2), Seed merge, Sort+align (+PAF)."""
import os, re, glob, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STAGES = ["GDB", "GIX", "Seed merge", "Sort+align"]
C = {"GDB": "#9e9e9e", "GIX": "#ff9800", "Seed merge": "#1976d2", "Sort+align": "#2e7d32"}

def _sec(tok):
    v = tok[:-1]
    if ":" in v:
        p = v.split(":")
        return float(p[0]) * 60 + float(p[1]) if len(p) == 2 else float(p[0]) * 3600 + float(p[1]) * 60 + float(p[2])
    return float(v)

def res(line):
    """Parse 'Xu Ys Zw P%' -> (user, sys, wall) seconds. Returns None if absent."""
    m = re.search(r'([\d:.]+)u\s+([\d:.]+)s\s+([\d:.]+)w', line)
    return (_sec(m.group(1) + "x"), _sec(m.group(2) + "x"), _sec(m.group(3) + "x")) if m else None

def parse(path):
    """Per-stage (user, sys, wall) summed. GDB/GIX from Total Resources; FastGA phases from
    the ordered 'Resources for phase' lines (merge, sort+align, paf)."""
    sec = None; seen = False; phase = []
    g = {"gdb1": None, "gix1": None, "gdb2": None, "gix2": None}
    for ln in open(path):
        if "FAtoGDB" in ln and " -" in ln: sec = "gdb1" if not seen else "gdb2"
        elif "GIXmake" in ln and " -" in ln: sec = "gix1" if not seen else "gix2"
        elif re.search(r'FastGA\b.* -T', ln): sec = "fastga"
        elif "Total Resources" in ln and sec in g:
            g[sec] = res(ln);  seen = seen or (sec == "gix1")
        elif "Resources for phase" in ln and sec == "fastga":
            r = res(ln)
            if r: phase.append(r)
    def add(*xs):
        xs = [x for x in xs if x]
        return (sum(x[0] for x in xs), sum(x[1] for x in xs), sum(x[2] for x in xs)) if xs else (0, 0, 0)
    merge = phase[0] if phase else (0, 0, 0)
    align = add(*phase[1:]) if len(phase) > 1 else (0, 0, 0)
    return {"GDB": add(g["gdb1"], g["gdb2"]), "GIX": add(g["gix1"], g["gix2"]),
            "Seed merge": merge, "Sort+align": align}

reps = sorted(glob.glob(os.path.join(HERE, "perf_data", "logs", "rep*.Llog")))
data = {s: {"wall": [], "cpu": []} for s in STAGES}
for f in reps:
    d = parse(f)
    for s in STAGES:
        u, sy, w = d[s]
        data[s]["wall"].append(w)
        data[s]["cpu"].append((u + sy) / w * 100 if w else 0)
wall = {s: st.median(data[s]["wall"]) for s in STAGES}
cpu = {s: st.median(data[s]["cpu"]) for s in STAGES}
total = sum(wall[s] for s in STAGES)

# /usr/bin/time totals (median)
tw, trss = [], []
for tf in sorted(glob.glob(os.path.join(HERE, "perf_data", "time", "rep*.time"))):
    txt = open(tf).read()
    m = re.search(r'Elapsed .*?: *([\d:.]+)', txt)
    if m: tw.append(_sec(m.group(1) + "x"))
    m = re.search(r'Maximum resident set size.*?: *(\d+)', txt)
    if m: trss.append(int(m.group(1)) / 1048576)  # GB
end2end = st.median(tw) if tw else total
rss = st.median(trss) if trss else 0

print(f"### Human performance breakdown (GRCh38 × CHM13, T=32, median of {len(reps)} reps)")
print("| Stage | wall (s) | share | CPU% | threading |")
print("|---|--:|--:|--:|---|")
thr = {"GDB": "1× (serial)", "GIX": "multi", "Seed merge": "multi", "Sort+align": "multi"}
for s in STAGES:
    print(f"| {s} | {wall[s]:.1f} | {wall[s]/total*100:.0f}% | {cpu[s]:.0f}% | {thr[s]} |")
print(f"| **Total** | {total:.1f} | 100% | — | — |")
print(f"\nend-to-end wall (/usr/bin/time): {end2end:.0f}s   peak RSS: {rss:.1f} GB")

# ---- figure: wall breakdown (left) + CPU% per stage (right) ----
fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.4))
x = np.arange(len(STAGES)); colors = [C[s] for s in STAGES]
a.bar(x, [wall[s] for s in STAGES], color=colors, width=0.6)
for xi, s in zip(x, STAGES):
    a.text(xi, wall[s] + total * 0.01, f"{wall[s]:.0f}s\n{wall[s]/total*100:.0f}%", ha="center", fontsize=8.5, fontweight="bold")
a.set_xticks(x); a.set_xticklabels(STAGES, fontsize=9); a.set_ylabel("wall time (s)")
a.set_title(f"Per-stage runtime (T=32, total {total:.0f}s)", fontsize=11); a.grid(alpha=0.3, axis="y")
b.bar(x, [cpu[s] for s in STAGES], color=colors, width=0.6)
b.axhline(3200, ls=":", color="#888", lw=1); b.text(len(STAGES)-0.5, 3250, "32 cores (3200%)", fontsize=7.5, color="#888", ha="right")
for xi, s in zip(x, STAGES):
    b.text(xi, cpu[s] + 60, f"{cpu[s]:.0f}%", ha="center", fontsize=8.5, fontweight="bold")
b.set_xticks(x); b.set_xticklabels(STAGES, fontsize=9); b.set_ylabel("CPU utilisation (%)")
b.set_title("Per-stage CPU% (100% = 1 core)", fontsize=11); b.grid(alpha=0.3, axis="y")
fig.suptitle("FastGA human-genome performance — GRCh38×CHM13 (~3.1 Gbp), T=32", fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(HERE, "perf_breakdown.png"); fig.savefig(out, dpi=140); print("wrote", out)
