#!/usr/bin/env python3
"""Incremental optimization profiling on the human pair (T=32) — comparison figures.
Reads stage_data/<stage>/{timeline.tsv, run.Llog, run.time} for each cumulative stage.

Runtime is split into three phases that are comparable across single-pass AND chunked runs:
  GDB          = sum of FAtoGDB "Total Resources" wall (serial, ~constant)
  Index+merge  = everything up to sort+align (GIXmake + seed merge; for -C this is the whole
                 chunk loop: build chunk, merge chunk, delete, x K) = end-to-end - GDB - align
  Sort+align   = the final "Sorting and merging alignments" phase (+ PAF) — untouched by the
                 storage optimizations, so ~constant across stages
Writes storage_timelines.png, runtime_breakdown.png, and a markdown table on stdout.
"""
import os, re, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SD = os.path.join(HERE, "stage_data")
STAGES = [("baseline", "upstream ddeea32"), ("opt1", "+ Opt1 early GIX deletion"),
          ("opt3", "+ Opt3 drop mask byte"), ("opt4", "+ Opt4 drop LCP byte"),
          ("optC_C16", "+ Opt C chunked (-C16)")]
PHASES = ["GDB", "Index + merge", "Sort + align"]
PC = {"GDB": "#9e9e9e", "Index + merge": "#ff9800", "Sort + align": "#2e7d32"}
ORANGE, BLUE = "#f5a623", "#4a90d9"

def _sec(t):
    t = t[:-1] if t and t[-1] in "wus" else t
    if ":" in t:
        p = t.split(":"); return float(p[0]) * 60 + float(p[1]) if len(p) == 2 else float(p[0]) * 3600 + float(p[1]) * 60 + float(p[2])
    return float(t)

def wall(line):
    m = re.search(r'([\d:.]+)w', line); return _sec(m.group(1)) if m else None

def parse(path, total):
    """Return {GDB, Index+merge, Sort+align} wall seconds."""
    sec = None; gdb = 0.0; align = 0.0; grab_align = False
    for ln in open(path):
        if "FAtoGDB" in ln and " -" in ln: sec = "gdb"
        elif "GIXmake" in ln and " -" in ln: sec = "gix"
        elif re.search(r'FastGA\b.* -T', ln): sec = "fastga"
        elif "Total Resources" in ln and sec == "gdb":
            w = wall(ln); gdb += w or 0
        elif "Sorting and merging alignments" in ln or "Converting aln" in ln:
            grab_align = True
        elif "Resources for phase" in ln and grab_align:
            w = wall(ln)
            if w: align += w
            grab_align = False
    im = max(0.0, total - gdb - align)
    return {"GDB": gdb, "Index + merge": im, "Sort + align": align}

def read_tl(path):
    t, per, tmp = [], [], []
    with open(path) as f:
        r = csv.reader(f, delimiter='\t'); next(r)
        for row in r:
            if len(row) >= 3:
                t.append(float(row[0])); per.append(float(row[1]) / 1024); tmp.append(float(row[2]) / 1024)
    return np.array(t), np.array(per), np.array(tmp)

def e2e(path):
    m = re.search(r'Elapsed .*?: *([\d:.]+)', open(path).read()); return _sec(m.group(1)) if m else 0

data = {}
for name, _ in STAGES:
    d = os.path.join(SD, name)
    if not os.path.exists(os.path.join(d, "run.Llog")): continue
    total = e2e(os.path.join(d, "run.time"))
    ph = parse(os.path.join(d, "run.Llog"), total)
    t, per, tmp = read_tl(os.path.join(d, "timeline.tsv"))
    tot = per + tmp
    ipk = int(np.argmax(tot)) if len(tot) else 0
    data[name] = {"ph": ph, "t": t, "per": per, "tmp": tmp, "wall": total,
                  "peak": float(tot[ipk]) if len(tot) else 0,
                  "peak_per": float(per[ipk]) if len(tot) else 0,
                  "peak_tmp": float(tmp[ipk]) if len(tot) else 0}
present = [(n, l) for n, l in STAGES if n in data]

# ---- table ----
print("### Incremental optimization on human (GRCh38 × CHM13, T=32)")
print("| Stage | GDB | Index+merge | Sort+align | wall (s) | peak GB (persist+temp) |")
print("|---|--:|--:|--:|--:|--:|")
for n, lab in present:
    ph = data[n]["ph"]; D = data[n]
    print(f"| {lab} | {ph['GDB']:.0f} | {ph['Index + merge']:.0f} | {ph['Sort + align']:.0f}"
          f" | {D['wall']:.0f} | {D['peak']:.1f}  ({D['peak_per']:.1f}+{D['peak_tmp']:.1f}) |")

# ---- Fig 1: storage timelines (one panel per stage) ----
N = len(present)
fig, axes = plt.subplots(N, 1, figsize=(11, 2.0 * N))
if N == 1: axes = [axes]
for ax, (n, lab) in zip(axes, present):
    t, per, tmp = data[n]["t"], data[n]["per"], data[n]["tmp"]; tot = per + tmp
    ax.fill_between(t, 0, per, color=ORANGE, alpha=0.85)
    ax.fill_between(t, per, tot, color=BLUE, alpha=0.85)
    ax.plot(t, tot, color="#2c6cb0", lw=0.8)
    ax.set_ylabel("GB")
    ax.set_ylim(0, 80 if n != "optC_C16" else ((max(tot) * 1.2) if len(tot) and max(tot) > 0 else 1))
    ax.set_title(f"{lab}   —   peak {data[n]['peak']:.1f} GB "
                 f"({data[n]['peak_per']:.1f} persistent + {data[n]['peak_tmp']:.1f} temp), wall {data[n]['wall']:.0f} s",
                 fontsize=9.5, loc="left")
    ax.grid(alpha=0.25)
axes[-1].set_xlabel("elapsed time (s)")
axes[0].fill_between([], [], color=ORANGE, alpha=0.85, label="persistent (GIX+GDB)")
axes[0].fill_between([], [], color=BLUE, alpha=0.85, label="temp (seed pairs + chunk scratch)")
axes[0].legend(loc="upper right", fontsize=8)
fig.suptitle("Storage footprint over time by optimization stage — human, T=32", fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.99]); fig.savefig(os.path.join(HERE, "storage_timelines.png"), dpi=130)

# ---- Fig 2: runtime breakdown (stacked bars) ----
fig, ax = plt.subplots(figsize=(10, 5))
x = np.arange(N); bottom = np.zeros(N)
for p in PHASES:
    vals = np.array([data[n]["ph"][p] for n, _ in present])
    ax.bar(x, vals, bottom=bottom, color=PC[p], label=p, width=0.6); bottom += vals
for xi, (n, _) in zip(x, present):
    ax.text(xi, data[n]["wall"] + 10, f"{data[n]['wall']:.0f}s", ha="center", fontsize=9, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([l for _, l in present], fontsize=8, rotation=15, ha="right")
ax.set_ylabel("wall time (s)"); ax.set_title("Runtime breakdown by optimization stage — human, T=32", fontweight="bold")
ax.legend(); ax.grid(alpha=0.3, axis="y")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "runtime_breakdown.png"), dpi=130)
print("\nsaved: storage_timelines.png, runtime_breakdown.png")
