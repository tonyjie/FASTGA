#!/usr/bin/env python3
"""Performance profiling for current upstream FastGA (EXAMPLE dataset).

Reads results.tsv (end-to-end metrics) + logs/T{NN}_rep{R}.Llog (FastGA -L
per-stage resource logs, this dir), prints tables, and writes three figures to
the parent dir:
  overall_scaling.png   total wall + speedup vs threads
  stage_breakdown.png   per-stage wall time, stacked, vs threads
  stage_scaling.png     per-stage speedup vs threads (each stage scales differently)

Stages: GDB (FAtoGDB x2, serial), GIX (GIXmake x2), Seed merge, Sort+align (+PAF).
"""
import csv, os, re, glob, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)

def wsec(tok):
    """'43.119w' or '1:00.184w' -> seconds (float)."""
    v = tok[:-1] if tok.endswith("w") else tok
    if ":" in v:
        p = v.split(":")
        return (float(p[0]) * 60 + float(p[1])) if len(p) == 2 else \
               (float(p[0]) * 3600 + float(p[1]) * 60 + float(p[2]))
    return float(v)

def wall_of(line):
    m = re.search(r'(\d+(?::\d+)*\.\d+)w', line)
    return wsec(m.group(1)) if m else None

def parse_Llog(path):
    """Return dict of per-stage wall seconds for one run."""
    sec = None                      # current command section
    phase = []                      # FastGA 'Resources for phase' walls, in order
    g = {"gdb1": 0, "gix1": 0, "gdb2": 0, "gix2": 0}
    seen_g1 = False
    for ln in open(path):
        if "FAtoGDB" in ln and " -" in ln:
            sec = "gdb1" if not seen_g1 else "gdb2"
        elif "GIXmake" in ln and " -" in ln:
            sec = "gix1" if not seen_g1 else "gix2"; seen_g1 = True if sec == "gix1" else seen_g1
        elif re.search(r'FastGA\b.* -T', ln):
            sec = "fastga"
        elif "Total Resources" in ln and sec in ("gdb1", "gix1", "gdb2", "gix2"):
            g[sec] = wall_of(ln) or 0
            if sec == "gix1": seen_g1 = True
        elif "Resources for phase" in ln and sec == "fastga":
            w = wall_of(ln)
            if w is not None: phase.append(w)
    merge = phase[0] if len(phase) > 0 else 0
    align = (phase[1] if len(phase) > 1 else 0) + (phase[2] if len(phase) > 2 else 0)  # sort+align (+PAF)
    return {"GDB": g["gdb1"] + g["gdb2"], "GIX": g["gix1"] + g["gix2"],
            "Seed merge": merge, "Sort+align": align}

STAGES = ["GDB", "GIX", "Seed merge", "Sort+align"]

# ---- per-stage timing from -L logs (median across reps) ----
by_t = {}
for f in glob.glob(os.path.join(HERE, "logs", "T*_rep*.Llog")):
    T = int(re.search(r'T(\d+)_', os.path.basename(f)).group(1))
    if T > 32: continue
    by_t.setdefault(T, []).append(parse_Llog(f))
T = sorted(by_t)
stage = {s: {t: st.median(r[s] for r in by_t[t]) for t in T} for s in STAGES}
tot_stage = {t: sum(stage[s][t] for s in STAGES) for t in T}

# ---- end-to-end totals from results.tsv ----
tt = {}
with open(os.path.join(HERE, "results.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if int(r["threads"]) <= 32 and r["wall_s"]:
            tt.setdefault(int(r["threads"]), []).append(float(r["wall_s"]))
wall = {t: st.median(tt[t]) for t in sorted(tt)}
speed = {t: wall[T[0]] / wall[t] for t in T}

# ---- merged table: per-stage wall (all T) + speedup + share ----
t0, tN = T[0], T[-1]
print("### Per-stage runtime (s), speedup, and share of total (median)")
print("| Stage | " + " | ".join(f"T={t}" for t in T) + f" | speedup@{tN} | share T={t0}→{tN} |")
print("|---|" + "--:|" * len(T) + "--:|--:|")
def sh(s, t): return stage[s][t] / tot_stage[t] * 100
for s in STAGES:
    walls = " | ".join(f"{stage[s][t]:.1f}" for t in T)
    sp = f"{stage[s][t0]/stage[s][tN]:.2f}x" if stage[s][tN] else "—"
    print(f"| {s} | {walls} | {sp} | {sh(s,t0):.0f}% → {sh(s,tN):.0f}% |")
print(f"| **Total (s)** | " + " | ".join(f"{tot_stage[t]:.1f}" for t in T) +
      f" | {tot_stage[t0]/tot_stage[tN]:.2f}x | 100% |")

C = {"GDB": "#9e9e9e", "GIX": "#ff9800", "Seed merge": "#1976d2", "Sort+align": "#2e7d32"}

# ---- Fig 1: overall scaling ----
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8))
a1.plot(T, [wall[t] for t in T], "o-", color="#1976d2", lw=2, label="upstream (ddeea32)")
a1.set_xscale("log", base=2); a1.set_yscale("log", base=2); a1.set_xticks(T); a1.set_xticklabels(T)
a1.set_xlabel("threads (-T)"); a1.set_ylabel("wall-clock (s, log2)")
a1.set_title("Overall end-to-end runtime"); a1.grid(True, which="both", alpha=0.25); a1.legend(frameon=False)
a2.plot(T, [speed[t] for t in T], "o-", color="#1976d2", lw=2, label="measured")
a2.plot(T, T, ":", color="#cccccc", lw=1.6, label="ideal linear")
a2.set_xscale("log", base=2); a2.set_yscale("log", base=2); a2.set_xticks(T); a2.set_xticklabels(T)
a2.set_yticks(T); a2.set_yticklabels(T)
a2.set_xlabel("threads (-T)"); a2.set_ylabel("speedup vs T1"); a2.set_title("Overall parallel speedup")
a2.grid(True, which="both", alpha=0.25); a2.legend(frameon=False)
fig.suptitle("FastGA performance — overall (EXAMPLE HAP1×HAP2)", fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.94]); fig.savefig(os.path.join(PARENT, "overall_scaling.png"), dpi=140)

# ---- Fig 2: stage breakdown (stacked) ----
fig, ax = plt.subplots(figsize=(9, 5))
x = np.arange(len(T)); bottom = np.zeros(len(T))
for s in STAGES:
    vals = np.array([stage[s][t] for t in T])
    ax.bar(x, vals, bottom=bottom, color=C[s], label=s, width=0.6)
    bottom += vals
for xi, t in zip(x, T):
    ax.text(xi, tot_stage[t] + tot_stage[T[0]] * 0.01, f"{tot_stage[t]:.0f}s", ha="center", fontsize=8.5, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([f"T{t}" for t in T])
ax.set_xlabel("threads (-T)"); ax.set_ylabel("wall time (s)")
ax.set_title("Per-stage time breakdown (EXAMPLE, median of reps)")
ax.legend(); ax.grid(alpha=0.3, axis="y")
fig.tight_layout(); fig.savefig(os.path.join(PARENT, "stage_breakdown.png"), dpi=140)

# ---- Fig 3: per-stage scaling ----
fig, ax = plt.subplots(figsize=(9, 5))
for s in STAGES:
    sp = [stage[s][T[0]] / stage[s][t] if stage[s][t] else np.nan for t in T]
    ax.plot(T, sp, "o-", color=C[s], lw=2, label=s)
ax.plot(T, T, ":", color="#cccccc", lw=1.6, label="ideal linear")
ax.set_xscale("log", base=2); ax.set_yscale("log", base=2); ax.set_xticks(T); ax.set_xticklabels(T)
ax.set_yticks(T); ax.set_yticklabels(T)
ax.set_xlabel("threads (-T)"); ax.set_ylabel("speedup vs T1")
ax.set_title("Per-stage thread scaling — each stage scales differently")
ax.legend(); ax.grid(True, which="both", alpha=0.25)
fig.tight_layout(); fig.savefig(os.path.join(PARENT, "stage_scaling.png"), dpi=140)
print("\nsaved: overall_scaling.png, stage_breakdown.png, stage_scaling.png")
