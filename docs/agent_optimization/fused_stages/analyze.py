#!/usr/bin/env python3
"""Three-way comparison: baseline (ddeea32, no -C) vs Opt C (agent-optimization, per-chunk
re-scan) vs fused-B (fused-scan-once, -R scan-once + tmpfs ktab).

Reads stage_data/<dataset>/T<T>/<cfg>/{timeline.tsv, run.Llog, run.time, md5.txt, count.txt}
for each config at each thread count and prints a markdown table + saves a figure.

Runtime is split into three phases, comparable across single-pass AND chunked runs (same
technique as human_stages/analyze.py):
  GDB          = sum of FAtoGDB "Total Resources" wall (serial, ~constant)
  Index+merge  = end-to-end - GDB - align (for -C configs this is the whole chunk loop:
                 build chunk, merge chunk, delete, x N; for fusedB the per-chunk GIXmake
                 calls carry -R so they skip re-scanning)
  Sort+align   = the final "Sorting and merging alignments" phase (+ PAF) -- untouched by
                 the scan-once/tmpfs changes, so ~constant across configs

Usage: analyze.py <dataset>   (dataset = EXAMPLE | human)
"""
import os, re, sys, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATASET = sys.argv[1] if len(sys.argv) > 1 else "EXAMPLE"
SD = os.path.join(HERE, "stage_data", DATASET)

CONFIGS = [("baseline", "baseline (ddeea32, no -C)"),
           ("optC_C4", "Opt C -C4 (re-scan/chunk)"),
           ("optC_C8", "Opt C -C8 (re-scan/chunk)"),
           ("fusedB_C4", "fused-B -C4 (-R scan-once, tmpfs)"),
           ("fusedB_C8", "fused-B -C8 (-R scan-once, tmpfs)")]
PHASES = ["GDB", "Index + merge", "Sort + align"]
PC = {"GDB": "#9e9e9e", "Index + merge": "#ff9800", "Sort + align": "#2e7d32"}
REAL, TMPFS, RSS_C = "#4a90d9", "#8e5ae2", "#e05c5c"


def _sec(t):
    t = t[:-1] if t and t[-1] in "wus" else t
    if ":" in t:
        p = t.split(":")
        return float(p[0]) * 60 + float(p[1]) if len(p) == 2 else float(p[0]) * 3600 + float(p[1]) * 60 + float(p[2])
    return float(t)


def wall(line):
    m = re.search(r'([\d:.]+)w', line)
    return _sec(m.group(1)) if m else None


def parse_llog(path, total):
    """Return {GDB, Index+merge, Sort+align} wall seconds from a FastGA -L log."""
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
    t, real, tmpfs, rss = [], [], [], []
    with open(path) as f:
        r = csv.reader(f, delimiter='\t'); next(r)
        for row in r:
            if len(row) >= 4:
                t.append(float(row[0])); real.append(float(row[1]) / 1024)
                tmpfs.append(float(row[2]) / 1024); rss.append(float(row[3]) / 1024)
    return np.array(t), np.array(real), np.array(tmpfs), np.array(rss)


def e2e(path):
    m = re.search(r'Elapsed .*?: *([\d:.]+)', open(path).read())
    return _sec(m.group(1)) if m else 0


def peak_rss_gb(path):
    m = re.search(r'Maximum resident set size \(kbytes\): *(\d+)', open(path).read())
    return int(m.group(1)) / 1024 / 1024 if m else 0.0


def read_txt(path, default=""):
    try:
        return open(path).read().strip()
    except FileNotFoundError:
        return default


Ts = sorted(int(d[1:]) for d in os.listdir(SD) if d.startswith("T") and
            os.path.isdir(os.path.join(SD, d))) if os.path.isdir(SD) else []

data = {}
for T in Ts:
    for name, _ in CONFIGS:
        d = os.path.join(SD, f"T{T}", name)
        if not os.path.exists(os.path.join(d, "run.Llog")):
            continue
        total = e2e(os.path.join(d, "run.time"))
        ph = parse_llog(os.path.join(d, "run.Llog"), total)
        t, real, tmpfs, rss = read_tl(os.path.join(d, "timeline.tsv"))
        tot = real + tmpfs
        ipk = int(np.argmax(tot)) if len(tot) else 0
        data[(T, name)] = {
            "ph": ph, "t": t, "real": real, "tmpfs": tmpfs, "rss": rss, "wall": total,
            "peak_real": float(real[ipk]) if len(tot) else 0.0,
            "peak_tmpfs": float(tmpfs[ipk]) if len(tot) else 0.0,
            "peak_rss_tl": float(rss.max()) if len(rss) else 0.0,
            "peak_rss": peak_rss_gb(os.path.join(d, "run.time")),
            "md5": read_txt(os.path.join(d, "md5.txt"), "MISSING")[:12],
            "count": read_txt(os.path.join(d, "count.txt"), "MISSING"),
        }

present = [(T, n, lab) for T in Ts for n, lab in CONFIGS if (T, n) in data]
if not present:
    print(f"no stage_data found under {SD} -- run run_fused.sh {DATASET} <T> first")
    sys.exit(0)

# ---- table ----
print(f"### fused-B (scan-once + tmpfs) vs baseline / Opt C -- {DATASET}\n")
ref_md5 = {}
for T in Ts:
    if (T, "baseline") in data:
        ref_md5[T] = data[(T, "baseline")]["md5"]

print("| T | Config | GDB (s) | Index+merge (s) | Sort+align (s) | wall (s) | real-disk peak (GB) | tmpfs peak (GB) | peak RSS (GB) | md5 | count | bit-exact |")
print("|--:|---|--:|--:|--:|--:|--:|--:|--:|---|--:|:--:|")
for T, n, lab in present:
    D = data[(T, n)]; ph = D["ph"]
    exact = "-" if n == "baseline" else ("YES" if D["md5"] == ref_md5.get(T) and D["count"] == data.get((T, "baseline"), {}).get("count") else "NO")
    print(f"| {T} | {lab} | {ph['GDB']:.0f} | {ph['Index + merge']:.0f} | {ph['Sort + align']:.0f}"
          f" | {D['wall']:.0f} | {D['peak_real']:.2f} | {D['peak_tmpfs']:.2f} | {D['peak_rss']:.2f}"
          f" | {D['md5']} | {D['count']} | {exact} |")

# ---- key checks ----
print("\n### Gating checks")
for T in Ts:
    if (T, "baseline") not in data:
        continue
    base = data[(T, "baseline")]
    print(f"\n**T={T}**")
    for n in ("optC_C4", "optC_C8", "fusedB_C4", "fusedB_C8"):
        if (T, n) not in data:
            continue
        D = data[(T, n)]
        ok = "OK" if D["md5"] == base["md5"] and D["count"] == base["count"] else "MISMATCH"
        print(f"- {n}: bit-exact vs baseline: **{ok}** (md5 {D['md5']} vs {base['md5']}, count {D['count']} vs {base['count']})")
    for cN, fN in (("4", "4"), ("8", "8")):
        oc, fb = (T, f"optC_C{cN}"), (T, f"fusedB_C{fN}")
        if oc in data and fb in data:
            im_oc, im_fb = data[oc]["ph"]["Index + merge"], data[fb]["ph"]["Index + merge"]
            rel = "<=" if im_fb <= im_oc else ">"
            print(f"- Index+merge: fusedB_C{fN} ({im_fb:.0f}s) {rel} optC_C{cN} ({im_oc:.0f}s)")
        if fb in data:
            print(f"- fusedB_C{fN} real-disk peak: {data[fb]['peak_real']:.2f} GB (tmpfs peak {data[fb]['peak_tmpfs']:.2f} GB)")

# ---- figure: runtime breakdown + footprint, one row per T ----
N = len(Ts)
if N:
    fig, axes = plt.subplots(N, 2, figsize=(14, 4.2 * N), squeeze=False)
    for row, T in enumerate(Ts):
        cfgs = [(n, lab) for n, lab in CONFIGS if (T, n) in data]
        x = np.arange(len(cfgs)); bottom = np.zeros(len(cfgs))
        axL = axes[row][0]
        for p in PHASES:
            vals = np.array([data[(T, n)]["ph"][p] for n, _ in cfgs])
            axL.bar(x, vals, bottom=bottom, color=PC[p], label=p, width=0.6); bottom += vals
        for xi, (n, _) in zip(x, cfgs):
            axL.text(xi, data[(T, n)]["wall"] + 1, f"{data[(T, n)]['wall']:.0f}s", ha="center", fontsize=8, fontweight="bold")
        axL.set_xticks(x); axL.set_xticklabels([lab for _, lab in cfgs], fontsize=7.5, rotation=20, ha="right")
        axL.set_ylabel("wall time (s)"); axL.set_title(f"T={T}: runtime breakdown", fontsize=10)
        axL.grid(alpha=0.3, axis="y")
        if row == 0: axL.legend(fontsize=8)

        axR = axes[row][1]
        w = 0.35
        real_v = np.array([data[(T, n)]["peak_real"] for n, _ in cfgs])
        tmpfs_v = np.array([data[(T, n)]["peak_tmpfs"] for n, _ in cfgs])
        rss_v = np.array([data[(T, n)]["peak_rss"] for n, _ in cfgs])
        axR.bar(x - w/2, real_v, width=w, color=REAL, label="real-disk peak")
        axR.bar(x - w/2, tmpfs_v, width=w, bottom=real_v, color=TMPFS, label="tmpfs peak")
        axR.bar(x + w/2, rss_v, width=w, color=RSS_C, label="peak RSS")
        axR.set_xticks(x); axR.set_xticklabels([lab for _, lab in cfgs], fontsize=7.5, rotation=20, ha="right")
        axR.set_ylabel("GB"); axR.set_title(f"T={T}: footprint (disk vs tmpfs vs RSS)", fontsize=10)
        axR.grid(alpha=0.3, axis="y")
        if row == 0: axR.legend(fontsize=8)
    fig.suptitle(f"fused-B (scan-once + tmpfs) vs baseline / Opt C -- {DATASET}", fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    outpng = os.path.join(HERE, f"fused_comparison_{DATASET}.png")
    fig.savefig(outpng, dpi=130)
    print(f"\nsaved: {os.path.relpath(outpng, HERE)}")
