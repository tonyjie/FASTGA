#!/usr/bin/env python3
"""Aggregate FastGA -L per-phase logs across a dataset matrix.
Parser reproduces the human baseline breakdown (docs/benchmark_baseline/human)."""
import re

STAGES = ["GDB", "GIX", "Seed merge", "Sort+align"]

def sec(tok):
    """'6.706' | '2:06.205' | '8:14.618' -> seconds."""
    v = tok.split(":")
    if len(v) == 1:
        return float(v[0])
    if len(v) == 2:
        return float(v[0]) * 60 + float(v[1])
    return float(v[0]) * 3600 + float(v[1]) * 60 + float(v[2])

def res(line):
    """'... Xu Ys Zw P%' -> (user, sys, wall) or None."""
    m = re.search(r'([\d:.]+)u\s+([\d:.]+)s\s+([\d:.]+)w', line)
    return (sec(m.group(1)), sec(m.group(2)), sec(m.group(3))) if m else None

def _cpu(line):
    m = re.search(r'([\d.]+)%', line)
    return float(m.group(1)) if m else None

def parse_llog(path):
    """Per-phase wall seconds. GDB/GIX summed over the two genomes; FastGA phases from the
    last FastGA invocation. Also returns cpu%, peak RSS (MB), n_aln, ave_len."""
    out = {s: 0.0 for s in STAGES}
    cpu = {s: [] for s in STAGES}
    rss_mb = 0.0; n_aln = 0; ave_len = 0
    lines = open(path).read().splitlines()

    # First pass: find positions of all FastGA invocations and prep commands
    fastga_indices = []
    gdb_indices = []
    gix_indices = []
    for i, ln in enumerate(lines):
        if re.search(r'FastGA\b.* -T', ln):
            fastga_indices.append(i)
        elif re.search(r'FAtoGDB\b.* -', ln):
            gdb_indices.append(i)
        elif re.search(r'GIXmake\b.* -', ln):
            gix_indices.append(i)

    # Determine which prep cycle to use: the single FAtoGDB x2 + GIXmake x2
    # cycle that immediately precedes the FINAL FastGA invocation. Both scans
    # are bounded to indices before last_fastga_idx so a stray GDB/GIX header
    # occurring after the final FastGA run (if one ever existed) is ignored.
    last_fastga_idx = fastga_indices[-1] if fastga_indices else float('inf')
    gdb_before = [i for i in gdb_indices if i < last_fastga_idx]
    gix_before = [i for i in gix_indices if i < last_fastga_idx]
    # The last prep pair consists of the last 2 GDBs before last_fastga_idx, so
    # its start is the 2nd-to-last GDB header. Only set this if there is also
    # at least one GIXmake before last_fastga_idx (i.e. a prep cycle exists).
    last_prep_gdb_start = None
    if gix_before and len(gdb_before) >= 2:
        last_prep_gdb_start = gdb_before[-2]

    # Second pass: accumulate GDB/GIX from the target range
    sect = None
    fastga_phase = 0
    gdb_total = 0.0
    gix_total = 0.0

    for i, ln in enumerate(lines):
        if i >= last_fastga_idx:
            # Handle FastGA section
            if re.search(r'FastGA\b.* -T', ln):
                sect = "fastga"
                fastga_phase = 0
                out["GDB"] = gdb_total
                out["GIX"] = gix_total
            elif "Total Resources" in ln and sect == "fastga":
                r = res(ln); c = _cpu(ln)
                m = re.search(r'([\d.]+)MB', ln)
                if m: rss_mb = float(m.group(1))
            elif "Resources for phase" in ln and sect == "fastga":
                r = res(ln); c = _cpu(ln)
                if r is None: continue
                if fastga_phase == 0:
                    out["Seed merge"] = r[2]; cpu["Seed merge"] = [c]; fastga_phase = 1
                elif fastga_phase == 1:
                    out["Sort+align"] = r[2]; cpu["Sort+align"] = [c]; fastga_phase = 2
            elif "non-redundant aln" in ln:
                m = re.search(r'([\d]+)\s+non-redundant aln.s of ave len\s+([\d]+)', ln)
                if m: n_aln = int(m.group(1)); ave_len = int(m.group(2))
        elif last_prep_gdb_start is not None and i >= last_prep_gdb_start and i < last_fastga_idx:
            # Handle prep section (only the last 2 GDB+GIX pairs). An earlier,
            # non-final FastGA invocation can fall inside this index window
            # (e.g. repeated FastGA runs sharing one prep cycle); reset sect
            # so that block's own "Total Resources" (peak RSS, not GDB/GIX
            # wall time) is never attributed to GDB/GIX.
            if re.search(r'FastGA\b.* -T', ln):
                sect = None
            elif re.search(r'FAtoGDB\b.* -', ln):
                sect = "gdb"
            elif re.search(r'GIXmake\b.* -', ln):
                sect = "gix"
            elif "Total Resources" in ln:
                r = res(ln); c = _cpu(ln)
                if sect == "gdb" and r:
                    gdb_total += r[2]
                    cpu["GDB"].append(c)
                elif sect == "gix" and r:
                    gix_total += r[2]
                    cpu["GIX"].append(c)

    out["cpu"] = {s: (sum(v)/len(v) if v and all(x is not None for x in v) else None)
                  for s, v in cpu.items()}
    out["rss_mb"] = rss_mb; out["n_aln"] = n_aln; out["ave_len"] = ave_len
    return out

import os, glob, statistics as _st, sys

def point_median(logs_dir):
    reps = [parse_llog(p) for p in sorted(glob.glob(os.path.join(logs_dir, "rep*.Llog")))]
    if not reps:
        raise FileNotFoundError(f"no rep*.Llog in {logs_dir}")
    med = {}
    for s in STAGES:
        med[s] = _st.median(r[s] for r in reps)
    med["rss_mb"] = _st.median(r["rss_mb"] for r in reps)
    med["n_aln"]  = int(_st.median(r["n_aln"] for r in reps))
    med["ave_len"] = int(_st.median(r["ave_len"] for r in reps))
    med["total"]  = sum(med[s] for s in STAGES)
    return med

def aggregate(base_dir, points):
    rows = []
    for label, rank in sorted(points, key=lambda x: x[1]):
        m = point_median(os.path.join(base_dir, label, "logs"))
        rows.append((label, rank, m))
    tsv = os.path.join(base_dir, "results.tsv")
    with open(tsv, "w") as f:
        f.write("label\trank\t" + "\t".join(STAGES) +
                "\ttotal_s\tsort_align_share\trss_mb\tn_aln\tave_len\n")
        for label, rank, m in rows:
            f.write(f"{label}\t{rank}\t" + "\t".join(f"{m[s]:.1f}" for s in STAGES) +
                    f"\t{m['total']:.1f}\t{m['Sort+align']/m['total']:.3f}"
                    f"\t{m['rss_mb']:.0f}\t{m['n_aln']}\t{m['ave_len']}\n")
    _plot(base_dir, rows)
    return tsv

def _plot(base_dir, rows):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, numpy as np
    C = {"GDB": "#9e9e9e", "GIX": "#ff9800", "Seed merge": "#1976d2", "Sort+align": "#2e7d32"}
    labels = [r[0] for r in rows]
    x = np.arange(len(rows)); bottom = np.zeros(len(rows))
    fig, ax = plt.subplots(figsize=(1.6 * len(rows) + 3, 4.5))
    for s in STAGES:
        share = np.array([r[2][s] / r[2]["total"] for r in rows])
        ax.bar(x, share, bottom=bottom, color=C[s], label=s, width=0.62)
        bottom += share
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("share of runtime"); ax.set_ylim(0, 1)
    ax.set_title("FastGA per-phase share vs divergence (most similar → most divergent)")
    ax.legend(loc="upper right", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(os.path.join(base_dir, "divergence_phase_share.png"), dpi=140)

if __name__ == "__main__":
    base = sys.argv[1]
    pts = [(a.split(":")[0], int(a.split(":")[1])) for a in sys.argv[2:]]
    print("wrote", aggregate(base, pts))
