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

    # Determine which prep to use: the one immediately before the last FastGA
    last_fastga_idx = fastga_indices[-1] if fastga_indices else float('inf')
    # Find the last GIXmake before last FastGA (it marks the end of the last prep pair)
    last_prep_gix = None
    for gix_idx in reversed(gix_indices):
        if gix_idx < last_fastga_idx:
            last_prep_gix = gix_idx
            break
    # Find the start of this prep pair (the last GDB at least one pair back)
    last_prep_gdb_start = None
    if last_prep_gix is not None:
        # The last prep pair consists of the last 2 GDBs before last_fastga_idx
        # So we want to include lines starting from the 3rd-to-last GDB
        gdb_count_from_end = 0
        for gdb_idx in reversed(gdb_indices):
            gdb_count_from_end += 1
            if gdb_count_from_end == 2:  # The 3rd GDB from the end marks the start of last pair
                last_prep_gdb_start = gdb_idx
                break

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
            # Handle prep section (only the last 2 GDB+GIX pairs)
            if re.search(r'FAtoGDB\b.* -', ln):
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
