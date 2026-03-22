#!/usr/bin/env python3
"""Compare storage timeline: baseline vs early GIX deletion (T=32, EXAMPLE dataset)."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import csv
import os

BENCH_DIR = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/benchmarks"
AUDIT_DIR = os.path.join(BENCH_DIR, "storage_audit_v2")
DOCS_DIR = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/docs"

# Phase timing for T=32 (from previous benchmark)
PHASES = [
    (0, 0.5, "GDB1", "#E8E8E8"),
    (0.5, 4.3, "GIX1", "#FFE0B2"),
    (4.3, 4.8, "GDB2", "#E8E8E8"),
    (4.8, 8.7, "GIX2", "#FFE0B2"),
    (8.7, 14.4, "Seed Merge", "#BBDEFB"),
    (14.4, 21.5, "Sort + Align", "#C8E6C9"),
]

GIX_PER_GENOME = 930  # MB
GDB_PER_GENOME = 21   # MB


def read_monitor(filepath):
    times, mbs = [], []
    with open(filepath) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for row in reader:
            if len(row) >= 2:
                t = float(row[0])
                b = float(row[1])
                times.append(t)
                mbs.append(b / 1048576)
    return np.array(times), np.array(mbs)


def build_persistent_model(t_fine, keep_gix=True):
    """Model persistent file sizes over time."""
    persistent = np.zeros_like(t_fine)
    for i, t in enumerate(t_fine):
        p = 0.0
        if t >= 0.5:
            # GDB1 done
            p += GDB_PER_GENOME
            # GIX1 building from 0.5 to 4.3
            frac = min(1.0, max(0, (t - 0.5) / 3.8))
            if keep_gix or t < 14.4:  # GIX deleted at ~14.4s (after seed merge)
                p += GIX_PER_GENOME * frac
            elif t >= 14.4:
                p += 0  # GIX deleted
        if t >= 4.8:
            p += GDB_PER_GENOME
            # GIX2 building from 4.8 to 8.7
            frac = min(1.0, max(0, (t - 4.8) / 3.9))
            if keep_gix or t < 14.4:
                p += GIX_PER_GENOME * frac
            elif t >= 14.4:
                p += 0
        persistent[i] = p
    return persistent


def plot_panel(ax, times, temp_mb, persistent, phases, title, ylim):
    # Phase backgrounds
    for start, end, label, color in phases:
        ax.axvspan(start, end, alpha=0.15, color=color, zorder=0)
        ax.axvline(start, color='gray', linewidth=0.5, linestyle=':', alpha=0.5)
        mid = (start + end) / 2
        if end - start > 2:
            ax.text(mid, ylim * 0.93, label, ha='center', va='center',
                    fontsize=7, color='#555', style='italic')

    # Stacked areas
    ax.fill_between(times, persistent, alpha=0.5, color='#FFA726', zorder=1)
    ax.fill_between(times, persistent, persistent + temp_mb, alpha=0.5, color='#42A5F5', zorder=2)
    ax.plot(times, persistent + temp_mb, color='#1565C0', linewidth=0.8, zorder=3)
    ax.plot(times, persistent, color='#E65100', linewidth=0.8, zorder=3, linestyle='--')

    # Peak annotation
    total = persistent + temp_mb
    peak_idx = np.argmax(total)
    peak_val = total[peak_idx]
    peak_time = times[peak_idx]

    ax.annotate(f'Peak: {peak_val:.0f} MB',
                xy=(peak_time, peak_val),
                xytext=(peak_time + 3, min(peak_val + 200, ylim * 0.85)),
                arrowprops=dict(arrowstyle='->', color='red', lw=1.2),
                fontsize=9, color='red', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='red', alpha=0.9))

    ax.set_ylabel('Disk Usage (MB)', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold', loc='left')
    ax.set_xlim(0, max(times) * 1.05)
    ax.set_ylim(0, ylim)
    ax.tick_params(labelsize=8)
    ax.grid(axis='y', alpha=0.3)


def main():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # Shared y-axis limit for fair comparison
    ylim = 2800

    # --- Panel 1: Baseline (with -k, GIX kept) ---
    baseline_file = os.path.join(AUDIT_DIR, "T32_baseline_tmpdir_monitor.tsv")
    times_b, temp_b = read_monitor(baseline_file)
    persistent_b = build_persistent_model(times_b, keep_gix=True)
    plot_panel(ax1, times_b, temp_b, persistent_b, PHASES,
               "Baseline (original FastGA, T=32) — GIX kept throughout", ylim)

    # --- Panel 2: Early GIX deletion (without -k) ---
    early_file = os.path.join(AUDIT_DIR, "T32_early_gix_tmpdir_monitor.tsv")
    times_e, temp_e = read_monitor(early_file)
    persistent_e = build_persistent_model(times_e, keep_gix=False)
    plot_panel(ax2, times_e, temp_e, persistent_e, PHASES,
               "With Early GIX Deletion (T=32) — GIX freed after seed merge", ylim)

    # Add "GIX deleted here" annotation on panel 2
    ax2.annotate('GIX deleted\n(~1,860 MB freed)',
                 xy=(14.4, 42), xytext=(17, 800),
                 arrowprops=dict(arrowstyle='->', color='green', lw=2),
                 fontsize=10, color='green', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='#E8F5E9', edgecolor='green', alpha=0.9))

    ax2.set_xlabel('Time (seconds)', fontsize=11)

    # Shared legend
    legend_elements = [
        mpatches.Patch(facecolor='#FFA726', alpha=0.5, label='Persistent (GIX + GDB)'),
        mpatches.Patch(facecolor='#42A5F5', alpha=0.5, label='Temp files (seed pairs, alignment)'),
    ]
    fig.legend(handles=legend_elements, loc='upper center', ncol=2, fontsize=10,
               bbox_to_anchor=(0.5, 1.01))

    fig.suptitle('Storage Timeline Comparison: Baseline vs Early GIX Deletion\n(EXAMPLE dataset ~86 Mbp per genome, T=32)',
                 fontsize=13, fontweight='bold', y=1.05)
    fig.tight_layout()

    outpath = os.path.join(DOCS_DIR, 'storage_timeline_early_gix_comparison.png')
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Saved: {outpath}")


if __name__ == '__main__':
    main()
