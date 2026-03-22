#!/usr/bin/env python3
"""Compare storage timeline: baseline vs Opt1 vs Opt1+Opt3 (T=32, EXAMPLE dataset)."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import csv
import os

BENCH_DIR = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/benchmarks"
AUDIT_DIR = os.path.join(BENCH_DIR, "storage_audit_v2")
DOCS_DIR = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/docs/storage_optimization"

# Phase timing for T=32 (from previous benchmarks)
PHASES = [
    (0, 0.5, "GDB1", "#E8E8E8"),
    (0.5, 4.3, "GIX1", "#FFE0B2"),
    (4.3, 4.8, "GDB2", "#E8E8E8"),
    (4.8, 8.7, "GIX2", "#FFE0B2"),
    (8.7, 14.4, "Seed Merge", "#BBDEFB"),
    (14.4, 21.5, "Sort + Align", "#C8E6C9"),
]

GDB_PER_GENOME = 21  # MB


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


def build_persistent(t_fine, gix_per_genome, keep_gix=True, delete_after=14.4):
    """Model persistent file sizes over time."""
    persistent = np.zeros_like(t_fine)
    for i, t in enumerate(t_fine):
        p = 0.0
        if t >= 0.5:
            p += GDB_PER_GENOME
            frac = min(1.0, max(0, (t - 0.5) / 3.8))
            if keep_gix or t < delete_after:
                p += gix_per_genome * frac
        if t >= 4.8:
            p += GDB_PER_GENOME
            frac = min(1.0, max(0, (t - 4.8) / 3.9))
            if keep_gix or t < delete_after:
                p += gix_per_genome * frac
        persistent[i] = p
    return persistent


def build_temp(t_fine, peak_mb=438, merge_start=8.7, merge_end=14.4, drain_end=21.5):
    """Model temp file usage."""
    temp = np.zeros_like(t_fine)
    for i, t in enumerate(t_fine):
        if merge_start <= t < merge_end:
            frac = (t - merge_start) / (merge_end - merge_start)
            temp[i] = peak_mb * frac
        elif merge_end <= t < drain_end:
            frac = 1.0 - (t - merge_end) / (drain_end - merge_end)
            temp[i] = peak_mb * max(0, frac)
    return temp


def plot_panel(ax, t_fine, persistent, temp, phases, title, ylim, annotate_peak=True, annotate_delete=False):
    for start, end, label, color in phases:
        ax.axvspan(start, end, alpha=0.15, color=color, zorder=0)
        ax.axvline(start, color='gray', linewidth=0.5, linestyle=':', alpha=0.5)
        mid = (start + end) / 2
        if end - start > 2:
            ax.text(mid, ylim * 0.93, label, ha='center', va='center',
                    fontsize=7, color='#555', style='italic')

    ax.fill_between(t_fine, persistent, alpha=0.5, color='#FFA726', zorder=1)
    ax.fill_between(t_fine, persistent, persistent + temp, alpha=0.5, color='#42A5F5', zorder=2)
    ax.plot(t_fine, persistent + temp, color='#1565C0', linewidth=0.8, zorder=3)
    ax.plot(t_fine, persistent, color='#E65100', linewidth=0.8, zorder=3, linestyle='--')

    total = persistent + temp
    peak_idx = np.argmax(total)
    peak_val = total[peak_idx]
    peak_time = t_fine[peak_idx]

    if annotate_peak:
        ax.annotate(f'Peak: {peak_val:.0f} MB',
                    xy=(peak_time, peak_val),
                    xytext=(peak_time + 2, min(peak_val + 200, ylim * 0.85)),
                    arrowprops=dict(arrowstyle='->', color='red', lw=1.2),
                    fontsize=9, color='red', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='red', alpha=0.9))

    if annotate_delete:
        # Find where persistent drops
        for i in range(1, len(t_fine)):
            if persistent[i] < persistent[i-1] - 100:
                drop_val = persistent[i-1] - persistent[i]
                ax.annotate(f'GIX deleted\n(-{drop_val:.0f} MB)',
                            xy=(t_fine[i], persistent[i] + temp[i]),
                            xytext=(t_fine[i] + 2.5, persistent[i] + temp[i] + 500),
                            arrowprops=dict(arrowstyle='->', color='green', lw=1.5),
                            fontsize=8, color='green', fontweight='bold',
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='#E8F5E9', edgecolor='green', alpha=0.9))
                break

    ax.set_ylabel('Disk Usage (MB)', fontsize=9)
    ax.set_title(title, fontsize=10, fontweight='bold', loc='left')
    ax.set_xlim(0, 23)
    ax.set_ylim(0, ylim)
    ax.tick_params(labelsize=8)
    ax.grid(axis='y', alpha=0.3)


def main():
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    ylim = 2800
    t_fine = np.linspace(0, 22, 2000)

    GIX_BASELINE = 930  # MB per genome (with mask byte)
    GIX_OPT3 = 868      # MB per genome (without mask byte)

    # --- Panel 1: Baseline (GIX kept, with mask byte) ---
    persistent_base = build_persistent(t_fine, GIX_BASELINE, keep_gix=True)
    temp_base = build_temp(t_fine)
    plot_panel(ax1, t_fine, persistent_base, temp_base, PHASES,
               "Baseline (original FastGA) — GIX with mask byte, kept throughout", ylim)

    # --- Panel 2: Opt1 only (early GIX deletion, with mask byte) ---
    persistent_opt1 = build_persistent(t_fine, GIX_BASELINE, keep_gix=False)
    temp_opt1 = build_temp(t_fine)
    plot_panel(ax2, t_fine, persistent_opt1, temp_opt1, PHASES,
               "Opt 1 (Early GIX Deletion) — GIX with mask byte, freed after seed merge", ylim,
               annotate_delete=True)

    # --- Panel 3: Opt1 + Opt3 (early deletion + no mask byte) ---
    persistent_opt13 = build_persistent(t_fine, GIX_OPT3, keep_gix=False)
    temp_opt13 = build_temp(t_fine)
    plot_panel(ax3, t_fine, persistent_opt13, temp_opt13, PHASES,
               "Opt 1 + Opt 3 (Early Deletion + No Mask Byte) — smaller GIX, freed after merge", ylim,
               annotate_delete=True)

    ax3.set_xlabel('Time (seconds)', fontsize=11)

    # Shared legend
    legend_elements = [
        mpatches.Patch(facecolor='#FFA726', alpha=0.5, label='Persistent (GIX + GDB)'),
        mpatches.Patch(facecolor='#42A5F5', alpha=0.5, label='Temp files (seed pairs, alignment)'),
    ]
    fig.legend(handles=legend_elements, loc='upper center', ncol=2, fontsize=10,
               bbox_to_anchor=(0.5, 1.01))

    fig.suptitle('Storage Timeline: Baseline → Opt 1 → Opt 1+3\n(EXAMPLE dataset ~86 Mbp per genome, T=32)',
                 fontsize=13, fontweight='bold', y=1.05)
    fig.tight_layout()

    outpath = os.path.join(DOCS_DIR, 'storage_timeline_opt3_comparison.png')
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Saved: {outpath}")


if __name__ == '__main__':
    main()
