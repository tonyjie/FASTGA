#!/usr/bin/env python3
"""Plot storage usage timeline for FastGA benchmark runs."""

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

# Phase boundaries from the T08 verbose log (approximate wall clock offsets):
# GDB HAP1: 0 - 0.5s
# GIX HAP1: 0.5 - 3.6s
# GDB HAP2: 3.6 - 4.1s
# GIX HAP2: 4.1 - 7.2s
# Seed Merge: 7.2 - 15.5s
# Sort + Align: 15.5 - 32.1s
# PAF Conv: 32.1 - 32.3s

PHASES_T08 = [
    (0, 0.5, "GDB\nHAP1", "#E8E8E8"),
    (0.5, 3.6, "GIX Build\nHAP1", "#FFE0B2"),
    (3.6, 4.1, "GDB\nHAP2", "#E8E8E8"),
    (4.1, 7.2, "GIX Build\nHAP2", "#FFE0B2"),
    (7.2, 15.5, "Seed Merge", "#BBDEFB"),
    (15.5, 32.1, "Sort + Align", "#C8E6C9"),
    (32.1, 32.5, "PAF", "#F3E5F5"),
]


def read_monitor(filepath):
    """Read monitor TSV and return (time, mb) arrays."""
    times, mbs = [], []
    with open(filepath) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip header
        for row in reader:
            if len(row) >= 2:
                t = float(row[0])
                b = float(row[1])
                times.append(t)
                mbs.append(b / 1048576)  # bytes to MB
    return np.array(times), np.array(mbs)


def plot_single_thread(thread_count, ax, phases, show_xlabel=True):
    """Plot storage timeline for a single thread count."""
    tdir = f"T{thread_count:02d}"
    monfile = os.path.join(AUDIT_DIR, f"{tdir}_tmpdir_monitor.tsv")
    if not os.path.exists(monfile):
        return

    times, mbs = read_monitor(monfile)

    # Draw phase backgrounds
    for start, end, label, color in phases:
        ax.axvspan(start, end, alpha=0.3, color=color, zorder=0)

    # Plot temp file usage
    ax.fill_between(times, mbs, alpha=0.6, color='#1976D2', zorder=2)
    ax.plot(times, mbs, color='#0D47A1', linewidth=0.8, zorder=3)

    # Add persistent file line (GIX appears during GIX build phase)
    # GIX HAP1 appears ~0.5s, HAP2 ~4.1s, full GIX by ~7.2s
    persistent = np.zeros_like(times)
    gix_per_genome = 930  # MB
    gdb_per_genome = 21   # MB
    for i, t in enumerate(times):
        if t < 0.5:
            persistent[i] = 0
        elif t < 3.6:
            # GDB1 done, GIX1 building
            persistent[i] = gdb_per_genome + min(gix_per_genome, gix_per_genome * (t - 0.5) / 3.1)
        elif t < 4.1:
            persistent[i] = gdb_per_genome + gix_per_genome
        elif t < 7.2:
            # GDB2 done, GIX2 building
            persistent[i] = 2 * gdb_per_genome + gix_per_genome + min(gix_per_genome, gix_per_genome * (t - 4.1) / 3.1)
        else:
            persistent[i] = 2 * (gdb_per_genome + gix_per_genome)

    ax.fill_between(times, persistent + mbs, persistent, alpha=0.4, color='#FF7043', zorder=1, label='Temp files (tmpdir)')
    ax.fill_between(times, persistent, alpha=0.4, color='#FFA726', zorder=1, label='Persistent (GIX+GDB)')

    # Redraw temp on top of persistent for stacked view
    ax.fill_between(times, 0, 0, alpha=0, zorder=0)  # dummy

    # Actually let's do a clean stacked area
    ax.clear()

    # Phase backgrounds
    for start, end, label, color in phases:
        ax.axvspan(start, end, alpha=0.15, color=color, zorder=0)
        mid = (start + end) / 2
        ax.text(mid, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 2400, label,
                ha='center', va='top', fontsize=6, color='#555555', zorder=5)

    # Stacked: persistent on bottom, temp on top
    ax.fill_between(times, persistent, alpha=0.5, color='#FFA726', zorder=1)
    ax.fill_between(times, persistent, persistent + mbs, alpha=0.5, color='#42A5F5', zorder=2)
    ax.plot(times, persistent + mbs, color='#1565C0', linewidth=0.8, zorder=3)
    ax.plot(times, persistent, color='#E65100', linewidth=0.8, zorder=3, linestyle='--')

    # Peak annotation
    total = persistent + mbs
    peak_idx = np.argmax(total)
    peak_val = total[peak_idx]
    peak_time = times[peak_idx]
    ax.annotate(f'Peak: {peak_val:.0f} MB\n({mbs[peak_idx]:.0f} MB temp)',
                xy=(peak_time, peak_val), xytext=(peak_time + 3, peak_val + 100),
                arrowprops=dict(arrowstyle='->', color='red', lw=1.2),
                fontsize=7, color='red', fontweight='bold', zorder=6)

    ax.set_ylabel('Disk Usage (MB)', fontsize=9)
    if show_xlabel:
        ax.set_xlabel('Time (seconds)', fontsize=9)
    ax.set_title(f'T={thread_count}', fontsize=10, fontweight='bold', loc='left')
    ax.set_xlim(0, max(times) * 1.05)
    ax.set_ylim(0, 2800)
    ax.tick_params(labelsize=8)
    ax.grid(axis='y', alpha=0.3)


def main():
    # --- Figure 1: Detailed T=8 timeline ---
    fig1, ax1 = plt.subplots(1, 1, figsize=(12, 5))

    tdir = "T08"
    monfile = os.path.join(AUDIT_DIR, f"{tdir}_tmpdir_monitor.tsv")
    times, mbs = read_monitor(monfile)

    # Persistent file model
    gix_per_genome = 930
    gdb_per_genome = 21
    persistent = np.zeros_like(times)
    for i, t in enumerate(times):
        if t < 0.5:
            persistent[i] = 0
        elif t < 3.6:
            frac = min(1.0, (t - 0.5) / 3.1)
            persistent[i] = gdb_per_genome + gix_per_genome * frac
        elif t < 4.1:
            persistent[i] = gdb_per_genome + gix_per_genome
        elif t < 7.2:
            frac = min(1.0, (t - 4.1) / 3.1)
            persistent[i] = 2 * gdb_per_genome + gix_per_genome + gix_per_genome * frac
        else:
            persistent[i] = 2 * (gdb_per_genome + gix_per_genome)

    # Phase backgrounds with labels
    for start, end, label, color in PHASES_T08:
        ax1.axvspan(start, end, alpha=0.15, color=color, zorder=0)
        mid = (start + end) / 2
        if end - start > 1.5:  # only label wide phases
            ax1.text(mid, 2650, label.replace('\n', ' '), ha='center', va='center',
                     fontsize=8, color='#444444', style='italic', zorder=5)

    # Stacked areas
    ax1.fill_between(times, persistent, alpha=0.5, color='#FFA726', zorder=1, label='Persistent files (GIX + GDB)')
    ax1.fill_between(times, persistent, persistent + mbs, alpha=0.5, color='#42A5F5', zorder=2, label='Temp files (seed pairs, alignment)')
    ax1.plot(times, persistent + mbs, color='#1565C0', linewidth=1.0, zorder=3)
    ax1.plot(times, persistent, color='#E65100', linewidth=1.0, zorder=3, linestyle='--')

    # Peak annotation
    total = persistent + mbs
    peak_idx = np.argmax(total)
    peak_val = total[peak_idx]
    peak_time = times[peak_idx]
    ax1.annotate(f'Peak total: {peak_val:.0f} MB\n(Persistent: {persistent[peak_idx]:.0f} MB + Temp: {mbs[peak_idx]:.0f} MB)',
                 xy=(peak_time, peak_val), xytext=(peak_time + 4, peak_val + 200),
                 arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
                 fontsize=9, color='red', fontweight='bold', zorder=6,
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='red', alpha=0.9))

    # Phase boundary lines
    for start, end, label, color in PHASES_T08:
        ax1.axvline(start, color='gray', linewidth=0.5, linestyle=':', alpha=0.5, zorder=0)

    ax1.set_xlabel('Time (seconds)', fontsize=11)
    ax1.set_ylabel('Disk Usage (MB)', fontsize=11)
    ax1.set_title('FastGA Storage Usage Timeline (T=8, EXAMPLE dataset ~86 Mbp per genome)',
                  fontsize=12, fontweight='bold')
    ax1.set_xlim(0, max(times) * 1.02)
    ax1.set_ylim(0, 2900)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.grid(axis='y', alpha=0.3)

    fig1.tight_layout()
    fig1.savefig(os.path.join(DOCS_DIR, 'storage_timeline_T08.png'), dpi=150, bbox_inches='tight')
    print(f"Saved: {os.path.join(DOCS_DIR, 'storage_timeline_T08.png')}")

    # --- Figure 2: Comparison across thread counts ---
    thread_counts = [1, 4, 8, 32]
    fig2, axes = plt.subplots(len(thread_counts), 1, figsize=(12, 3 * len(thread_counts)), sharex=False)

    # Phase timing per thread count (approximate from logs)
    phase_timings = {
        1:  [(0, 0.5), (0.5, 9.8), (9.8, 10.3), (10.3, 19.9), (19.9, 49.9), (49.9, 133.0), (133.0, 133.5)],
        4:  [(0, 0.5), (0.5, 4.4), (4.4, 4.9), (4.9, 8.5), (8.5, 19.7), (19.7, 48.3), (48.3, 48.5)],
        8:  [(0, 0.5), (0.5, 3.6), (3.6, 4.1), (4.1, 7.2), (7.2, 15.5), (15.5, 32.1), (32.1, 32.3)],
        32: [(0, 0.5), (0.5, 4.3), (4.3, 4.8), (4.8, 8.7), (8.7, 14.4), (14.4, 21.5), (21.5, 21.7)],
    }
    phase_colors = ["#E8E8E8", "#FFE0B2", "#E8E8E8", "#FFE0B2", "#BBDEFB", "#C8E6C9", "#F3E5F5"]
    phase_labels = ["GDB1", "GIX1", "GDB2", "GIX2", "Seed Merge", "Sort+Align", "PAF"]

    for idx, (tc, ax) in enumerate(zip(thread_counts, axes)):
        tdir = f"T{tc:02d}"
        monfile = os.path.join(AUDIT_DIR, f"{tdir}_tmpdir_monitor.tsv")
        if not os.path.exists(monfile):
            continue

        times, mbs = read_monitor(monfile)
        timings = phase_timings[tc]

        # Persistent model
        persistent = np.zeros_like(times)
        gix_start1, gix_end1 = timings[1]
        gix_start2, gix_end2 = timings[3]
        for i, t in enumerate(times):
            if t < gix_start1:
                persistent[i] = 0
            elif t < gix_end1:
                frac = min(1.0, (t - gix_start1) / (gix_end1 - gix_start1))
                persistent[i] = gdb_per_genome + gix_per_genome * frac
            elif t < gix_start2:
                persistent[i] = gdb_per_genome + gix_per_genome
            elif t < gix_end2:
                frac = min(1.0, (t - gix_start2) / (gix_end2 - gix_start2))
                persistent[i] = 2 * gdb_per_genome + gix_per_genome + gix_per_genome * frac
            else:
                persistent[i] = 2 * (gdb_per_genome + gix_per_genome)

        # Phase backgrounds
        for (start, end), color, label in zip(timings, phase_colors, phase_labels):
            ax.axvspan(start, end, alpha=0.15, color=color, zorder=0)
            if end - start > (max(times) * 0.05):
                ax.text((start + end) / 2, 2650, label, ha='center', va='center',
                        fontsize=7, color='#555', style='italic')

        # Stacked areas
        ax.fill_between(times, persistent, alpha=0.5, color='#FFA726', zorder=1)
        ax.fill_between(times, persistent, persistent + mbs, alpha=0.5, color='#42A5F5', zorder=2)
        ax.plot(times, persistent + mbs, color='#1565C0', linewidth=0.8, zorder=3)
        ax.plot(times, persistent, color='#E65100', linewidth=0.8, zorder=3, linestyle='--')

        total = persistent + mbs
        peak_idx = np.argmax(total)
        ax.annotate(f'Peak: {total[peak_idx]:.0f} MB',
                    xy=(times[peak_idx], total[peak_idx]),
                    xytext=(times[peak_idx] + max(times) * 0.08, total[peak_idx] + 150),
                    arrowprops=dict(arrowstyle='->', color='red', lw=1),
                    fontsize=8, color='red', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='red', alpha=0.8))

        ax.set_ylabel('MB', fontsize=9)
        ax.set_title(f'T={tc} (wall: {max(times):.0f}s)', fontsize=10, fontweight='bold', loc='left')
        ax.set_ylim(0, 2900)
        ax.set_xlim(0, max(times) * 1.05)
        ax.tick_params(labelsize=8)
        ax.grid(axis='y', alpha=0.3)

    axes[-1].set_xlabel('Time (seconds)', fontsize=11)

    # Shared legend
    legend_elements = [
        mpatches.Patch(facecolor='#FFA726', alpha=0.5, label='Persistent (GIX + GDB)'),
        mpatches.Patch(facecolor='#42A5F5', alpha=0.5, label='Temp files (seed pairs, alignment)'),
    ]
    fig2.legend(handles=legend_elements, loc='upper center', ncol=2, fontsize=10,
                bbox_to_anchor=(0.5, 1.02))

    fig2.suptitle('Storage Usage Timeline Across Thread Counts\n(EXAMPLE dataset ~86 Mbp per genome)',
                  fontsize=13, fontweight='bold', y=1.06)
    fig2.tight_layout()
    fig2.savefig(os.path.join(DOCS_DIR, 'storage_timeline_comparison.png'), dpi=150, bbox_inches='tight')
    print(f"Saved: {os.path.join(DOCS_DIR, 'storage_timeline_comparison.png')}")


if __name__ == '__main__':
    main()
