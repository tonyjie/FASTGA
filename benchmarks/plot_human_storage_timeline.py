#!/usr/bin/env python3
"""Plot storage usage timeline for FastGA human genome benchmark (GRCh38 vs CHM13, T=32).

Since the tmpdir monitor didn't capture temp files for this run, we reconstruct
the timeline from:
1. Per-phase wall clock times from verbose log → when GIX files are built
2. Disk free space polling (30s interval) → actual disk consumption over time
3. Known persistent file sizes → GIX/GDB for each genome
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

DOCS_DIR = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/docs"

# === Phase timing from verbose log (wall clock seconds) ===
# GDB GRCh38:     0.0 - 9.2
# GIX GRCh38:     9.2 - 48.3   (39.1s)
# GDB CHM13:     48.3 - 57.4
# GIX CHM13:     57.4 - 100.1  (42.7s)
# Seed merge:   100.1 - 105.4  (5.3s)
# Sort+align:   105.4 - 587.3  (481.9s)
# (PAF not run — output was .1aln)

PHASES = [
    (0, 9.2, "GDB\nGRCh38", "#E8E8E8"),
    (9.2, 48.3, "GIX Build\nGRCh38", "#FFE0B2"),
    (48.3, 57.4, "GDB\nCHM13", "#E8E8E8"),
    (57.4, 100.1, "GIX Build\nCHM13", "#FFE0B2"),
    (100.1, 105.4, "Seed\nMerge", "#BBDEFB"),
    (105.4, 587.3, "Sort + Align (66 parts)", "#C8E6C9"),
]

# === Persistent file sizes (measured) ===
GIX_GRCH38 = 32.49  # GB
GIX_CHM13 = 30.18   # GB
GDB_GRCH38 = 0.748  # GB
GDB_CHM13 = 0.743   # GB

# === Disk free space polling data (from benchmark output) ===
# timestamp_offset(s), free_GB
disk_free_data = [
    (0, 308),
    (30, 290),
    (60, 272),
    (90, 250),
    (120, 237),
    (150, 239),
    (180, 239),
    (210, 239),
    (240, 239),
    (270, 240),
    (300, 241),
    (330, 242),
    (360, 243),
    (390, 243),
    (420, 243),
    (450, 243),
    (480, 243),
    (510, 243),
    (540, 243),
    (570, 243),
    (588, 244),  # end
]

BASELINE_FREE = 308  # GB free at start


def main():
    fig, ax = plt.subplots(1, 1, figsize=(14, 6))

    # --- Build persistent storage model (smooth) ---
    # GIX builds linearly during its phase
    t_fine = np.linspace(0, 590, 2000)
    persistent = np.zeros_like(t_fine)
    temp = np.zeros_like(t_fine)

    for i, t in enumerate(t_fine):
        p = 0.0
        if t >= 0:
            # GDB GRCh38 created at ~9.2s (instant, small)
            if t >= 9.2:
                p += GDB_GRCH38
        if t >= 9.2:
            # GIX GRCh38 builds from 9.2 to 48.3
            frac = min(1.0, max(0, (t - 9.2) / (48.3 - 9.2)))
            p += GIX_GRCH38 * frac
        if t >= 57.4:
            # GDB CHM13
            p += GDB_CHM13
            # GIX CHM13 builds from 57.4 to 100.1
            frac = min(1.0, max(0, (t - 57.4) / (100.1 - 57.4)))
            p += GIX_CHM13 * frac
        elif t >= 48.3:
            p += GDB_CHM13 if t >= 48.3 else 0

        persistent[i] = p

    # --- Compute temp storage from disk free data ---
    # Total consumed = BASELINE_FREE - free_at_time
    # Temp = total_consumed - persistent_at_time
    disk_times = np.array([d[0] for d in disk_free_data])
    disk_free = np.array([d[1] for d in disk_free_data])
    disk_consumed = BASELINE_FREE - disk_free

    # Interpolate persistent at disk poll times
    persistent_at_polls = np.interp(disk_times, t_fine, persistent)
    temp_at_polls = np.maximum(0, disk_consumed - persistent_at_polls)

    # Interpolate temp to fine time grid
    temp = np.interp(t_fine, disk_times, temp_at_polls)
    # Smooth out negative values
    temp = np.maximum(0, temp)

    # --- Phase backgrounds ---
    for start, end, label, color in PHASES:
        ax.axvspan(start, end, alpha=0.15, color=color, zorder=0)
        ax.axvline(start, color='gray', linewidth=0.5, linestyle=':', alpha=0.5, zorder=0)
        mid = (start + end) / 2
        if end - start > 20:
            ax.text(mid, 74, label.replace('\n', ' '), ha='center', va='center',
                    fontsize=8, color='#444444', style='italic', zorder=5)

    # --- Stacked areas ---
    ax.fill_between(t_fine, persistent, alpha=0.5, color='#FFA726', zorder=1,
                    label='Persistent files (GIX + GDB)')
    ax.fill_between(t_fine, persistent, persistent + temp, alpha=0.5, color='#42A5F5', zorder=2,
                    label='Temp files (seed pairs, alignment)')
    ax.plot(t_fine, persistent + temp, color='#1565C0', linewidth=1.0, zorder=3)
    ax.plot(t_fine, persistent, color='#E65100', linewidth=1.0, zorder=3, linestyle='--')

    # --- Peak annotation ---
    total = persistent + temp
    peak_idx = np.argmax(total)
    peak_val = total[peak_idx]
    peak_time = t_fine[peak_idx]
    peak_persistent = persistent[peak_idx]
    peak_temp = temp[peak_idx]

    ax.annotate(
        f'Peak total: {peak_val:.0f} GB\n(Persistent: {peak_persistent:.0f} GB + Temp: {peak_temp:.0f} GB)',
        xy=(peak_time, peak_val), xytext=(peak_time + 80, peak_val + 5),
        arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
        fontsize=9, color='red', fontweight='bold', zorder=6,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='red', alpha=0.9))

    # --- Disk free scatter points ---
    ax.scatter(disk_times, disk_consumed, color='#D32F2F', s=15, zorder=7,
               label='Measured disk consumed (30s polling)', edgecolors='white', linewidth=0.5)

    # --- Labels ---
    ax.set_xlabel('Time (seconds)', fontsize=11)
    ax.set_ylabel('Disk Usage (GB)', fontsize=11)
    ax.set_title('FastGA Storage Usage Timeline (T=32, Human Genome: GRCh38 vs CHM13, ~3.1 Gbp each)',
                 fontsize=12, fontweight='bold')
    ax.set_xlim(0, 600)
    ax.set_ylim(0, 80)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    # Add time axis in minutes too
    ax2 = ax.twiny()
    ax2.set_xlim(0, 600/60)
    ax2.set_xlabel('Time (minutes)', fontsize=10)
    ax2.tick_params(labelsize=8)

    fig.tight_layout()
    outpath = f"{DOCS_DIR}/storage_timeline_human_T32.png"
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Saved: {outpath}")


if __name__ == '__main__':
    main()
