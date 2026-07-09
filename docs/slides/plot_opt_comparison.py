#!/usr/bin/env python3
"""Compare storage timeline: Opt 1 vs Opt 1 + Opt 2 (EXAMPLE dataset, T=32).
   Note: "Opt 2" in this figure corresponds to "Opt 3 (Eliminate Mask Byte)"
   in the codebase, renumbered for presentation clarity."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Phase timing for T=32 (from benchmarks)
PHASES = [
    (0, 0.5, "GDB1", "#E8E8E8"),
    (0.5, 4.3, "GIX1", "#FFE0B2"),
    (4.3, 4.8, "GDB2", "#E8E8E8"),
    (4.8, 8.7, "GIX2", "#FFE0B2"),
    (8.7, 14.4, "Seed Merge", "#BBDEFB"),
    (14.4, 21.5, "Sort + Align", "#C8E6C9"),
]

GDB_PER_GENOME = 21  # MB


def build_persistent(t_fine, gix_per_genome, delete_after=14.4):
    """Model persistent file sizes — GIX always deleted after seed merge."""
    persistent = np.zeros_like(t_fine)
    for i, t in enumerate(t_fine):
        p = 0.0
        if t >= 0.5:
            p += GDB_PER_GENOME
            frac = min(1.0, max(0, (t - 0.5) / 3.8))
            if t < delete_after:
                p += gix_per_genome * frac
        if t >= 4.8:
            p += GDB_PER_GENOME
            frac = min(1.0, max(0, (t - 4.8) / 3.9))
            if t < delete_after:
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


def plot_panel(ax, t_fine, persistent, temp, phases, title, ylim,
               annotate_peak=True, annotate_delete=True):
    # Phase background bands
    for start, end, label, color in phases:
        ax.axvspan(start, end, alpha=0.15, color=color, zorder=0)
        ax.axvline(start, color='gray', linewidth=0.5, linestyle=':', alpha=0.5)
        mid = (start + end) / 2
        if end - start > 2:
            ax.text(mid, ylim * 0.93, label, ha='center', va='center',
                    fontsize=8, color='#555', style='italic')

    # Fill areas
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
                    fontsize=10, color='red', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                              edgecolor='red', alpha=0.9))

    if annotate_delete:
        for i in range(1, len(t_fine)):
            if persistent[i] < persistent[i-1] - 100:
                drop_val = persistent[i-1] - persistent[i]
                ax.annotate(f'GIX deleted\n(-{drop_val:.0f} MB)',
                            xy=(t_fine[i], persistent[i] + temp[i]),
                            xytext=(t_fine[i] + 2.5, persistent[i] + temp[i] + 500),
                            arrowprops=dict(arrowstyle='->', color='green', lw=1.5),
                            fontsize=9, color='green', fontweight='bold',
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='#E8F5E9',
                                      edgecolor='green', alpha=0.9))
                break

    # Sort+align disk annotation
    sort_align_start = 14.4
    for i, t in enumerate(t_fine):
        if t >= sort_align_start + 0.5:
            sa_disk = persistent[i] + temp[i]
            ax.annotate(f'Sort+align: {sa_disk:.0f} MB',
                        xy=(t, sa_disk),
                        xytext=(t + 1, sa_disk + 300),
                        arrowprops=dict(arrowstyle='->', color='#1565C0', lw=1),
                        fontsize=9, color='#1565C0',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                  edgecolor='#1565C0', alpha=0.8))
            break

    ax.set_ylabel('Disk Usage (MB)', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold', loc='left')
    ax.set_xlim(0, 23)
    ax.set_ylim(0, ylim)
    ax.tick_params(labelsize=9)
    ax.grid(axis='y', alpha=0.3)


def main():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True)
    ylim = 2800
    t_fine = np.linspace(0, 22, 2000)

    GIX_WITH_MASK = 930     # MB per genome (original, with mask byte)
    GIX_WITHOUT_MASK = 868  # MB per genome (without mask byte)

    # --- Panel 1: Opt 1 (early GIX deletion, with mask byte) ---
    persistent_opt1 = build_persistent(t_fine, GIX_WITH_MASK)
    temp_opt1 = build_temp(t_fine)
    plot_panel(ax1, t_fine, persistent_opt1, temp_opt1, PHASES,
               "Opt 1 (Early GIX Deletion) \u2014 GIX freed after seed merge", ylim)

    # --- Panel 2: Opt 1 + Opt 2 (early deletion + no mask byte) ---
    persistent_opt12 = build_persistent(t_fine, GIX_WITHOUT_MASK)
    temp_opt12 = build_temp(t_fine)
    plot_panel(ax2, t_fine, persistent_opt12, temp_opt12, PHASES,
               "Opt 1 + Opt 2 (Early Deletion + Eliminate Mask Byte) \u2014 smaller GIX, freed after merge",
               ylim)

    ax2.set_xlabel('Time (seconds)', fontsize=11)

    # Shared legend
    legend_elements = [
        mpatches.Patch(facecolor='#FFA726', alpha=0.5, label='Persistent (GIX + GDB)'),
        mpatches.Patch(facecolor='#42A5F5', alpha=0.5, label='Temp files (seed pairs, alignment)'),
    ]
    fig.legend(handles=legend_elements, loc='upper center', ncol=2, fontsize=11,
               bbox_to_anchor=(0.5, 1.01))

    fig.suptitle('Storage Optimization Comparison: Opt 1 vs Opt 1 + Opt 2\n'
                 '(EXAMPLE dataset ~86 Mbp per genome, T=32)',
                 fontsize=14, fontweight='bold', y=1.07)
    fig.tight_layout()

    outpath = os.path.join(OUTPUT_DIR, 'opt_comparison.png')
    fig.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")


if __name__ == '__main__':
    main()
