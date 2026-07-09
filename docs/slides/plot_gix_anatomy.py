#!/usr/bin/env python3
"""Generate a clean figure showing the anatomy of a single GIX ktab entry."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

NAVY = '#1B2A4A'
DARK_GRAY = '#7F8C8D'


def main():
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 8)
    ax.axis('off')

    # Title
    ax.text(8, 7.3, 'Anatomy of a .ktab Entry (15 bytes per k-mer)',
            ha='center', fontsize=17, fontweight='bold', color=NAVY)
    ax.text(8, 6.6, '2.3 billion entries per human genome  \u2192  ~32 GB per genome  \u2192  ~63 GB for a pair',
            ha='center', fontsize=12, color=DARK_GRAY)

    # Entry bar
    entry_y = 3.8
    entry_h = 2.0
    fields = [
        ('K-mer Suffix',    7, '#3498DB', '7 bytes\n28 bases (2-bit encoded)'),
        ('Mask',            1, '#E74C3C', '1 byte'),
        ('LCP',             1, '#9B59B6', '1 byte'),
        ('Position',        4, '#E67E22', '4 bytes\noffset in contig'),
        ('Contig+Strand',   2, '#27AE60', '2 bytes\nwhich contig, which strand'),
    ]

    total_bytes = sum(f[1] for f in fields)
    bar_w = 13.0
    scale = bar_w / total_bytes
    x = 1.5

    for name, nbytes, color, detail in fields:
        w = nbytes * scale
        rect = plt.Rectangle((x, entry_y), w, entry_h,
                              facecolor=color, edgecolor='white',
                              linewidth=2.5, alpha=0.9)
        ax.add_patch(rect)
        # Field name inside bar
        fontsize = 13 if nbytes >= 4 else 11 if nbytes >= 2 else 9
        ax.text(x + w/2, entry_y + entry_h/2 + 0.15, name,
                ha='center', va='center', fontsize=fontsize,
                fontweight='bold', color='white')
        # Byte count inside bar
        ax.text(x + w/2, entry_y + entry_h/2 - 0.35, f'{nbytes}B',
                ha='center', va='center', fontsize=10, color='white', alpha=0.9)
        # Detail below bar
        ax.text(x + w/2, entry_y - 0.4, detail,
                ha='center', va='top', fontsize=9, color=color, fontweight='bold')
        x += w

    # Bracket above showing total
    ax.annotate('', xy=(1.5, entry_y + entry_h + 0.3),
                xytext=(1.5 + bar_w, entry_y + entry_h + 0.3),
                arrowprops=dict(arrowstyle='|-|', color=NAVY, lw=1.5))
    ax.text(8, entry_y + entry_h + 0.55, '15 bytes total',
            ha='center', fontsize=12, color=NAVY, fontweight='bold')

    # Percentage labels at bottom
    pct_y = 1.2
    x = 1.5
    for name, nbytes, color, _ in fields:
        w = nbytes * scale
        pct = nbytes / total_bytes * 100
        ax.text(x + w/2, pct_y, f'{pct:.0f}%',
                ha='center', va='center', fontsize=11, color=color, fontweight='bold')
        x += w

    fig.tight_layout()
    outpath = os.path.join(OUTPUT_DIR, 'gix_anatomy.png')
    fig.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")


if __name__ == '__main__':
    main()
