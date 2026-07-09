#!/usr/bin/env python3
"""Generate figure comparing hash table (cache-unfriendly) vs sorted merge (cache-friendly)."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patches as FancyArrowPatch
import numpy as np
import os

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Colors
NAVY = '#1B2A4A'
RED = '#E74C3C'
GREEN = '#27AE60'
BLUE = '#3498DB'
ORANGE = '#E67E22'
LIGHT_GRAY = '#ECF0F1'
DARK_GRAY = '#7F8C8D'
CACHE_HIT = '#2ECC71'
CACHE_MISS = '#E74C3C'


def draw_memory_block(ax, x, y, w, h, n_cells=16, label="", highlight_cells=None,
                      cell_color=LIGHT_GRAY, highlight_color=ORANGE):
    """Draw a memory block as a row of cells."""
    cell_w = w / n_cells
    for i in range(n_cells):
        color = highlight_color if highlight_cells and i in highlight_cells else cell_color
        rect = plt.Rectangle((x + i * cell_w, y), cell_w, h,
                              facecolor=color, edgecolor='#BDC3C7', linewidth=0.5)
        ax.add_patch(rect)
    if label:
        ax.text(x + w / 2, y + h + 0.15, label, ha='center', va='bottom',
                fontsize=10, fontweight='bold', color=NAVY)


def draw_arrow(ax, x1, y1, x2, y2, color=RED, style='->', lw=1.5):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw))


def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7.5))

    # ================================================================
    # LEFT PANEL: Hash Table (Cache-Unfriendly)
    # ================================================================
    ax1.set_xlim(-0.5, 11)
    ax1.set_ylim(-1, 10)
    ax1.set_aspect('equal')
    ax1.axis('off')
    ax1.set_title('Hash Table Lookup\n(e.g., minimap2)', fontsize=15,
                  fontweight='bold', color=RED, pad=15)

    # Draw genome A k-mers on the left
    ax1.text(0.5, 9.3, 'Genome A k-mers', fontsize=11, fontweight='bold', color=NAVY)
    kmer_labels = ['k-mer 1', 'k-mer 2', 'k-mer 3', 'k-mer 4', 'k-mer 5']
    kmer_y = [8.3, 7.3, 6.3, 5.3, 4.3]
    for i, (label, y) in enumerate(zip(kmer_labels, kmer_y)):
        rect = plt.Rectangle((0.2, y), 2.0, 0.7, facecolor='#D5E8D4',
                              edgecolor='#82B366', linewidth=1.5, zorder=2)
        ax1.add_patch(rect)
        ax1.text(1.2, y + 0.35, label, ha='center', va='center', fontsize=9, zorder=3)

    # Draw hash function arrow
    ax1.text(3.0, 8.8, 'hash()', fontsize=10, color=DARK_GRAY, style='italic',
             rotation=0, ha='center')

    # Draw hash table (genome B) - large block representing memory
    ax1.text(7.0, 9.3, 'Hash Table (Genome B)\n~6-8 GB in RAM',
             fontsize=11, fontweight='bold', color=NAVY, ha='center')
    # Draw as a grid of cells
    table_x, table_y = 4.5, 0.5
    table_w, table_h = 5.0, 8.0
    n_rows, n_cols = 16, 8

    cell_w = table_w / n_cols
    cell_h = table_h / n_rows
    for r in range(n_rows):
        for c in range(n_cols):
            rect = plt.Rectangle((table_x + c * cell_w, table_y + r * cell_h),
                                 cell_w, cell_h, facecolor=LIGHT_GRAY,
                                 edgecolor='#D5D8DC', linewidth=0.3)
            ax1.add_patch(rect)

    # Draw random access arrows from k-mers to scattered hash table cells
    target_cells = [(2, 14), (6, 3), (1, 9), (5, 11), (3, 6)]  # (col, row) - scattered
    arrow_colors = [RED, '#C0392B', '#E74C3C', '#CB4335', '#B03A2E']

    for i, ((col, row), y_start) in enumerate(zip(target_cells, kmer_y)):
        tx = table_x + col * cell_w + cell_w / 2
        ty = table_y + row * cell_h + cell_h / 2
        # Highlight the target cell
        rect = plt.Rectangle((table_x + col * cell_w, table_y + row * cell_h),
                              cell_w, cell_h, facecolor=CACHE_MISS,
                              edgecolor=RED, linewidth=1.5, zorder=2)
        ax1.add_patch(rect)
        # Draw arrow
        draw_arrow(ax1, 2.2, y_start + 0.35, tx, ty, color=arrow_colors[i], lw=1.3)

    # Cache miss annotations
    for i, (col, row) in enumerate(target_cells):
        tx = table_x + col * cell_w + cell_w / 2
        ty = table_y + row * cell_h + cell_h / 2
        ax1.text(tx, ty, '\u2717', ha='center', va='center',
                 fontsize=8, color='white', fontweight='bold', zorder=3)

    # Bottom annotation
    ax1.text(5.5, -0.5,
             'Each lookup jumps to a random location\n'
             '\u2717 Cache miss every time (~200 cycles)',
             ha='center', va='top', fontsize=11, color=RED,
             fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#FDEDEC', edgecolor=RED, alpha=0.9))

    # ================================================================
    # RIGHT PANEL: Sorted Merge (Cache-Friendly)
    # ================================================================
    ax2.set_xlim(-0.5, 11)
    ax2.set_ylim(-1, 10)
    ax2.set_aspect('equal')
    ax2.axis('off')
    ax2.set_title('Sorted Linear Merge\n(FastGA)', fontsize=15,
                  fontweight='bold', color=GREEN, pad=15)

    # Draw sorted table A (left column)
    ax2.text(1.3, 9.3, 'Genome A\n(sorted k-mers)', fontsize=10,
             fontweight='bold', color=NAVY, ha='center')
    n_entries = 12
    entry_h = 0.55
    start_y = 8.0
    for i in range(n_entries):
        y = start_y - i * entry_h
        # Color: green gradient for current scan position
        if i < 5:
            color = '#D5F5E3' if i < 4 else CACHE_HIT
            edge = '#27AE60' if i == 4 else '#82B366'
            lw = 2.0 if i == 4 else 0.8
        else:
            color = LIGHT_GRAY
            edge = '#D5D8DC'
            lw = 0.5
        rect = plt.Rectangle((0.2, y), 2.2, entry_h - 0.05,
                              facecolor=color, edgecolor=edge, linewidth=lw)
        ax2.add_patch(rect)
        ax2.text(1.3, y + entry_h/2 - 0.02, f'ACGT...{i:02d}', ha='center', va='center',
                 fontsize=7, fontfamily='monospace', color=NAVY)

    # Current position pointer for A
    ax2.annotate('', xy=(0.15, start_y - 4 * entry_h + entry_h/2),
                 xytext=(-0.3, start_y - 4 * entry_h + entry_h/2),
                 arrowprops=dict(arrowstyle='->', color=GREEN, lw=2.5))
    ax2.text(-0.35, start_y - 4 * entry_h + entry_h/2, '\u25B6',
             fontsize=14, color=GREEN, ha='right', va='center')

    # Draw sorted table B (right column)
    ax2.text(5.3, 9.3, 'Genome B\n(sorted k-mers)', fontsize=10,
             fontweight='bold', color=NAVY, ha='center')
    for i in range(n_entries):
        y = start_y - i * entry_h
        if i < 5:
            color = '#D6EAF8' if i < 4 else '#5DADE2'
            edge = '#2E86C1' if i == 4 else '#85C1E9'
            lw = 2.0 if i == 4 else 0.8
        else:
            color = LIGHT_GRAY
            edge = '#D5D8DC'
            lw = 0.5
        rect = plt.Rectangle((4.2, y), 2.2, entry_h - 0.05,
                              facecolor=color, edgecolor=edge, linewidth=lw)
        ax2.add_patch(rect)
        ax2.text(5.3, y + entry_h/2 - 0.02, f'ACGT...{i:02d}', ha='center', va='center',
                 fontsize=7, fontfamily='monospace', color=NAVY)

    # Current position pointer for B
    ax2.annotate('', xy=(4.15, start_y - 4 * entry_h + entry_h/2),
                 xytext=(3.7, start_y - 4 * entry_h + entry_h/2),
                 arrowprops=dict(arrowstyle='->', color=BLUE, lw=2.5))
    ax2.text(3.65, start_y - 4 * entry_h + entry_h/2, '\u25B6',
             fontsize=14, color=BLUE, ha='right', va='center')

    # Match arrow between current entries
    y_match = start_y - 4 * entry_h + entry_h/2
    draw_arrow(ax2, 2.45, y_match, 4.15, y_match, color=GREEN, lw=2.5)
    ax2.text(3.3, y_match + 0.25, 'Match!', ha='center', va='bottom',
             fontsize=10, color=GREEN, fontweight='bold')

    # Direction arrows showing sequential scan
    for table_x_pos in [1.3, 5.3]:
        ax2.annotate('', xy=(table_x_pos, start_y - 5.5 * entry_h),
                     xytext=(table_x_pos, start_y - 3 * entry_h),
                     arrowprops=dict(arrowstyle='->', color=DARK_GRAY, lw=1.5,
                                     linestyle='--'))
    ax2.text(1.3, start_y - 6 * entry_h, 'scan \u2193', ha='center',
             fontsize=9, color=DARK_GRAY, style='italic')
    ax2.text(5.3, start_y - 6 * entry_h, 'scan \u2193', ha='center',
             fontsize=9, color=DARK_GRAY, style='italic')

    # CPU cache diagram on the right
    cache_x = 7.8
    ax2.text(cache_x + 1.2, 9.3, 'CPU Cache', fontsize=11,
             fontweight='bold', color=NAVY, ha='center')

    # Cache levels
    cache_levels = [
        ('L1', 0.8, 7.8, CACHE_HIT, '~4 cycles'),
        ('L2', 1.2, 6.6, '#58D68D', '~12 cycles'),
        ('L3', 1.6, 5.2, '#82E0AA', '~40 cycles'),
    ]
    for label, w_scale, cy, color, latency in cache_levels:
        bw = 2.4 * w_scale
        rect = plt.Rectangle((cache_x + 1.2 - bw/2, cy), bw, 0.9,
                              facecolor=color, edgecolor='#1E8449',
                              linewidth=1.2, alpha=0.6, zorder=1)
        ax2.add_patch(rect)
        ax2.text(cache_x + 1.2, cy + 0.45, f'{label}\n{latency}',
                 ha='center', va='center', fontsize=8, fontweight='bold', zorder=2)

    # RAM at bottom
    rect = plt.Rectangle((cache_x - 0.5, 3.5), 3.4, 1.0,
                          facecolor=LIGHT_GRAY, edgecolor=DARK_GRAY,
                          linewidth=1.2, alpha=0.6)
    ax2.add_patch(rect)
    ax2.text(cache_x + 1.2, 4.0, 'RAM\n~200 cycles', ha='center', va='center',
             fontsize=8, color=DARK_GRAY)

    # Arrow showing prefetcher
    ax2.annotate('Prefetcher loads\nnext entries ahead',
                 xy=(cache_x + 1.2, 7.7), xytext=(cache_x + 1.2, 5.0),
                 arrowprops=dict(arrowstyle='->', color=GREEN, lw=2,
                                 connectionstyle='arc3,rad=0.3'),
                 fontsize=9, color=GREEN, ha='center', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='#EAFAF1',
                           edgecolor=GREEN, alpha=0.9))

    # Bottom annotation
    ax2.text(5.5, -0.5,
             'Both tables scanned left-to-right, one entry at a time\n'
             '\u2713 Sequential access \u2192 prefetcher keeps cache warm',
             ha='center', va='top', fontsize=11, color=GREEN,
             fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#EAFAF1',
                       edgecolor=GREEN, alpha=0.9))

    # ================================================================
    # Overall title
    # ================================================================
    fig.suptitle('Why Cache Coherence Matters: Hash Table vs Sorted Merge',
                 fontsize=17, fontweight='bold', color=NAVY, y=1.02)

    fig.tight_layout()
    outpath = os.path.join(OUTPUT_DIR, 'cache_comparison.png')
    fig.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")


if __name__ == '__main__':
    main()
