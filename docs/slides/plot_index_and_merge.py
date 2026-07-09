#!/usr/bin/env python3
"""Generate figures explaining:
1. How FastGA builds the sorted index (MSD radix sort)
2. How the linear merge finds seeds (adaptamer matching)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

NAVY = '#1B2A4A'
RED = '#E74C3C'
GREEN = '#27AE60'
BLUE = '#3498DB'
ORANGE = '#E67E22'
PURPLE = '#8E44AD'
LIGHT_GRAY = '#ECF0F1'
DARK_GRAY = '#7F8C8D'
YELLOW = '#F1C40F'


# ================================================================
# FIGURE 1: Building the Sorted Index (MSD Radix Sort)
# ================================================================

def plot_index_building():
    fig, axes = plt.subplots(1, 4, figsize=(18, 8),
                              gridspec_kw={'width_ratios': [2.5, 3, 3, 2.5]})

    # ---- Panel 1: Genome scan → k-mers ----
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis('off')
    ax.set_title('Step 1: Extract K-mers\nfrom Genome', fontsize=13, fontweight='bold',
                 color=NAVY, pad=10)

    # Draw genome as a long bar
    genome_y = 10.5
    rect = plt.Rectangle((0.5, genome_y), 9, 0.6, facecolor='#D5E8D4',
                          edgecolor='#82B366', linewidth=1.5)
    ax.add_patch(rect)
    ax.text(5, genome_y + 0.3, 'Genome: ...ACGTTAGCACGTCCGA...', ha='center', va='center',
            fontsize=8, fontfamily='monospace')

    # Show sliding window extracting k-mers
    kmers = [
        ('ACGTTAGCAC...', True),
        ('CGTTAGCACG...', False),
        ('GTTAGCACGT...', True),
        ('TTAGCACGTC...', False),
        ('TAGCACGTCC...', True),
        ('AGCACGTCCG...', True),
        ('GCACGTCCGA...', False),
    ]

    ax.text(5, 9.5, 'Syncmer filter:\nkeep ~75% of positions',
            ha='center', va='center', fontsize=9, color=DARK_GRAY, style='italic')

    y = 8.2
    for i, (kmer, keep) in enumerate(kmers):
        color = '#D5E8D4' if keep else '#FADBD8'
        edge = '#82B366' if keep else '#E74C3C'
        alpha = 1.0 if keep else 0.4
        rect = plt.Rectangle((1.0, y), 8, 0.55, facecolor=color,
                              edgecolor=edge, linewidth=1, alpha=alpha)
        ax.add_patch(rect)
        symbol = '\u2713' if keep else '\u2717'
        symbol_color = GREEN if keep else RED
        ax.text(0.6, y + 0.27, symbol, ha='center', va='center',
                fontsize=10, color=symbol_color, fontweight='bold')
        ax.text(5, y + 0.27, kmer, ha='center', va='center',
                fontsize=7.5, fontfamily='monospace', alpha=alpha)
        y -= 0.7

    ax.text(5, 2.5, f'Output: ~2.3B\nk-mers per\nhuman genome',
            ha='center', va='center', fontsize=10, color=NAVY, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#EBF5FB', edgecolor=BLUE))

    # Big arrow
    ax.annotate('', xy=(5, 3.5), xytext=(5, 4.2),
                arrowprops=dict(arrowstyle='->', color=NAVY, lw=2))

    # ---- Panel 2: Unsorted → buckets (MSD first pass) ----
    ax = axes[1]
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 12)
    ax.axis('off')
    ax.set_title('Step 2: Distribute by\nFirst Byte (256 Buckets)', fontsize=13,
                 fontweight='bold', color=NAVY, pad=10)

    # Unsorted pile
    unsorted_kmers = [
        ('TTAGCACG...', '#FDEDEC'),
        ('ACGTTAGC...', '#D4EFDF'),
        ('GCACGTCC...', '#F9E79F'),
        ('ACGTCCGA...', '#D4EFDF'),
        ('GTTAGCAC...', '#F9E79F'),
        ('TAGCACGT...', '#FDEDEC'),
        ('AGCACGTC...', '#D4EFDF'),
        ('CCGATAGC...', '#D6EAF8'),
    ]

    ax.text(3.5, 11.2, 'Unsorted k-mers', fontsize=10, fontweight='bold',
            color=DARK_GRAY, ha='center')
    y = 10.2
    for kmer, color in unsorted_kmers:
        rect = plt.Rectangle((0.5, y), 6, 0.45, facecolor=color,
                              edgecolor='#BDC3C7', linewidth=0.8)
        ax.add_patch(rect)
        ax.text(3.5, y + 0.22, kmer, ha='center', va='center',
                fontsize=7.5, fontfamily='monospace')
        y -= 0.55

    # Arrow to buckets
    ax.annotate('', xy=(8.5, 8.5), xytext=(7, 8.5),
                arrowprops=dict(arrowstyle='->', color=NAVY, lw=2))
    ax.text(7.7, 9.0, 'Sort by\n1st byte', ha='center', fontsize=9,
            color=NAVY, fontweight='bold')

    # Buckets
    bucket_data = [
        ('AA..', '#D4EFDF', ['ACGTTAGC...', 'ACGTCCGA...', 'AGCACGTC...']),
        ('CC..', '#D6EAF8', ['CCGATAGC...']),
        ('GG..', '#F9E79F', ['GCACGTCC...', 'GTTAGCAC...']),
        ('TT..', '#FDEDEC', ['TTAGCACG...', 'TAGCACGT...']),
    ]

    y = 7.5
    for prefix, color, entries in bucket_data:
        bh = len(entries) * 0.45 + 0.3
        rect = plt.Rectangle((7.5, y - bh + 0.3), 4, bh, facecolor=color,
                              edgecolor='#7F8C8D', linewidth=1.2, alpha=0.5)
        ax.add_patch(rect)
        ax.text(7.8, y + 0.1, prefix, fontsize=8, fontweight='bold', color=NAVY)
        ey = y - 0.3
        for entry in entries:
            ax.text(9.5, ey, entry, ha='center', va='center',
                    fontsize=6.5, fontfamily='monospace')
            ey -= 0.45
        y -= bh + 0.3

    ax.text(6, 1.5, 'Each bucket processed\nindependently\n(parallelizable!)',
            ha='center', va='center', fontsize=9, color=GREEN, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#EAFAF1', edgecolor=GREEN))

    # ---- Panel 3: Recursive sort within buckets ----
    ax = axes[2]
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 12)
    ax.axis('off')
    ax.set_title('Step 3: Recursively Sort\nEach Bucket (by next byte...)', fontsize=13,
                 fontweight='bold', color=NAVY, pad=10)

    # Show the AA bucket being recursively sorted
    ax.text(6, 11.2, 'Bucket "AA..." → sort by 2nd byte', fontsize=10,
            fontweight='bold', color=NAVY, ha='center')

    sub_buckets = [
        ('ACGT...', '#C8E6C9', ['ACGTCCGA...', 'ACGTTAGC...']),
        ('AGCA...', '#A5D6A7', ['AGCACGTC...']),
    ]

    y = 10.0
    for prefix, color, entries in sub_buckets:
        bh = len(entries) * 0.5 + 0.3
        rect = plt.Rectangle((1, y - bh + 0.3), 10, bh, facecolor=color,
                              edgecolor='#4CAF50', linewidth=1, alpha=0.5)
        ax.add_patch(rect)
        ax.text(1.3, y + 0.05, prefix, fontsize=8, fontweight='bold', color=NAVY)
        ey = y - 0.3
        for entry in entries:
            ax.text(6, ey, entry, ha='center', va='center',
                    fontsize=7.5, fontfamily='monospace')
            ey -= 0.5
        y -= bh + 0.4

    # Show recursion continuing
    ax.annotate('', xy=(6, 7.5), xytext=(6, 8.0),
                arrowprops=dict(arrowstyle='->', color=NAVY, lw=1.5))
    ax.text(6, 7.0, 'Continue sorting by 3rd, 4th, ... byte\nuntil bucket is small enough,\nthen Shell sort',
            ha='center', va='center', fontsize=9, color=DARK_GRAY, style='italic')

    # Show final sorted result
    ax.text(6, 5.5, 'Final sorted order (with LCP):', fontsize=10,
            fontweight='bold', color=NAVY, ha='center')

    sorted_entries = [
        ('ACGTCCGA...', 'LCP=0', '#C8E6C9'),
        ('ACGTTAGC...', 'LCP=32', '#A5D6A7'),
        ('AGCACGTC...', 'LCP=4', '#81C784'),
    ]

    y = 4.5
    for entry, lcp, color in sorted_entries:
        rect = plt.Rectangle((1.5, y), 6.5, 0.55, facecolor=color,
                              edgecolor='#388E3C', linewidth=1)
        ax.add_patch(rect)
        ax.text(4.75, y + 0.27, entry, ha='center', va='center',
                fontsize=8, fontfamily='monospace')
        ax.text(9, y + 0.27, lcp, ha='center', va='center',
                fontsize=8, color=PURPLE, fontweight='bold')
        y -= 0.7

    ax.text(6, 2.0, 'LCP = how many bases\nshared with previous entry\n(free byproduct of sorting!)',
            ha='center', va='center', fontsize=9, color=PURPLE, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#F4ECF7', edgecolor=PURPLE))

    # ---- Panel 4: Write to disk ----
    ax = axes[3]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis('off')
    ax.set_title('Step 4: Write\nSorted Index to Disk', fontsize=13,
                 fontweight='bold', color=NAVY, pad=10)

    # .gix file
    rect = plt.Rectangle((1, 9.5), 8, 1.5, facecolor='#FCF3CF',
                          edgecolor=ORANGE, linewidth=1.5)
    ax.add_patch(rect)
    ax.text(5, 10.5, '.gix (128 MB)', fontsize=10, fontweight='bold',
            ha='center', color=NAVY)
    ax.text(5, 9.9, 'Prefix jump table\n16M entries \u00d7 8 bytes', fontsize=8,
            ha='center', color=DARK_GRAY)

    # .ktab files
    ktab_y = 8.2
    for i in range(4):
        y = ktab_y - i * 1.5
        rect = plt.Rectangle((1, y), 8, 1.2, facecolor='#D6EAF8',
                              edgecolor=BLUE, linewidth=1.2)
        ax.add_patch(rect)
        if i < 3:
            ax.text(5, y + 0.6, f'.ktab.{i+1}', fontsize=9, fontweight='bold',
                    ha='center', color=NAVY)
            ax.text(5, y + 0.15, f'Sorted entries: suffix + LCP + pos + contig',
                    fontsize=7, ha='center', color=DARK_GRAY)
        else:
            ax.text(5, y + 0.6, f'... .ktab.N', fontsize=9, fontweight='bold',
                    ha='center', color=NAVY)
            ax.text(5, y + 0.15, f'N = 8\u201364 partitions', fontsize=7,
                    ha='center', color=DARK_GRAY)

    ax.text(5, 1.8, '~32 GB per\nhuman genome\n(~10 GB/Gbp)', fontsize=11,
            ha='center', va='center', color=RED, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#FDEDEC', edgecolor=RED))

    fig.suptitle('How FastGA Builds the Sorted K-mer Index (GIX)',
                 fontsize=17, fontweight='bold', color=NAVY, y=1.02)
    fig.tight_layout()
    outpath = os.path.join(OUTPUT_DIR, 'index_building.png')
    fig.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")


# ================================================================
# FIGURE 2: Linear Merge (Adaptamer Matching)
# ================================================================

def plot_linear_merge():
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 15)
    ax.axis('off')

    # Title
    ax.text(12, 14.5, 'How the Linear Merge Finds Seeds (Adaptamer Matching)',
            ha='center', fontsize=17, fontweight='bold', color=NAVY)

    # ---- Draw Table A (left) ----
    ta_x, ta_y = 1, 1.5
    ta_w, entry_h = 5.5, 0.65
    ax.text(ta_x + ta_w/2, 13.5, 'Genome A (sorted)', fontsize=13,
            fontweight='bold', ha='center', color=NAVY)
    ax.text(ta_x + ta_w/2, 13.0, 'Pointer \u25B6 advances top-to-bottom',
            fontsize=9, ha='center', color=DARK_GRAY, style='italic')

    entries_a = [
        ('AACCCGGGTTTAACCC...', 0, False, ''),
        ('AACCCGGGTTTAGCCC...', 36, False, ''),
        ('AACCCTTTTTTAAGGG...', 20, False, ''),
        ('ACGTTAGCACGTCCGA...', 0, True, '\u25B6 Current'),
        ('ACGTTAGCACGTTTAG...', 36, False, ''),
        ('ACGTTAGCCCGTAAAC...', 32, False, ''),
        ('AGCACGTCCGATAGCG...', 0, False, ''),
        ('CCGATAGCGTTTACCC...', 0, False, ''),
        ('GCACGTCCGATAGCGT...', 0, False, ''),
        ('GTTAGCACGTCCGATA...', 0, False, ''),
        ('TTAGCACGTCCGATAG...', 0, False, ''),
    ]

    y = 12.2
    for i, (kmer, lcp, is_current, label) in enumerate(entries_a):
        if is_current:
            color = '#A9DFBF'
            edge = GREEN
            lw = 2.5
        elif i < 3:
            color = '#E8E8E8'
            edge = '#CCC'
            lw = 0.5
        else:
            color = '#D5E8D4'
            edge = '#82B366'
            lw = 0.8
        rect = plt.Rectangle((ta_x, y), ta_w, entry_h - 0.05,
                              facecolor=color, edgecolor=edge, linewidth=lw)
        ax.add_patch(rect)
        # Show k-mer text, highlight matching prefix for current
        ax.text(ta_x + 0.15, y + entry_h/2 - 0.02, kmer[:18] + '...',
                va='center', fontsize=7.5, fontfamily='monospace')
        if label:
            ax.text(ta_x - 0.3, y + entry_h/2, label, ha='right', va='center',
                    fontsize=9, color=GREEN, fontweight='bold')
        y -= entry_h

    # ---- Draw Table B (right) ----
    tb_x = 10.5
    ax.text(tb_x + ta_w/2, 13.5, 'Genome B (sorted)', fontsize=13,
            fontweight='bold', ha='center', color=NAVY)
    ax.text(tb_x + ta_w/2, 13.0, 'Scanned to find matching block',
            fontsize=9, ha='center', color=DARK_GRAY, style='italic')

    entries_b = [
        ('AACCCGGGTTTAACCC...', 0, False, '', None),
        ('AACCCTTTTTTAAGGG...', 20, False, '', None),
        ('ACGTTAGCACGTAACC...', 0, True, 'LCP=40', GREEN),
        ('ACGTTAGCACGTCCTT...', 36, True, 'LCP=36', '#2ECC71'),
        ('ACGTTAGCACTTGGGA...', 32, True, 'LCP=32', '#58D68D'),
        ('ACGTTAGCCCAATTTG...', 28, False, 'LCP=28\n(block ends)', ORANGE),
        ('AGCACGTCCGATAGCG...', 0, False, '', None),
        ('CCGATAGCGTTTACCC...', 0, False, '', None),
        ('GCACGTCCGATAGCGT...', 0, False, '', None),
        ('GTTAGCACGTCCGATA...', 0, False, '', None),
        ('TTAGCACGTCCGATAG...', 0, False, '', None),
    ]

    y = 12.2
    for i, (kmer, lcp, is_match, label, label_color) in enumerate(entries_b):
        if is_match:
            color = '#AED6F1'
            edge = BLUE
            lw = 2.0
        elif i < 2:
            color = '#E8E8E8'
            edge = '#CCC'
            lw = 0.5
        else:
            color = '#D6EAF8'
            edge = '#85C1E9'
            lw = 0.8
        rect = plt.Rectangle((tb_x, y), ta_w, entry_h - 0.05,
                              facecolor=color, edgecolor=edge, linewidth=lw)
        ax.add_patch(rect)
        ax.text(tb_x + 0.15, y + entry_h/2 - 0.02, kmer[:18] + '...',
                va='center', fontsize=7.5, fontfamily='monospace')
        if label and label_color:
            ax.text(tb_x + ta_w + 0.2, y + entry_h/2, label,
                    va='center', fontsize=8.5, color=label_color, fontweight='bold')
        y -= entry_h

    # ---- Draw matching arrows ----
    a_current_y = 12.2 - 3 * entry_h + entry_h/2
    b_match_ys = [12.2 - i * entry_h + entry_h/2 for i in [2, 3, 4]]

    for by in b_match_ys:
        ax.annotate('', xy=(tb_x, by), xytext=(ta_x + ta_w, a_current_y),
                    arrowprops=dict(arrowstyle='->', color=GREEN, lw=1.5,
                                    connectionstyle='arc3,rad=0.1'))

    # ---- Matching prefix highlight box ----
    ax.text(8.2, 13.5, 'Match!', fontsize=14, fontweight='bold',
            color=GREEN, ha='center')

    # ---- Explanation box: How LCP drives the merge ----
    box_x, box_y = 17.5, 12.5
    explanation = [
        ("1. A\u2019s pointer is at ACGTTAGCACGT...", NAVY),
        ("", NAVY),
        ("2. Scan B for entries sharing prefix:", NAVY),
        ("   B[2]: ACGTTAGCACGT... \u2192 LCP=40 (exact!)", GREEN),
        ("   B[3]: ACGTTAGCACGT... \u2192 LCP=36", '#2ECC71'),
        ("   B[4]: ACGTTAGCACTT... \u2192 LCP=32", '#58D68D'),
        ("   B[5]: ACGTTAGCCC...   \u2192 LCP=28 \u2190 STOP", ORANGE),
        ("", NAVY),
        ("3. Block = B[2..4], longest prefix = 40", NAVY),
        ("   Block size = 3 entries (\u2264 freq 10)", NAVY),
        ("", NAVY),
        ("4. Emit 3 seed pairs:", GREEN),
        ("   (A:pos100,chr1) \u2194 (B:pos200,chr3)", GREEN),
        ("   (A:pos100,chr1) \u2194 (B:pos800,chr5)", GREEN),
        ("   (A:pos100,chr1) \u2194 (B:pos150,chr2)", GREEN),
        ("", NAVY),
        ("5. Advance A\u2019s pointer to next entry", NAVY),
    ]

    for i, (text, color) in enumerate(explanation):
        ax.text(box_x, box_y - i * 0.55, text, fontsize=8.5,
                fontfamily='monospace', color=color, va='center')

    # Border around explanation
    rect = plt.Rectangle((17, 3.3), 6.8, 10.0, facecolor='#FAFAFA',
                          edgecolor=DARK_GRAY, linewidth=1, linestyle='--',
                          zorder=0)
    ax.add_patch(rect)
    ax.text(20.4, 13.1, 'Merge Logic', fontsize=11, fontweight='bold',
            color=NAVY, ha='center')

    # ---- Bottom: key insight box ----
    ax.text(8, 0.5,
            'Key: Both tables are sorted, so matching entries form a contiguous block in B.\n'
            'LCP values mark block boundaries \u2014 no byte-by-byte comparison needed.\n'
            'One forward pass through both tables finds ALL seeds. Total: 5.3s for two human genomes.',
            ha='center', va='center', fontsize=11, color=NAVY,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#EBF5FB',
                      edgecolor=BLUE, alpha=0.9))

    fig.tight_layout()
    outpath = os.path.join(OUTPUT_DIR, 'linear_merge.png')
    fig.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")


if __name__ == '__main__':
    plot_index_building()
    plot_linear_merge()
    print("\nAll figures generated!")
