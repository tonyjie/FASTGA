#!/usr/bin/env python3
"""Generate benchmark performance figures for FastGA slides."""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

OUTPUT_DIR = "/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/docs/slides"

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 13,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'legend.fontsize': 12,
    'figure.facecolor': 'white',
    'axes.facecolor': '#FAFAFA',
    'axes.grid': True,
    'grid.alpha': 0.3,
})

COLORS = {
    'gdb':        '#95A5A6',  # gray
    'gix':        '#3498DB',  # blue
    'seed_merge': '#E67E22',  # orange
    'sort_align': '#E74C3C',  # red
    'paf':        '#2ECC71',  # green
    'total':      '#2C3E50',  # dark
    'ideal':      '#BDC3C7',  # light gray dashed
}


# ── Data ───────────────────────────────────────────────────────────────────

threads = [1, 2, 4, 8, 16, 32]

# Overall
wall_overall = [133.5, 73.3, 49.2, 32.6, 25.1, 22.0]
speedup_overall = [1.00, 1.82, 2.71, 4.09, 5.32, 6.07]

# Per-phase wall clock (rep1)
wall_gdb       = [1.0, 1.1, 1.0, 1.0, 1.0, 1.0]
wall_gix       = [18.9, 11.1, 7.5, 6.2, 6.2, 7.7]
wall_merge     = [30.0, 16.5, 11.2, 8.3, 7.1, 5.7]
wall_sort_aln  = [83.1, 43.9, 28.7, 16.6, 9.8, 7.1]
wall_paf       = [0.5, 0.3, 0.2, 0.2, 0.2, 0.2]

# Per-phase speedup
speedup_gix       = [1.00, 1.70, 2.53, 3.04, 3.06, 2.46]
speedup_merge     = [1.00, 1.83, 2.69, 3.64, 4.23, 5.25]
speedup_sort_aln  = [1.00, 1.90, 2.90, 5.01, 8.46, 11.79]

# Peak RSS (MB)
rss = [513, 522, 552, 653, 682, 740]


# ── Figure 1: Overall Thread Scaling (Wall Clock + Speedup) ───────────────

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

# Left: Wall clock
ax1.plot(threads, wall_overall, 'o-', color=COLORS['total'], linewidth=2.5,
         markersize=8, label='Wall clock')
for t, w in zip(threads, wall_overall):
    ax1.annotate(f'{w:.1f}s', (t, w), textcoords="offset points",
                 xytext=(0, 12), ha='center', fontsize=10)
ax1.set_xlabel('Threads')
ax1.set_ylabel('Wall Clock Time (seconds)')
ax1.set_title('Overall Wall Clock vs Threads')
ax1.set_xscale('log', base=2)
ax1.set_xticks(threads)
ax1.set_xticklabels(threads)
ax1.set_ylim(0, 150)

# Right: Speedup
ax2.plot(threads, speedup_overall, 'o-', color=COLORS['sort_align'], linewidth=2.5,
         markersize=8, label='Actual speedup')
ax2.plot(threads, threads, '--', color=COLORS['ideal'], linewidth=1.5,
         label='Ideal (linear)')
for t, s in zip(threads, speedup_overall):
    ax2.annotate(f'{s:.1f}x', (t, s), textcoords="offset points",
                 xytext=(0, 12), ha='center', fontsize=10)
ax2.set_xlabel('Threads')
ax2.set_ylabel('Speedup')
ax2.set_title('Overall Speedup vs Threads')
ax2.set_xscale('log', base=2)
ax2.set_xticks(threads)
ax2.set_xticklabels(threads)
ax2.set_ylim(0, 35)
ax2.legend(loc='upper left')

fig.suptitle('FastGA Thread Scaling (EXAMPLE Dataset, ~86 Mbp)', fontsize=17, y=1.02)
fig.tight_layout()
fig.savefig(f'{OUTPUT_DIR}/thread_scaling_overall.png', dpi=200, bbox_inches='tight')
print("Saved thread_scaling_overall.png")


# ── Figure 2: Per-Phase Stacked Bar Chart ─────────────────────────────────

fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(threads))
width = 0.6

bottoms = np.zeros(len(threads))
for label, data, color in [
    ('GDB', wall_gdb, COLORS['gdb']),
    ('GIX Build', wall_gix, COLORS['gix']),
    ('Seed Merge', wall_merge, COLORS['seed_merge']),
    ('Sort + Align', wall_sort_aln, COLORS['sort_align']),
    ('PAF Conv', wall_paf, COLORS['paf']),
]:
    bars = ax.bar(x, data, width, bottom=bottoms, label=label, color=color,
                  edgecolor='white', linewidth=0.5)
    bottoms += np.array(data)

# Add total labels on top
for i, total in enumerate(bottoms):
    ax.text(i, total + 1.5, f'{total:.0f}s', ha='center', fontsize=11, fontweight='bold')

ax.set_xlabel('Threads')
ax.set_ylabel('Wall Clock Time (seconds)')
ax.set_title('Per-Phase Runtime Breakdown by Thread Count')
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.legend(loc='upper right', framealpha=0.9)
ax.set_ylim(0, 155)

fig.tight_layout()
fig.savefig(f'{OUTPUT_DIR}/thread_scaling_stacked.png', dpi=200, bbox_inches='tight')
print("Saved thread_scaling_stacked.png")


# ── Figure 3: Per-Phase Speedup Comparison ────────────────────────────────

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(threads, threads, '--', color=COLORS['ideal'], linewidth=1.5, label='Ideal (linear)')
ax.plot(threads, speedup_gix, 's-', color=COLORS['gix'], linewidth=2, markersize=7,
        label='GIX Build')
ax.plot(threads, speedup_merge, '^-', color=COLORS['seed_merge'], linewidth=2, markersize=7,
        label='Seed Merge')
ax.plot(threads, speedup_sort_aln, 'o-', color=COLORS['sort_align'], linewidth=2.5,
        markersize=8, label='Sort + Align')

# Annotate sort+align at T=32
ax.annotate('11.8x', (32, 11.79), textcoords="offset points",
            xytext=(-30, 8), fontsize=11, fontweight='bold', color=COLORS['sort_align'])
# Annotate seed merge at T=32
ax.annotate('5.3x', (32, 5.25), textcoords="offset points",
            xytext=(-25, 8), fontsize=11, color=COLORS['seed_merge'])
# Annotate GIX regression at T=32
ax.annotate('2.5x\n(regresses!)', (32, 2.46), textcoords="offset points",
            xytext=(-50, -30), fontsize=10, color=COLORS['gix'],
            arrowprops=dict(arrowstyle='->', color=COLORS['gix'], lw=1.2))

ax.set_xlabel('Threads')
ax.set_ylabel('Speedup')
ax.set_title('Per-Phase Speedup Comparison')
ax.set_xscale('log', base=2)
ax.set_xticks(threads)
ax.set_xticklabels(threads)
ax.set_ylim(0, 35)
ax.legend(loc='upper left')

fig.tight_layout()
fig.savefig(f'{OUTPUT_DIR}/thread_scaling_per_phase.png', dpi=200, bbox_inches='tight')
print("Saved thread_scaling_per_phase.png")


# ── Figure 4: Phase Proportion Shift ──────────────────────────────────────

fig, ax = plt.subplots(figsize=(10, 5.5))

# Percentages
pct_gdb      = [0.7, 1.5, 2.0, 3.0, 4.1, 4.6]
pct_gix      = [14.1, 15.3, 15.5, 19.3, 25.4, 35.4]
pct_merge    = [22.5, 22.6, 23.0, 25.6, 29.2, 26.4]
pct_sort_aln = [62.3, 60.3, 59.1, 51.4, 40.4, 32.6]
pct_paf      = [0.4, 0.4, 0.4, 0.6, 0.8, 1.0]

x = np.arange(len(threads))
width = 0.6

bottoms = np.zeros(len(threads))
for label, data, color in [
    ('GDB', pct_gdb, COLORS['gdb']),
    ('GIX Build', pct_gix, COLORS['gix']),
    ('Seed Merge', pct_merge, COLORS['seed_merge']),
    ('Sort + Align', pct_sort_aln, COLORS['sort_align']),
    ('PAF Conv', pct_paf, COLORS['paf']),
]:
    ax.bar(x, data, width, bottom=bottoms, label=label, color=color,
           edgecolor='white', linewidth=0.5)
    bottoms += np.array(data)

ax.set_xlabel('Threads')
ax.set_ylabel('Percentage of Total Wall Time')
ax.set_title('Phase Proportion Shift with Thread Count')
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.set_ylim(0, 105)
ax.legend(loc='upper right', framealpha=0.9)
ax.yaxis.set_major_formatter(ticker.PercentFormatter())

fig.tight_layout()
fig.savefig(f'{OUTPUT_DIR}/thread_scaling_proportions.png', dpi=200, bbox_inches='tight')
print("Saved thread_scaling_proportions.png")


# ── Figure 5: EXAMPLE vs Human Comparison ─────────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Human per-phase breakdown (pie chart)
human_phases = ['GDB\n18s (3%)', 'GIX Build\n82s (14%)',
                'Seed Merge\n5.3s (1%)', 'Sort + Align\n482s (82%)']
human_vals = [18.3, 81.8, 5.3, 481.9]
human_colors = [COLORS['gdb'], COLORS['gix'], COLORS['seed_merge'], COLORS['sort_align']]

axes[0].pie(human_vals, labels=human_phases, colors=human_colors,
            autopct='', startangle=140, textprops={'fontsize': 11})
axes[0].set_title('Human Genome\nPhase Breakdown (T=32)', fontsize=14)

# EXAMPLE per-phase breakdown (pie chart)
example_phases = ['GDB\n1.0s (5%)', 'GIX Build\n7.7s (35%)',
                  'Seed Merge\n5.7s (26%)', 'Sort + Align\n7.1s (33%)']
example_vals = [1.0, 7.7, 5.7, 7.1]

axes[1].pie(example_vals, labels=example_phases, colors=human_colors,
            autopct='', startangle=140, textprops={'fontsize': 11})
axes[1].set_title('EXAMPLE Dataset\nPhase Breakdown (T=32)', fontsize=14)

# Scale comparison (bar chart)
categories = ['Wall\nClock', 'Seeds', 'Peak\nRSS', 'Align-\nments']
example_vals_bar = [22.0, 51, 0.74, 323]
human_vals_bar = [587.8, 1298, 19.0, 518]
scale_labels = ['27x', '25x', '26x', '1.6x']

x = np.arange(len(categories))
w = 0.35
bars1 = axes[2].bar(x - w/2, example_vals_bar, w, label='EXAMPLE (86 Mbp)',
                     color=COLORS['gix'], alpha=0.8)
bars2 = axes[2].bar(x + w/2, human_vals_bar, w, label='Human (3.1 Gbp)',
                     color=COLORS['sort_align'], alpha=0.8)

# Scale annotations
for i, (ev, hv, sl) in enumerate(zip(example_vals_bar, human_vals_bar, scale_labels)):
    axes[2].annotate(sl, (i + w/2, hv), textcoords="offset points",
                     xytext=(0, 5), ha='center', fontsize=10, fontweight='bold')

axes[2].set_ylabel('Value (mixed units)')
axes[2].set_title('EXAMPLE vs Human\nScale Comparison', fontsize=14)
axes[2].set_xticks(x)
axes[2].set_xticklabels(categories)
axes[2].legend(fontsize=10)
axes[2].set_yscale('log')

fig.suptitle('FastGA Performance: EXAMPLE vs Human Genome', fontsize=17, y=1.02)
fig.tight_layout()
fig.savefig(f'{OUTPUT_DIR}/benchmark_comparison.png', dpi=200, bbox_inches='tight')
print("Saved benchmark_comparison.png")

print("\nAll figures generated!")
