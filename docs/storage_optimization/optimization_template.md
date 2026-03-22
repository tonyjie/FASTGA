# Optimization N: [Name]

## Summary

[One paragraph: what the optimization does and why it reduces storage.]

| Property | Value |
|---|---|
| **Type** | [Zero-cost / Trade-off] |
| **Quality tier** | [Tier 1 (Bit-exact) / Tier 2 (Comparable)] |
| **Target** | [What storage component is reduced] |
| **Code changes** | [Which files, rough scope] |

## What Changed

[Description of the code change. Include key design decisions and why.]

## Evaluation Results

### 1. Output Verification

**For Tier 1 (Bit-exact):**

| Check | Result |
|---|---|
| `.1aln` binary diff | [identical / differs by N bytes (header only)] |
| ONEview text diff (skip header) | [zero differences / N differences] |
| Alignment count | [baseline] vs [optimized] |

**For Tier 2 (Comparable):**

| Metric | Baseline | Optimized | Change |
|---|---|---|---|
| Non-redundant alignments | | | |
| Average alignment length | | | |
| Genome coverage (bp) | | | |
| [Other relevant metrics] | | | |

### 2. Storage Impact

**EXAMPLE dataset (~86 Mbp per genome, T=32):**

| Metric | Baseline | Optimized | Change |
|---|---:|---:|---|
| Peak total disk | | | |
| Persistent files (after run) | | | |
| Peak temp files | | | |
| Disk during sort+align | | | |

[Include storage timeline comparison figure if applicable.]

**Human genome (GRCh38 vs CHM13, T=32) — if tested:**

| Metric | Baseline | Optimized | Change |
|---|---:|---:|---|
| Peak total disk | | | |
| Persistent files (after run) | | | |
| Peak temp files | | | |
| Disk during sort+align | | | |

### 3. Performance Impact

**EXAMPLE dataset (T=32):**

| Metric | Baseline | Optimized | Change |
|---|---:|---:|---|
| Wall clock | | | |
| User CPU | | | |
| System CPU | | | |
| Peak RSS | | | |

Per-phase wall clock:

| Phase | Baseline | Optimized | Change |
|---|---:|---:|---|
| GDB (both) | | | |
| GIX build (both) | | | |
| Seed merge | | | |
| Sort + align | | | |

**Human genome (T=32) — if tested:**

| Metric | Baseline | Optimized | Change |
|---|---:|---:|---|
| Wall clock | | | |
| User CPU | | | |
| Peak RSS | | | |

## Git Info

| | Commit | Description |
|---|---|---|
| Docs | `<hash>` | Optimization N documentation |
| Code | `<hash>` | Code change description |

To revert this optimization's code change while keeping docs:
```bash
git revert <code_commit_hash>
```

## Checklist

- [ ] Code change implemented
- [ ] Builds without new warnings
- [ ] Output verification: Tier 1 bit-exact OR Tier 2 comparable quality
- [ ] Storage timeline comparison figure generated
- [ ] Performance comparison (runtime) measured
- [ ] Tested on EXAMPLE dataset
- [ ] Tested on human genome dataset (if applicable)
- [ ] Results documented in this file
- [ ] README.md status table updated
