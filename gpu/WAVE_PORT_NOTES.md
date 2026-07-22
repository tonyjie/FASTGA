# WAVE_PORT_NOTES: `align.c`'s wave contract for the GPU port

This is the oracle document for porting FastGA's CPU wave aligner (`forward_wave` /
`reverse_wave` in `align.c`) to a GPU kernel. Every claim below cites an exact `align.c`
line (or, where the wave depends on data conventions defined elsewhere, the defining line
in `align.h`, `GDB.c`, `gene_core.c`, or `FastGA.c`). All line numbers are against the
`agent-optimization-wt` worktree as of this writing; re-verify against `main` before using
this as ground truth for a port against a different branch.

The two functions are near-mirror images of each other:

- `forward_wave` (`align.c:385–916`) extends **rightward** from the seed, growing
  `(aepos,bepos)` — the suffix of the local alignment.
- `reverse_wave` (`align.c:919–1459`) extends **leftward** from the seed, growing
  `(abpos,bbpos)` — the prefix of the local alignment.

They are called back-to-back by `Local_Alignment` (`align.c:1464–1620`), which is the only
public entry point covered here (`align.h:205–236` is its contract).

---

## Step 1 — the CPU wave contract

### (a) Inputs

`forward_wave`'s signature (`align.c:385–386`):

```c
static int forward_wave(_Work_Data *work, _Align_Spec *spec, Alignment *align,
                        int *mind, int maxd, int mida, int minp, int maxp, int aoff)
```

`reverse_wave`'s signature is identical except `mind` is passed **by value**, not by
pointer (`align.c:919–920`), because only the forward pass needs to hand its final low
diagonal back to the caller (see "outputs" below).

- **`aseq` / `bseq`** (`align.c:387–388` forward, `align.c:921–922` reverse): numeric
  (0=A,1=C,2=G,3=T) sequence buffers pulled off `align->aseq` / `align->bseq`. Both buffers
  are **sentinel-terminated with the byte value `4`** at both ends — this is the
  `New_Contig_Buffer` (`GDB.c:1722–1731`, allocates `maxctg+8` and returns `contig+1` so
  `buffer[-1]` is addressable) + `Get_Contig` (`GDB.c:1743–1842`) convention:
  `buffer[-1] = 4` is set at `GDB.c:1807/1831/1910` and the trailing `buffer[len] = 4` is
  set inside `Uncompress_Read` at `gene_core.c:397` (also directly at `GDB.c:1898` for the
  `Get_Contig_Piece` path). `reverse_wave` additionally offsets both pointers by `-1`
  (`align.c:921–922`, `aseq = align->aseq - 1`) because it walks backward and indexes with
  `x -= 1`. The `4` sentinel is what the `c==4`/`d==4` clip checks below detect.
- **Band `[*mind, maxd]`** — the diagonal band to search, in `low`/`hgh` (`align.c:417–418`
  forward, `align.c:951–952` reverse). Diagonal `k = a - b` (position space, not index
  space).
- **Seed anti-diagonal `mida`** (`align.c:386,417`) — `a + b` at the seed; the 0-wave is
  computed by scanning every diagonal `k` in `[low,hgh]` from this single anti-diagonal
  (`align.c:462–466`: `x = (mida+k)>>1`, i.e. the seed point restricted to diagonal `k`).
- **`aoff`** (`align.c:386`, used at `align.c:482` forward / `align.c:1010` reverse) —
  trace-point boundary phase offset. Computed once in `Local_Alignment`
  (`align.c:1524–1527`): `aoff = align->alen % spec->trace_space` if `ACOMP(align->flags)`
  else `0`. It shifts where the `tspace`-periodic checkpoint boundaries (`na`) land when the
  A-sequence has been complemented, since boundaries are anchored to `align->alen`, not to
  0.
- **`minp`/`maxp`** (`align.c:386`) — hard anti-diagonal-space clamp on how far the band can
  widen (`align.c:655–667` forward, `align.c:1179–1191` reverse); set by `Local_Alignment`
  (`align.c:1507–1522`) mainly to stop a self-comparison (`selfie = align->aseq==align->bseq`,
  `align.c:1502`) from reaching across the `k=0` diagonal.
- **`Align_Spec`** fields actually read by the waves (`align.c:403–407` forward,
  `align.c:937–941` reverse): `tspace = spec->trace_space` (checkpoint/trace-point period —
  FastGA always builds it as `100`, see `FastGA.c:87` `#define TSPACE 100` and
  `FastGA.c:3997` `New_Align_Spec(1.-ALIGN_RATE,100,gdb1->freq,0)`, so `reach=0` in
  practice); `PATH_AVE = spec->ave_path` (min recent-match count gate for trimming,
  computed at `align.c:255` from `ave_corr`); `REACH = spec->reach`; `SCORE`/`TABLE =
  spec->score/spec->table` (suffix-positivity lookup tables built by `set_table`,
  `align.c:211–222`, and populated in `New_Align_Spec`, `align.c:255–269`).

### (b) Per-diagonal furthest-reach inner loop

The 0-wave (initial furthest-reach per diagonal from the seed anti-diagonal) is
`align.c:462–548` (forward) / `align.c:990–1076` (reverse). Steady-state waves repeat the
same slide inside the `dif`-loop at `align.c:723–743` (forward) / `align.c:1243–1263`
(reverse).

The slide itself (forward, `align.c:494–511`):

```c
while (1)
  { c = bs[x];
    if (c == 4) { more = 0; if (bclip < k) bclip = k; break; }   // B-sequence end reached
    d = aseq[x];
    if (c != d) { if (d == 4) { more = 0; aclip = k; } break; }   // mismatch, or A-end
    x += 1;
  }
```

i.e. "slide while `c==d`" (match), stop on mismatch, and specially flag exhaustion of one
sequence via the `4` sentinel (`c==4` → B-clip at this diagonal `bclip`; `d==4` on a
mismatch path → A-clip `aclip`). `bs = bseq - hgh` is re-based per outer `k` so that
`bs[x] == bseq[x-k]` for the current diagonal (`align.c:465,547` — `bs += 1` after each `k`
as `k` decrements). The reverse version slides `x -= 1` instead (`align.c:1037`,
`align.c:1259`) and tests `bclip`/`aclip` with the opposite comparison direction
(`align.c:1025`, `align.c:1247`) since it is extending toward decreasing coordinates.

`c = (x<<1) - k` (`align.c:512,744`; reverse identical form at `align.c:1039,1264`) converts
the slide's stopping `x` (an A-position) back into the `besta`-comparable "anti-diagonal
value" `V[k]` used for furthest-reach comparisons.

### (c) X-drop / trim rule

Constants (`align.c:166–180`):

- `TRIM_LEN = 15`(the trim decision window is 2× this = 30 columns; align.c:166-168
  comment), `PATH_LEN = 60` (align.c:172, the full recent-alignment bitvector `BVEC T`
  tracked per diagonal), `TRIM_MASK = 0x7fff` (`(1<<TRIM_LEN)-1`, align.c:178).
- **`TRIM_MLAG = 250`** (`align.c:179`) — how far the last-accepted-trim point (`lasta`) is
  allowed to lag behind the best point (`besta`) before the wave loop gives up. Used
  verbatim as the wave's continuation condition:
  - forward: `while (more && lasta >= besta - TRIM_MLAG)` at **`align.c:583`**
  - reverse: `while (more && lasta <= besta + TRIM_MLAG)` at **`align.c:1108`**
    (mirrored inequality since reverse's `besta` decreases as it extends).
- `WAVE_LAG = 70` (`align.c:180`) is a *separate* per-wave band-pruning threshold, not the
  trim rule: at the end of each wave iteration, diagonals whose furthest-reach `V[k]` has
  fallen more than `WAVE_LAG` behind `besta` are dropped from the live band
  (`align.c:823–831` forward, `align.c:1343–1351` reverse).

The trim *decision* — whether the current furthest-reach point is accepted as a new
"trim point" (`trima/trimx/trimd/trimha`, the coordinates that will actually be reported as
the alignment's endpoint, as opposed to `besta`, which is just the raw furthest-reach
optimum) — is the suffix-positivity test done at every new `besta` (forward,
`align.c:770–783`; reverse, `align.c:1290–1303`):

```c
if (c > besta)
  { besta = c; bestx = x;
    if (m >= PATH_AVE)                                    // enough recent matches
      { lasta = c;
        if (TABLE[b & TRIM_MASK] >= 0)                    // last 15 cols suffix-positive
          if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)  // last 30
            { trima = c; trimx = x; trimd = dif; trimha = ha; }
```

`b` is the 60-bit shifting match/mismatch bitvector `BVEC T` for this diagonal (1 =
match column, `align.c:718–720,740–742`); `m` is `M[k]`, a running match-quality count
gated against `PATH_AVE = spec->ave_path` (`align.c:255`). `TABLE`/`SCORE` are the
suffix-positivity lookup tables from `set_table` (`align.c:211–222`), built so that
`TABLE[x]` is non-negative iff the 15-column suffix encoded by bit pattern `x` has
non-negative score under match=+`mscore`/mismatch=-`dscore` — this operationalizes the
"last `2*TRIM_LEN` edits are prefix/suffix-positive" comment at `align.c:166–168`.

If the wave loop exits with `more==0` and `REACH` set (`align.c:852` forward,
`align.c:1372` reverse) and a valid `morem` was recorded, the reported endpoint is instead
overridden with the `morea/morex/mored/moreha` "extend-to-boundary" point recorded when a
clip (`aclip`/`bclip`, one sequence ran out) was hit — see `align.c:551–574,796–821` and
`align.c:1078–1101,1316–1341` for where `more{a,x,d,ha}` get set.

### (d) `Pebble{ptr,diag,diff,mark}` checkpoints and the backward walk

```c
typedef struct
  { int ptr; int diag; int diff; int mark; } Pebble;   // align.c:347–352
```

`cells[]` (`Pebble *cells`, e.g. `align.c:400`) is a growable array of checkpoints, one
pushed per diagonal each time its furthest-reach point crosses a `tspace`-periodic
A-coordinate boundary (`na`). Fields: `ptr` = index of the *previous* checkpoint on the same
diagonal's lineage (`-1` = none, i.e. linked-list back-pointer); `diag` = the diagonal `k` at
checkpoint time; `diff` = the wave number `dif` (cumulative edit distance) at checkpoint
time; `mark` = the boundary A-value (`na`) this checkpoint represents.

Checkpoints are created in two places:

- **0-wave root checkpoint**, one per diagonal, *before* the slide, recording the boundary
  just below the seed (`align.c:482–491` forward: `na = ((x+(tspace-aoff))/tspace-1)*tspace
  +aoff; pb->mark = na;`); then **0-wave crossing checkpoints**, one per boundary the slide
  actually crosses (`align.c:514–533`, `while (x >= na) { ...push...; na += tspace; }`).
  Reverse computes the same `na` boundary at `align.c:1010` but its root checkpoint diverges
  from forward's: `align.c:1018` sets `pb->mark = x;` — the raw seed A-position, not the
  `na` boundary value — while everything else (both passes' crossing checkpoints and
  steady-state pushes: `align.c:490` forward root, `align.c:1057` reverse crossings) uses
  `na`/`NA[k]`. Crossings are at `align.c:1041–1060`, with `na` counting down.
- **Steady-state wave checkpoints** (`align.c:746–768` forward, `align.c:1266–1288`
  reverse): same `while (x >= NA[k])` boundary-crossing push, but guarded by
  `cells[ha].mark < NA[k]` (forward) / `cells[ha].mark > NA[k]` (reverse) so a checkpoint is
  only pushed if the inherited chain doesn't already record this boundary — avoiding
  duplicate pushes when a diagonal's lineage is unchanged from the previous wave.

`cells[]` grows dynamically via `Realloc` whenever `avail >= cmax` (forward:
`align.c:473–480,515–521,748–756`; reverse: `align.c:1001–1008,1042–1049,1268–1276`) — on
the CPU this is an unbounded amortized-growth heap array shared across the whole wave loop.

**Backward walk emitting trace-points** (forward, `align.c:860–895`): first the
pointer chain rooted at `trimha` is reversed in place so it can be walked *forward* in time
(`align.c:863–869`):

```c
a = -1;
for (h = trimha; h >= 0; h = b)
  { b = cells[h].ptr; cells[h].ptr = a; a = h; }
h = a;                                    // h now = earliest checkpoint in the lineage
```

then walked from earliest to latest, emitting one `(diff-delta, b-delta)` pair per
checkpoint into `apath->trace` (`align.c:871–891`):

```c
k = cells[h].diag; b = (mida-k)>>1; e = 0; low = k;
for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
  { k = cells[h].diag; a = cells[h].mark - k; d = cells[h].diff;
    atrace[atlen++] = (uint16)(d-e);      // diffs since previous trace point
    atrace[atlen++] = (uint16)(a-b);      // B-delta since previous trace point
    b = a; e = d;
  }
```

with a final partial segment out to `(trimx,trimy)` appended/merged at `align.c:892–905`.
This matches the trace-point semantics documented at `align.h:58–77`. Reverse
(`align.c:1366–1449`) does the mirror-image walk but writes **backward** into the trace
buffer using `--atlen` (`align.c:1412–1428`) because it must prepend to the trace segment
the forward pass already wrote — `apath->trace` was originally centered at
`work->points + maxtp` (`align.c:1495`) precisely so the reverse pass has room to grow
backward.

### (e) Outputs written to `apath`

Forward (`align.c:907–913`, `Path *apath = align->path`, `align.c:389`):
`apath->aepos = trimx; apath->bepos = trimy; apath->diffs = trimd; apath->tlen = atlen;`
and `*mind = low` (`align.c:913`) — the caller-visible band low is updated to the final
diagonal reached, which `Local_Alignment` then feeds as the *single* diagonal
(`low,low`) for the subsequent `reverse_wave` call (`align.c:1543`:
`reverse_wave(work,spec,align,low,low,anti,minp,maxp,aoff)`).

Reverse (`align.c:1451–1456`): `apath->abpos = trimx; apath->bbpos = trimy; apath->diffs =
apath->diffs + trimd;` (accumulates onto the forward pass's diff count) `apath->tlen =
apath->tlen - atlen; apath->trace = atrace + atlen;` (extends the trace buffer backward and
re-points `trace` at the new, earlier start).

`Path`/`Alignment` struct field layout: `align.h:89–95` (`Path`), `align.h:145–152`
(`Alignment`). The public contract for `Local_Alignment` itself is `align.h:205–233`; the
prototype is `align.h:235–236`.

---

## Step 2 — GPU port mapping

**Warp mapping.** One warp (32 lanes) is assigned to one active band `[low,hgh]` at a given
wave (`dif`). Lanes tile the live diagonals: the CPU already keeps the band narrow via the
`WAVE_LAG=70` pruning at the end of every wave (`align.c:823–831` / `align.c:1343–1351`), so
in practice `hgh-low` stays small (bounded roughly by a small multiple of `WAVE_LAG`) —
comfortably ≤ ~256 diagonals, i.e. ≤ 8 diagonals per lane, matching the brief's sizing.
Each lane owns its slice of `V/M/T/HA/NA` (the per-diagonal furthest-reach state,
`align.c:391–399` forward / `align.c:925–933` reverse) for its assigned diagonals and runs
the `align.c:494–511`-style slide independently and in parallel across the band.

**The `dif` loop stays sequential.** The outer `while (more && lasta ...)` loop
(`align.c:583` / `align.c:1108`) cannot be parallelized across `dif`: computing diagonal
`k`'s new furthest-reach at wave `dif` reads `V[k-1],V[k],V[k+1]` from wave `dif-1`
(`align.c:687–716` forward, `align.c:1207–1236` reverse) — a strict antidiagonal DP
recurrence. Each `dif` step is one kernel-internal warp-synchronous round; only the
per-diagonal work *within* a `dif` step is parallel across lanes.

**Warp argmax replaces scalar best/trim tracking.** The CPU's `besta`/`trima` (and their
`x`/`diff`/`ha` companions) are single scalars updated sequentially as the `for k = hgh
downto low` loop visits diagonals one at a time (forward, `align.c:680,770`; reverse,
`align.c:1200,1290`). On GPU this becomes a warp-wide argmax/reduction over each lane's
local `c` value once per `dif` step, since all diagonals in the band are evaluated
concurrently rather than in a fixed sequential order.

**`cells[]` become a bounded per-warp checkpoint buffer.** The CPU's `Pebble cells[]`
(`align.c:347–352`) is an unbounded, amortized-`Realloc`-growth array shared across the
whole wave (growth sites: `align.c:473–480,515–521,748–756` forward;
`align.c:1001–1008,1042–1049,1268–1276` reverse). Checkpoint density is bounded by
`tspace=100` (one push per diagonal per 100bp crossed, `align.c:403,482,492,746,767`), and
the band width and total wave count (`dif`) are themselves bounded (by `WAVE_LAG` pruning
and `TRIM_MLAG`, respectively), so the *total* checkpoint count per alignment call is
bounded. The GPU port should size a fixed per-warp checkpoint buffer from these bounds
rather than reproducing the CPU's dynamic-growth strategy, and replace the CPU's pointer-
chasing linked-list walk (`align.c:863–891` / `align.c:1383–1449`) with an equivalent
bounded, index-addressable structure.

### 3 expected GPU deviations (to attribute in later validation)

1. **Parallel band-order tie-breaking.** The CPU visits diagonals in a fixed sequential
   order (`for k = hgh; k >= low; k--` forward `align.c:680`; `for k = low; k <= hgh; k++`
   reverse `align.c:1200`) and only overwrites `besta` on a *strict* `>`/`<` comparison
   (`align.c:770` `if (c > besta)`, `align.c:1290` `if (c < besta)`), so among diagonals
   tied for the furthest reach, the first one visited in that fixed order wins. A GPU
   warp-argmax evaluates all diagonals in the band concurrently; unless its tie-break rule
   is deliberately made to match the CPU's `k`-order priority, ties can resolve to a
   different diagonal, changing which lineage's checkpoint chain (and hence which trace
   points) get reported for degenerate/tied wavefronts.
2. **Argmax vs. scalar first-best.** Related but distinct from (1): the CPU's trim
   acceptance (`align.c:770–783` / `align.c:1290–1303`) is gated on `c > besta` being true
   *at the moment that diagonal is processed*, so a later-processed diagonal with an equal
   `c` never becomes the new best even though it's evaluated in the same wave. A one-shot
   warp reduction naturally picks *a* maximum but has no inherent notion of "processed
   first," so the exact rule used to pick among equal-`c` lanes needs to be an explicit,
   documented convention (e.g. lowest diagonal index wins) if bit-exactness against the CPU
   trim point is required.
3. **`int` vs. native GPU position width.** The CPU uses 32-bit `int` throughout for `V`,
   `M`, `HA`, `NA`, and the `besta/trima/morea` family. Retired-diagonal sentinels are
   *directional*, not symmetric: forward maximizes reach as it extends rightward, so a
   retired diagonal must read as artificially *low* — it uses plain **`-1`**
   (`align.c:657` `V[low] = -1;`, `align.c:664` `V[hgh] = am = -1;`, `align.c:675` `ac =
   V[hgh+1] = V[low-1] = -1;`; note `align.c:667` `am = V[--hgh];` is the *keep* branch, not
   a sentinel). Reverse minimizes reach as it extends leftward, so a retired diagonal must
   read as artificially *high* — it uses **`INT32_MAX`** (`align.c:1181` `V[low] = ap =
   INT32_MAX;`, `align.c:1188` `V[hgh] = INT32_MAX;`, `align.c:1195` `ac = V[hgh+1] =
   V[low-1] = INT32_MAX;`). Trace deltas are separately truncated to `uint16` on write
   (`align.c:884–885,900–901,1412–1417,1426–1428`). A GPU port choosing a different native
   width for position/diagonal arithmetic (e.g. wider accumulators, or unsigned diagonal
   indices) must reproduce these exact sentinel values *and* their directional asymmetry, or
   a genome large enough to approach them will diverge from the CPU's clipping behavior.

---

## Verified but out of scope for this doc

`Local_Alignment` (`align.c:1464–1620`) is documented here only as much as needed to give
the waves their inputs/outputs; its DUB_TRIM-based short-alignment re-centering logic
(`align.c:1535,1551–1576`, `DUB_TRIM=45` at `align.c:170`) and ACOMP-flag trace-reversal
post-processing (`align.c:1578–1601`) are not part of the wave kernel itself and are not
covered here.
