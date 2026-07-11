# Opt C scan-once + tmpfs ktab (Approach B) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make FastGA's bilateral chunked build/merge (Opt C) scan each genome once instead of once per chunk, run its scratch on tmpfs so the ktab never hits real disk, and measure perf/memory/storage vs. baseline and Opt C on EXAMPLE + human — bit-exact.

**Architecture:** Add a GIXmake `-R` (reuse) mode: the initial `-n` scan build persists its pos-lists (root-named, not unlinked) plus a tiny `Buckets` counts sidecar; per-chunk sort builds pass `-R` to skip `distribute()` and reuse that state. FastGA's chunk loop gains `-R` on its per-chunk calls and deletes the persisted scratch afterward. tmpfs is a run-time choice (`-P /dev/shm`), not code.

**Tech Stack:** C (GIXmake.c, FastGA.c), pthreads; verification via `ONEview` + `md5sum` and `cmp` on `.ktab` files; bash measurement harness under `docs/agent_optimization/`.

## Global Constraints

- Behavior of the **non-chunked path must not change** — every new behavior is gated on `-n`/`-R`/`-C` (the cooperative-chunk signals FastGA already emits only in chunked mode).
- **Bit-exact** is the gating test: `ONEview <out.1aln> | tail -n +2 | md5sum` must equal the `ddeea32` baseline at matched thread count. EXAMPLE = **323,569** alignments; human = **518,037** alignments.
- No commits pushed without the user's OK. Work on branch `fused-scan-once` (already created off `agent-optimization`).
- Counts sidecar path: `<SORT_PATH>/.<ROOT>.gcnt`; persistent pos-lists: `<SORT_PATH>/.post.<ROOT>.<k>.idx` (k = 0 … NPARTS·NTHREADS−1).
- `Buckets` is one contiguous block: `Buckets[0] = Malloc(NTHREADS*NUM_BUCK*sizeof(int64))`, `Buckets[p]=Buckets[p-1]+NUM_BUCK` (`GIXmake.c:2092-2094`) — so `Buckets[0]` addresses all `NTHREADS*NUM_BUCK` int64s.

---

## File structure

- `GIXmake.c` — all index-build changes: new `REUSE` global + `-R` flag, root-named/persistent pos-lists, `write_counts_sidecar`/`read_counts_sidecar`, the `-R` skip-`distribute` branch in `make_index`/`main`.
- `FastGA.c` — chunk-loop only: add `-R` to the per-chunk `GIXmake` commands (`~5365-5391`); delete persisted pos-lists + `.gcnt` after the loop (`~5436`).
- `docs/agent_optimization/fused_stages/` (new) — measurement harness: `run_fused.sh`, `analyze.py`, `stage_data/`, `README.md`.
- `docs/agent_optimization/plan_improvement.md` — fold in measured numbers at the end.

---

### Task 1: GIXmake — `-R` flag, persistent root-named pos-lists, counts sidecar (write side)

**Files:**
- Modify: `GIXmake.c` — globals (~`66-68`), flag parse (`1745`, `1794`, `~1736`), `POST_NAME` (`1912`), pos-list open/unlink (`2105-2111`), add `write_counts_sidecar` near `651`, call it after `distribute()` (`~2127`).

**Interfaces:**
- Produces: global `static int REUSE;`; `static void write_counts_sidecar(char *tpath, char *troot);` writing `<tpath>/.<troot>.gcnt` = `int NTHREADS`, `int NPARTS`, then `NTHREADS*NUM_BUCK` `int64` from `Buckets[0]`. Persistent pos-list name `<SORT_PATH>/.post.<ROOT>.<k>.idx` when `STUB_ONLY||REUSE`.

- [ ] **Step 1: Add the `REUSE` global.** Next to `STUB_ONLY` (`GIXmake.c:67`):

```c
static int STUB_ONLY;    //  -n: write .gix stub only, no ktab partition files
static int REUSE;        //  -R: reuse persisted pos-lists (+.gcnt) from a prior -n scan, skip distribute()
```

- [ ] **Step 2: Parse `-R`.** In the flag block, change `ARG_FLAGS("vn")` (`GIXmake.c:1745`) to `ARG_FLAGS("vnR")`; init `REUSE = 0;` next to `STUB_ONLY = 0;` (`~1736`); set it next to `STUB_ONLY = flags['n'];` (`1794`):

```c
    STUB_ONLY = flags['n'];
    REUSE     = flags['R'];
```

- [ ] **Step 3: Root-name pos-lists in cooperative mode.** Replace the `POST_NAME` assignment (`GIXmake.c:1912`):

```c
    if (STUB_ONLY || REUSE)
      POST_NAME = Strdup(Catenate(SORT_PATH,"/.post.",TROOT,"."),"Allocating post name");
    else
      POST_NAME = Strdup(Catenate(SORT_PATH,"/.",Numbered_Suffix("post.",getpid(),"."),""),
                         "Allocating post name");
```

- [ ] **Step 4: Persist pos-lists (don't unlink) on the scan build.** Replace the create/unlink loop body (`GIXmake.c:2105-2111`):

```c
        { name = Numbered_Suffix(POST_NAME,k,".idx");
          Units[k] = open(name,O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
          if (Units[k] < 0)
            { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,name);
              exit (1);
            }
          if (!(STUB_ONLY || REUSE))   //  cooperative mode keeps pos-lists for the -R builds
            unlink(name);
        }
```

- [ ] **Step 5: Add `write_counts_sidecar`.** Insert immediately after `write_ksplit_sidecar` (`GIXmake.c:668`):

```c
//  Buckets counts sidecar (.gcnt): the per-(thread,bucket) entry counts that scan_thread
//  leaves in Buckets after distribute() — exactly what k_sort's finger step reads. Written by
//  the -n scan build right after distribute() (before k_sort mangles Buckets into fingers),
//  read by -R builds so they can k_sort without re-scanning.  Format: int NTHREADS, int NPARTS,
//  then NTHREADS*NUM_BUCK int64 from the contiguous Buckets[0].

static void write_counts_sidecar(char *tpath, char *troot)
{ char *sname;
  int   fd;

  sname = Malloc(strlen(tpath)+strlen(troot)+20,"Allocating gcnt path");
  if (sname == NULL) return;
  sprintf(sname,"%s/.%s.gcnt",tpath,troot);
  fd = open(sname,O_WRONLY|O_CREAT|O_TRUNC,0666);
  if (fd < 0) { free(sname); return; }
  write(fd,&NTHREADS,sizeof(int));
  write(fd,&NPARTS,sizeof(int));
  write(fd,Buckets[0],((int64) NTHREADS)*NUM_BUCK*sizeof(int64));
  close(fd);
  free(sname);
}
```

- [ ] **Step 6: Call it after `distribute()` on the scan build.** After `distribute(gdb);` (`GIXmake.c:2126`), before the `.split` write:

```c
  distribute(gdb);

  if (STUB_ONLY)               //  the scan build: snapshot Buckets before k_sort mangles it
    write_counts_sidecar(TPATH,TROOT);

  if (!CHUNKED)
    write_ksplit_sidecar(TPATH,TROOT);
```

- [ ] **Step 7: Build.**

Run: `make GIXmake 2>&1 | tail -5`
Expected: compiles cleanly, `GIXmake` binary produced.

- [ ] **Step 8: Verify persistence + sidecar (manual scan build).**

```bash
cd /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/EXAMPLE
../FAtoGDB HAP1.fasta.gz 2>/dev/null; rm -f .HAP1.post.* .HAP1.gcnt
../GIXmake -n -T4 -P. HAP1
ls -la .HAP1.gcnt .HAP1.post.*.idx | head
```
Expected: `.HAP1.gcnt` exists (size = 8 + NTHREADS·1024·8 bytes) and multiple `.HAP1.post.<k>.idx` files persist (not unlinked). Clean up: `rm -f .HAP1.post.* .HAP1.gcnt .HAP1.gix .HAP1.split`.

- [ ] **Step 9: Commit.**

```bash
git add GIXmake.c
git commit -m "GIXmake: persist root-named pos-lists + write .gcnt counts sidecar on -n scan"
```

---

### Task 2: GIXmake — `-R` reuse path (skip `distribute`, reconstruct state, sort)

**Files:**
- Modify: `GIXmake.c` — add `read_counts_sidecar` near `690`; the `-R` branch replacing `distribute()` (`~2126`); pos-list open already handled in Task 1 Step 4 (for `-R`, `O_CREAT` on an existing file is harmless but confirm read path).

**Interfaces:**
- Consumes: `.gcnt` and `.post.<ROOT>.*.idx` from Task 1; `REF_KSPLIT`/`REF_NPARTS`/`NPARTS` from `-X` (`GIXmake.c:2016`).
- Produces: `static int read_counts_sidecar(char *tpath, char *troot);` filling `Buckets[0]`; a populated `Ksplit[]` + `Select[]` when `REUSE`.

- [ ] **Step 1: Add `read_counts_sidecar`.** After `write_counts_sidecar`:

```c
static int read_counts_sidecar(char *tpath, char *troot)
{ char *sname;
  int   fd, nt, np;

  sname = Malloc(strlen(tpath)+strlen(troot)+20,"Allocating gcnt path");
  if (sname == NULL) return (-1);
  sprintf(sname,"%s/.%s.gcnt",tpath,troot);
  fd = open(sname,O_RDONLY);
  if (fd < 0)
    { fprintf(stderr,"%s: -R cannot open counts sidecar %s\n",Prog_Name,sname);
      free(sname); return (-1);
    }
  if (read(fd,&nt,sizeof(int)) != sizeof(int) || read(fd,&np,sizeof(int)) != sizeof(int)
      || nt != NTHREADS || np != NPARTS)
    { fprintf(stderr,"%s: -R counts sidecar %s mismatch (T/NPARTS)\n",Prog_Name,sname);
      close(fd); free(sname); return (-1);
    }
  if (read(fd,Buckets[0],((int64) NTHREADS)*NUM_BUCK*sizeof(int64))
        != (ssize_t)(((int64) NTHREADS)*NUM_BUCK*sizeof(int64)))
    { fprintf(stderr,"%s: -R counts sidecar %s truncated\n",Prog_Name,sname);
      close(fd); free(sname); return (-1);
    }
  close(fd); free(sname);
  return (0);
}
```

- [ ] **Step 2: Branch on `REUSE` instead of scanning.** Replace the `distribute(gdb); … write_ksplit_sidecar` block (`GIXmake.c:2126-2134`) with:

```c
  if (REUSE)
    { int i, p;                    //  reconstruct scan output without re-scanning
      for (i = 0; i <= NPARTS; i++)   //  Ksplit from -X (REF_KSPLIT already read)
        Ksplit[i] = REF_KSPLIT[i];
      p = 0;                          //  Select from Ksplit (mirrors GIXmake.c:733-737)
      for (i = 0; i < NUM_BUCK; i++)
        { while (p < NPARTS && Ksplit[p+1] <= i) p += 1;
          Select[i] = p;
        }
      if (read_counts_sidecar(TPATH,TROOT) < 0)   //  Buckets counts from .gcnt
        exit (1);
    }
  else
    { distribute(gdb);

      if (STUB_ONLY)
        write_counts_sidecar(TPATH,TROOT);

      if (!CHUNKED)
        write_ksplit_sidecar(TPATH,TROOT);
    }
```

(This subsumes the Task 1 Step 6 edit — the `else` branch is that code.) `-R` requires `-X`, so `REF_KSPLIT != NULL`; if a caller omits `-X`, `REF_KSPLIT` is NULL — add a guard: right after arg parse where `REUSE`/`REF_STUB` are known, `if (REUSE && REF_STUB == NULL) { fprintf(stderr,"%s: -R requires -X\n",Prog_Name); exit(1); }`.

- [ ] **Step 3: Build.**

Run: `make GIXmake 2>&1 | tail -5`
Expected: clean compile.

- [ ] **Step 4: Isolated bit-exact test — `-R` ktab == normal ktab.**

```bash
cd /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/EXAMPLE
../FAtoGDB HAP1.fasta.gz 2>/dev/null
# reference: normal whole build
rm -f .HAP1.ktab.* .HAP1.gix .HAP1.split .HAP1.gcnt .HAP1.post.*
../GIXmake -T4 -P. HAP1 && mkdir -p /tmp/ref && cp .HAP1.ktab.* /tmp/ref/
# scan-once: -n scan, then -R per-partition (NPARTS from .split)
rm -f .HAP1.ktab.* .HAP1.gix
../GIXmake -n -T4 -P. HAP1                       # persists pos-lists + .gcnt + .split
NP=$(od -An -i -N4 .HAP1.split | tr -d ' ')
for p in $(seq 1 $NP); do ../GIXmake -R -C$p:$p -X./.HAP1.split -T4 -P. HAP1; done
# compare: chunked build numbers ktab.1 per chunk -> rebuild names, cmp against ref
for p in $(seq 1 $NP); do cmp /tmp/ref/.HAP1.ktab.$p .HAP1.ktab.1 && echo "part $p OK"; rm -f .HAP1.ktab.1; ../GIXmake -R -C$p:$p -X./.HAP1.split -T4 -P. HAP1; done
```
Expected: every partition's `-R` ktab is byte-identical to the reference. (Note: chunked builds renumber ktab to `ktab.1`; the loop rebuilds one part at a time and `cmp`s.) If any differ, the reconstructed `Ksplit`/`Select`/`Buckets` diverged — stop and diff.

Clean up: `rm -rf /tmp/ref .HAP1.ktab.* .HAP1.gix .HAP1.split .HAP1.gcnt .HAP1.post.*`.

- [ ] **Step 5: Commit.**

```bash
git add GIXmake.c
git commit -m "GIXmake: add -R reuse mode (skip distribute, rebuild Ksplit/Select/Buckets from sidecars)"
```

---

### Task 3: FastGA — thread `-R` through the chunk loop + clean up persisted scratch

**Files:**
- Modify: `FastGA.c` — the four per-chunk `GIXmake` command `sprintf`s (`~5365-5391`); add pos-list + `.gcnt` deletion after the chunk loop (`~5436`).

**Interfaces:**
- Consumes: Task 1/2 GIXmake `-R`.
- Produces: chunked FastGA runs that scan each genome once.

- [ ] **Step 1: Add `-R` to the per-chunk commands.** In each of the four `sprintf(cmd,"GIXmake%s -C%d:%d -X%s/.%s.split -T%d -P%s %s/%s", …)` sites (g1 LOG/no-LOG, g2 LOG/no-LOG; `FastGA.c:5365,5370,5384,5389`), insert ` -R` after `GIXmake%s`:

```c
              sprintf(cmd,"GIXmake%s -R -C%d:%d -X%s/.%s.split -T%d -P%s %s/%s",
                      VERBOSE?" -v":"",cfirst,clast,PATH1,ROOT1,
                      NTHREADS,SORT_PATH,PATH1,ROOT1);
```
(Same insertion in all four; PATH/ROOT differ per genome as they already do.)

- [ ] **Step 2: Delete persisted pos-lists + `.gcnt` after the chunk loop.** After the loop closes (`FastGA.c:5436`, before/near the stub deletions):

```c
        //  Scan-once cleanup: the -n scan builds persisted pos-lists + .gcnt for the -R
        //  chunk builds to reuse; remove them now.
        { int p; char *gn = Malloc(strlen(SORT_PATH)+strlen(ROOT1)+strlen(ROOT2)+40,"cleanup");
          for (p = 0; p < g2_nparts*NTHREADS; p++)
            { sprintf(gn,"%s/.post.%s.%d.idx",SORT_PATH,ROOT1,p); unlink(gn);
              sprintf(gn,"%s/.post.%s.%d.idx",SORT_PATH,ROOT2,p); unlink(gn); }
          sprintf(gn,"%s/.%s.gcnt",SORT_PATH,ROOT1); unlink(gn);
          sprintf(gn,"%s/.%s.gcnt",SORT_PATH,ROOT2); unlink(gn);
          free(gn);
        }
```
(`g2_nparts` is already in scope, `= read_ksplit_nparts(PATH1,ROOT1)` at `FastGA.c:5328`; pos-lists live in `SORT_PATH` = the `-P` dir.)

- [ ] **Step 3: Build all.**

Run: `make FastGA GIXmake 2>&1 | tail -5`
Expected: clean compile.

- [ ] **Step 4: Bit-exact end-to-end on EXAMPLE (chunked).**

```bash
cd /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/EXAMPLE
# baseline reference md5 (build stock ddeea32 once if not cached), T=4:
#   ../FastGA -T4 -1:/tmp/base.1aln HAP1 HAP2   (with a ddeea32 binary)  -> BASE_MD5
../FastGA -T4 -C4 -P/dev/shm -1:/tmp/fused.1aln HAP1 HAP2
../ONEview /tmp/fused.1aln | tail -n +2 | md5sum
../ALNshow -w /tmp/fused.1aln HAP1 HAP2 2>/dev/null | grep -c '^' || true
```
Expected: alignment payload md5 equals the `ddeea32` baseline at T=4, and **323,569** non-redundant alignments. Also confirm no leftover `.post.*`/`.gcnt` in `/dev/shm` after exit.

- [ ] **Step 5: Commit.**

```bash
git add FastGA.c
git commit -m "FastGA: pass -R to per-chunk GIXmake (scan once) and clean up persisted pos-lists/.gcnt"
```

---

### Task 4: Measurement harness + EXAMPLE results (baseline / Opt C / fused-B)

**Files:**
- Create: `docs/agent_optimization/fused_stages/run_fused.sh`, `analyze.py`, `README.md`, `stage_data/`.

**Interfaces:**
- Consumes: the `fused-scan-once` FastGA/GIXmake; a `ddeea32` baseline binary and the current `agent-optimization` (Opt C) binary for comparison.
- Produces: per-config wall + 3-phase breakdown, real-disk/tmpfs/RSS footprint, and the bit-exact md5, on EXAMPLE at T=8 and T=32.

- [ ] **Step 1: Write `run_fused.sh`.** Model it on `docs/agent_optimization/human_stages/run_stages_human.sh` (worktree-build per config, `-L` log, `/usr/bin/time -v`, storage sampling via `du` + `/proc/PID/fd` `st_size`). Three configs: `baseline` (ddeea32, real disk), `optC` (`agent-optimization`, `-C4`/`-C8`, real disk), `fusedB` (`fused-scan-once`, `-C4`/`-C8`, `-P/dev/shm`). Also sample `/dev/shm` occupancy for `fusedB`. Emit `stage_data/<cfg>/{run.Llog,run.time,timeline.tsv}` and the ONEview md5.

- [ ] **Step 2: Run EXAMPLE T=8 and T=32.**

Run: `bash docs/agent_optimization/fused_stages/run_fused.sh EXAMPLE 8 && bash docs/agent_optimization/fused_stages/run_fused.sh EXAMPLE 32`
Expected: all three configs complete; `fusedB` md5 == `baseline` md5 (323,569 aln) at each T.

- [ ] **Step 3: Write `analyze.py` → table + figures.** 3-phase breakdown (GDB / Index+merge / Sort+align) per config; a footprint panel (real-disk vs tmpfs vs RSS). Model on `human_stages/analyze.py`.

Run: `python3 docs/agent_optimization/fused_stages/analyze.py EXAMPLE`
Expected: prints a markdown table; key check — `fusedB` Index+merge < `optC` Index+merge (re-scan removed) and ≈ baseline; `fusedB` real-disk ≈ 0.

- [ ] **Step 4: Commit.**

```bash
git add docs/agent_optimization/fused_stages
git commit -m "measure: fused-B (scan-once + tmpfs) vs baseline/Opt C on EXAMPLE (T=8,32)"
```

---

### Task 5: Human results (T=32)

**Files:**
- Modify: `docs/agent_optimization/fused_stages/README.md` + `stage_data/human/`.

**Interfaces:**
- Consumes: Task 4 harness; the working `5671357` FAtoGDB (ddeea32's segfaults on CHM13 — see `human_stages/`).
- Produces: human three-way results + bit-exact md5.

- [ ] **Step 1: Confirm tmpfs headroom.** `df -h /dev/shm` — need ≥ ~40 GB free for the human `fusedB` run. If short, note it and fall back to a large real-disk `-P` for `fusedB` (still validates scan-once; loses the tmpfs-zero-disk axis) — record which.

- [ ] **Step 2: Run human T=32 (all three configs).**

Run: `bash docs/agent_optimization/fused_stages/run_fused.sh human 32`  (~30–60 min)
Expected: `fusedB` md5 == baseline (**518,037** aln); `fusedB` Index+merge drops from Opt C's ~305 s toward baseline ~102 s.

- [ ] **Step 3: Analyze + record.**

Run: `python3 docs/agent_optimization/fused_stages/analyze.py human`
Expected: table with wall, 3-phase, real-disk peak, tmpfs occupancy, peak RSS, total-RAM, md5 per config.

- [ ] **Step 4: Commit.**

```bash
git add docs/agent_optimization/fused_stages
git commit -m "measure: fused-B human T=32 three-way results (bit-exact)"
```

---

### Task 6: Fold measured numbers into the docs

**Files:**
- Modify: `docs/agent_optimization/plan_improvement.md` (item 1 "Expected outcome" → "Measured"); `docs/agent_optimization/fused_stages/README.md` (summary + reproduce).

- [ ] **Step 1: Replace the item-1 estimate table** in `plan_improvement.md` with the measured EXAMPLE + human numbers, and note the lesser-variant (this is exactly the lesser variant realized — scan-once with ktab briefly on tmpfs, not the full in-RAM merge).

- [ ] **Step 2: Write `fused_stages/README.md`** — what was built (Approach B), the three-way tables, the correctness statement, caveats (tmpfs headroom; ktab still written, to RAM), and copy-paste reproduce.

- [ ] **Step 3: Commit.**

```bash
git add docs/agent_optimization/plan_improvement.md docs/agent_optimization/fused_stages/README.md
git commit -m "docs: record measured fused-B (scan-once + tmpfs) results; update plan_improvement item 1"
```

---

## Self-review

**Spec coverage:** GIXmake persistent pos-lists (Task 1) ✓; counts sidecar write (Task 1) ✓ + read (Task 2) ✓; `-R` skip-distribute + Ksplit/Select rebuild (Task 2) ✓; FastGA `-R` wiring + cleanup (Task 3) ✓; tmpfs run mode (Tasks 3–5, `-P /dev/shm`) ✓; bit-exact verification on both datasets (Tasks 3,4,5) ✓; three-way measurement (Tasks 4,5) ✓; docs (Task 6) ✓. Non-goals (no in-process merge refactor, seed temp stays) respected — no task touches `adaptamer_merge`'s data source.

**Placeholder scan:** the harness scripts (Task 4 Step 1, Task 5) are described by "model on `human_stages/…`" rather than inlined verbatim — acceptable because a working, line-for-line analog exists in-repo to copy; every GIXmake/FastGA code change is shown in full.

**Type consistency:** `REUSE` (int), `write_counts_sidecar(char*,char*)` / `read_counts_sidecar(char*,char*)` used identically across Tasks 1–2; `Buckets[0]` block size `((int64)NTHREADS)*NUM_BUCK*sizeof(int64)` identical in write/read; pos-list name `.post.<ROOT>.<k>.idx` and `.gcnt` name identical in GIXmake (create) and FastGA (cleanup).
