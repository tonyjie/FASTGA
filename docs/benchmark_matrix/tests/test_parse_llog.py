import os, re, sys, statistics
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from aggregate_matrix import sec, res, parse_llog, peak_rss_mb

FIX = os.path.join(HERE, "fixtures", "human", "logs")
FIX_TIME = os.path.join(HERE, "fixtures", "human", "time")
SYN = os.path.join(HERE, "fixtures", "synthetic")

def test_sec_formats():
    assert sec("6.706") == 6.706
    assert abs(sec("2:06.205") - 126.205) < 1e-6
    assert abs(sec("8:14.618") - 494.618) < 1e-6

def test_res_line():
    assert res("  Resources for phase:  40:05.580u  9:19.434s  8:14.618w  599.5%") == \
        (2405.58, 559.434, 494.618)
    assert res("no resources here") is None

def test_parse_human_shares():
    p = parse_llog(os.path.join(FIX, "rep1.Llog"))
    total = p["GDB"] + p["GIX"] + p["Seed merge"] + p["Sort+align"]
    # Known human breakdown: GDB ~3%, GIX ~15%, merge ~1%, sort+align ~81%
    assert 0.78 <= p["Sort+align"] / total <= 0.84
    assert p["Seed merge"] / total < 0.03
    assert 0.10 <= p["GIX"] / total <= 0.20
    assert p["n_aln"] == 518037
    assert p["rss_mb"] == 19.0

def test_intervening_fastga_not_counted():
    # Two FastGA invocations share a single prep cycle (FAtoGDB/GIXmake x2).
    # The FIRST (non-final) FastGA invocation's own "Total Resources" line
    # (520.000w, 15MB) must NOT be summed into GDB/GIX or reported as rss_mb;
    # only the SECOND (final) invocation's phases/rss/aln stats should surface.
    p = parse_llog(os.path.join(SYN, "two_fastga_shared_prep.Llog"))
    assert p["GDB"] == 18.0
    assert p["GIX"] == 87.0
    assert p["Seed merge"] == 7.0
    assert p["Sort+align"] == 480.0
    assert p["rss_mb"] == 19.0
    assert p["n_aln"] == 777
    assert p["ave_len"] == 300

def test_peak_rss_from_time():
    # Expected: median across whatever rep*.time fixtures exist, computed from
    # the actual file contents (not hardcoded), of "Maximum resident set size
    # (kbytes): N" converted kb -> MB.
    kb_vals = []
    for fn in sorted(os.listdir(FIX_TIME)):
        text = open(os.path.join(FIX_TIME, fn)).read()
        m = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', text)
        assert m, f"no RSS line in {fn}"
        kb_vals.append(int(m.group(1)))
    expected_mb = statistics.median(kb_vals) / 1024
    got = peak_rss_mb(FIX_TIME)
    assert got is not None
    assert abs(got - expected_mb) < 50

def test_peak_rss_missing_dir_returns_none():
    assert peak_rss_mb(os.path.join(HERE, "fixtures", "no_such_dir")) is None

def test_aggregate_one_point(tmp_path):
    from aggregate_matrix import aggregate
    d = tmp_path / "divergence" / "human" / "logs"
    d.mkdir(parents=True)
    for r in (1, 2, 3):
        (d / f"rep{r}.Llog").write_text(open(os.path.join(FIX, f"rep{r}.Llog")).read())
    t = tmp_path / "divergence" / "human" / "time"
    t.mkdir(parents=True)
    (t / "rep1.time").write_text(open(os.path.join(FIX_TIME, "rep1.time")).read())
    expected_mb = peak_rss_mb(str(t))
    tsv = aggregate(str(tmp_path / "divergence"), [("human", 0)])
    body = open(tsv).read()
    assert "human" in body and "Sort+align" in body
    assert os.path.exists(os.path.join(str(tmp_path / "divergence"), "divergence_phase_share.png"))
    header = body.splitlines()[0]
    assert "cpu_sortalign" in header
    line = [ln for ln in body.splitlines() if ln.startswith("human\t")][0]
    cells = line.split("\t")
    rss_cell = cells[-7]  # ...rss_mb, n_aln, ave_len, cpu_GDB, cpu_GIX, cpu_seed, cpu_sortalign
    assert rss_cell == f"{expected_mb:.0f}"
    assert rss_cell != "19"
    cpu_sortalign_cell = cells[-1]
    assert cpu_sortalign_cell != ""
    assert 400 <= float(cpu_sortalign_cell) <= 800
