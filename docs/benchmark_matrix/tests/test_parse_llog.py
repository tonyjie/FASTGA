import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from aggregate_matrix import sec, res, parse_llog

FIX = os.path.join(HERE, "fixtures", "human", "logs")
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

def test_aggregate_one_point(tmp_path):
    from aggregate_matrix import aggregate
    d = tmp_path / "divergence" / "human" / "logs"
    d.mkdir(parents=True)
    for r in (1, 2, 3):
        (d / f"rep{r}.Llog").write_text(open(os.path.join(FIX, f"rep{r}.Llog")).read())
    tsv = aggregate(str(tmp_path / "divergence"), [("human", 0)])
    body = open(tsv).read()
    assert "human" in body and "Sort+align" in body
    assert os.path.exists(os.path.join(str(tmp_path / "divergence"), "divergence_phase_share.png"))
