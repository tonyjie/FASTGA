import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from aggregate_matrix import sec, res, parse_llog

FIX = os.path.join(HERE, "fixtures", "human", "logs")

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
