from pathlib import Path
import gemmi
import sails
import pytest


@pytest.fixture(scope='session')
def data_base_path():
    return Path(__file__).parent / "test_data"


@pytest.fixture(scope='session')
def cglycan(data_base_path):
    s_path = data_base_path / "6PLH_deglycosylated.cif"
    m_path = data_base_path / "6PLH.mtz"

    s = gemmi.read_structure(str(s_path))
    m = gemmi.read_mtz_file(str(m_path))
    return s, m, 1, "FP", "SIGFP", sails.c_glycosylate_from_objects


def test_cglycosylation(cglycan):
    s, m, l = sails.glycosylate(*cglycan)
    assert s
    assert m
    assert l
    assert isinstance(s, gemmi.Structure)
    assert isinstance(m, gemmi.Mtz)

    assert 'cycles' in l
    assert 'date' in l

    cycles = l['cycles']
    assert cycles
    cycle = cycles[0]
    assert 'cycle' in cycle
    assert isinstance(cycle['cycle'], int)
    assert 'entries' in cycle
    entries = cycle['entries']

    expected_key = 'D-MAN-1'
    assert expected_key in entries
    assert len(entries.keys()) == 1
    sugar = entries[expected_key]

    rscc_key = 'rscc'
    rsr_key = 'rsr'
    dds_key = 'dds'

    assert rscc_key in sugar
    assert rsr_key in sugar
    assert dds_key in sugar

    rscc_score = sugar[rscc_key]
    rsr_score = sugar[rsr_key]
    dds_score = sugar[dds_key]

    assert rscc_score > 0.7
    assert rsr_score > 0.9
    assert dds_score < 0.75



