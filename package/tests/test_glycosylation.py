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


@pytest.fixture(scope='session')
def glycosylation_fixture(cglycan):
    return sails.glycosylate(*cglycan)


def test_structure_integrity(glycosylation_fixture):
    s, m, l = glycosylation_fixture
    assert s
    assert isinstance(s, gemmi.Structure)


def test_mtz_integrity(glycosylation_fixture):
    s, m, l = glycosylation_fixture
    assert m
    assert isinstance(m, gemmi.Mtz)


def test_output_integrity(glycosylation_fixture):
    s, m, l = glycosylation_fixture
    assert l
    assert 'cycles' in l
    assert 'date' in l
    cycles = l['cycles']
    assert cycles
    cycle = cycles[0]
    assert 'cycle' in cycle
    assert isinstance(cycle['cycle'], int)


def test_entries_completeness(glycosylation_fixture):
    s, m, l = glycosylation_fixture
    cycle = l['cycles'][0]
    assert 'entries' in cycle
    entries = cycle['entries']
    assert 'D-MAN-1' in entries


def test_entries_integrity(glycosylation_fixture):
    s, m, l = glycosylation_fixture
    entries = l['cycles'][0]['entries']
    sugar = entries['D-MAN-1']
    assert 'rscc' in sugar
    assert 'rsr' in sugar
    assert 'dds' in sugar


def test_sugar_scores(glycosylation_fixture):
    s, m, l = glycosylation_fixture
    sugar = l['cycles'][0]['entries']['D-MAN-1']
    assert sugar['rscc'] > 0.7
    assert sugar['rsr'] > 0.9
    assert sugar['dds'] < 0.75
