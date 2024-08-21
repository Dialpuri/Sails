from pathlib import Path
from pprint import pprint
import xml.etree.ElementTree as ET

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
    return s, m, 1, "FP", "SIGFP", "", "", sails.Type.c_glycosylate


def test_xtal_cglycosylation(cglycan):
    s, m, l, snfgs = sails.glycosylate_xtal(*cglycan)
    assert s
    assert m
    assert l
    assert snfgs

    assert isinstance(s, gemmi.Structure)
    assert isinstance(m, gemmi.Mtz)

    # test log file
    assert 'cycles' in l
    assert 'date' in l

    cycles = l['cycles']
    assert cycles
    cycle = cycles[0]
    assert 'cycle' in cycle
    assert isinstance(cycle['cycle'], int)
    assert 'entries' in cycle
    entries = cycle['entries']

    expected_key = 'D-AMAN-1'
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

    # test snfg output
    assert 1 in snfgs
    c1 = snfgs[1]
    expected_trp_key = "C-TRP-16"
    assert expected_trp_key in c1

    root = ET.fromstring(c1[expected_trp_key])
    assert root.tag == '{http://www.w3.org/2000/svg}svg'
    assert root.attrib['width'] == '1200'
    assert root.attrib['height'] == '800'

    line_count = len(root.findall('{http://www.w3.org/2000/svg}line'))
    rect_count = len(root.findall('{http://www.w3.org/2000/svg}rect'))
    circle_count = len(root.findall('{http://www.w3.org/2000/svg}circle'))
    text_count = len(root.findall('{http://www.w3.org/2000/svg}text'))
    tspan_count = len(root.findall('.//{http://www.w3.org/2000/svg}tspan'))

    # Assert the expected counts
    assert line_count == 3, f"Expected 3 <line> elements, but found {line_count}"
    assert rect_count == 1, f"Expected 1 <rect> element, but found {rect_count}"
    assert circle_count == 1, f"Expected 1 <circle> element, but found {circle_count}"
    assert text_count == 3, f"Expected 3 <text> elements, but found {text_count}"
    assert tspan_count == 1, f"Expected 1 <tspan> element, but found {tspan_count}"


