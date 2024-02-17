import pysails as sails
import pytest
def test_sails_input_exists():
    assert "Input" in dir(sails)

def test_sails_output_exists():
    assert "Output" in dir(sails)

def test_sails_find_exists():
    assert "Find" in dir(sails)

def test_sails_input_constructors():
    with pytest.raises(TypeError):
        i = sails.Input("", "", "", "", "", "", "", "", "")
    i = sails.Input("", "", "", "", "", "", "", 0, "")

def test_sails_output_constructors():
    o = sails.Output("")