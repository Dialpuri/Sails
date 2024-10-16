from sails import get_conformations


def run():
    get_conformations("NAG", "NAG", 4, 1, "src/sails/data", 6)
