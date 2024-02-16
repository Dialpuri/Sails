import pysails as sails

def glycosylate():
    input_data = ["tests/data/5fji_dg.pdb", "tests/data/5fji.mtz", "FP,SIGFP", "", "", "FWT,PHWT", "FREE", 1, ""]
    output = "tests/data/5fji_glycosylated.pdb"
    ip = sails.Input(*input_data)
    op = sails.Output(output)
    f = sails.Find(ip, op)


if __name__ == "__main__":
    glycosylate()