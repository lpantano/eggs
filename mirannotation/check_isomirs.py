import logging
from argparse import ArgumentParser

logger = logging.getLogger()


class realSeq:
    def __init__(self, seq):
        self.seq = seq
        self.miraligner = ""
        self.bench = ""


def _read_miraligner(in_file):
    """read miraligner results"""
    with open(in_file) as in_handle:
        in_handle.readline()
        for line in in_handle:
            cols = line.strip().split("\t")
            logger.debug(line)
            ok_mir, ok_add, ok_mut, ok_t5, ok_t3 = _check_pos(cols[1], cols[2], cols[5], cols[6], cols[7], cols[8])
            print "%s %s %s %s %s %s" % (cols[1], ok_mir, ok_mut, ok_add, ok_t5, ok_t3)
            logger.debug("res %s %s %s %s %s %s" % (cols[1], ok_mir, ok_mut, ok_add, ok_t5, ok_t3))


def _check_pos(name, chr, mut, add, t5, t3):
    """compare name of the seq with position map"""
    ok_mir, ok_add, ok_mut, ok_t5, ok_t3 = 0, 0, 0, 0, 0
    cref = name.split("_")[1].lower()
    cquery = chr.lower()
    logger.debug("%s %s" % (cref, cquery))
    t5_ref, t3_ref = name.split("_")[3].split(":")
    mut_ref = name.split("_")[4].split(":")[1]
    add_ref = name.split("_")[5].split(":")[1]
    if (len(t3)-2 == abs(int(t3_ref))) or (t3 == t3_ref):
        ok_t3 = True
    if (len(t5)-2 == abs(int(t5_ref))) or (t5 == t5_ref):
        ok_t5 = True
    if mut != "0" and mut_ref != "null":
        ok_mut = True
    if mut == "0" and mut_ref == "null":
        ok_mut = True
    if len(add)-2 == len(add_ref) or (add == "0" and add_ref == "null"):
        ok_add = True
    if cref == cquery:
        ok_mir = True
    return ok_mir, ok_add, ok_mut, ok_t5, ok_t3


if __name__ == "__main__":
        parser = ArgumentParser(description="check annotation")
        parser.add_argument("--log", help="debug mode", action='store_true')
        args = parser.parse_args()
        if not args.log:
            numeric_level = getattr(logging, "INFO", None)
        else:
            numeric_level = getattr(logging, "DEBUG", None)
        logging.basicConfig(level=numeric_level)
        _read_miraligner("miraligner/sim.20.hsa.mirna")
