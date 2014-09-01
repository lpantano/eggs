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
            cref, cquery, pref = _check_pos(cols[1], cols[2], cols[3], cols[4])
            print "%s %s %s %s %s %s %s %s" % (cols[1], cref, cquery, pref[0], pref[1], int(cols[3])-1, cols[4], cols[-1])


def _check_pos(name, chr, start, end):
    """compare name of the seq with position map"""
    cref = "-".join(name.split("_")[0].lower().split("-")[:3])
    cquery = "-".join(chr.lower().split("-")[:3])
    logger.debug("%s %s" % (cref, cquery))
    pref = name.split("_")[1].split(":")
    return cref, cquery, pref


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
