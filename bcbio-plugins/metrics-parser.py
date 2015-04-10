"""Parse final project yaml file to extract metrics and output a tsv file"""

import os.path as op
import glob
from argparse import ArgumentParser
import yaml


def _extract(data):
    res = {}
    for l in data['samples']:
        res[l['summary']['metrics']['Name']] = l['summary']['metrics']
    return res


def _header(h):
    fixed = [e.replace(" ", "_").replace("(", "").replace(")", " ").replace(" ", "") for e in h]
    return fixed


def _matrix(res):
    print "\t".join(_header(res[res.keys()[0]].keys()))
    for sample in res:
        print "\t".join(map(str, res[sample].values()))


def _get_files(path):
    return [e for e in glob.glob(path + '/*/qc/fastqc/*txt')]


def _get_quality(fn, section, name, col=2):
    sample = fn.split("/")[-4]
    with open(fn) as in_handle:
        for line in in_handle:
            if line.startswith(section):
                line = in_handle.next()
                line = in_handle.next()
                while not line.startswith('>>'):
                    values = line.split("\t")[0:col]
                    pos = values[0].split("-")
                    if len(pos) > 1:
                        for p in range(int(pos[0]), int(pos[1])):
                            print "\t".join([name, sample, str(p)] + values)
                    else:
                        print "\t".join([name, sample, pos[0]] + values)

                    line = in_handle.next().strip()


if __name__ == "__main__":
        parser = ArgumentParser(description="Run a single cell analysis.")
        parser.add_argument("--yaml", help="yaml file")
        parser.add_argument("--fastqc_read", action="store_true", help="join all fastqc results")
        parser.add_argument("--fastqc_content", action="store_true", help="join all fastqc results")
        parser.add_argument("--dir", help="final bcbio dir")
        args = parser.parse_args()

        if args.yaml:
            data = yaml.load(open(args.yaml, 'r'))
            metrics = _extract(data)
            _matrix(metrics)
        if args.fastqc_read:
            fastqc_output = _get_files(op.abspath(args.dir))
            [_get_quality(fn, '>>Per base sequence quality', 'nt_qual') for fn in fastqc_output]
            [_get_quality(fn, '>>Sequence Length Distribution', 'length') for fn in fastqc_output]
        if args.fastqc_content:
            fastqc_output = _get_files(op.abspath(args.dir))
            [_get_quality(fn, '>>Per base sequence content', 'content', 5) for fn in fastqc_output]
