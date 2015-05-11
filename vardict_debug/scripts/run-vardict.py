"""
calculate coverage across a list of regions
"""
import os
import os.path as op
import yaml

from argparse import ArgumentParser
from ichwrapper import cluster, arguments

from bcbio.install import _get_data_dir
from bcbio.utils import safe_makedir

from vardict.run import trymemory


def _update_algorithm(data, resources):
    """
    Update algorithm dict with new cores set
    """
    new_data = []
    for sample in data:
        sample[0]['config']['algorithm'] = resources
        new_data.append(sample)
    return new_data

def _prepare_samples(args):
    """
    create dict for each sample having all information
    """
    if args.galaxy:
        system_config = args.galaxy
    else:
        system_config = os.path.join(_get_data_dir(), "galaxy", "bcbio_system.yaml")
    config = yaml.load(open(system_config))
    config['algorithm'] = {}
    data = []
    for sample in args.files:
        dt = {}
        dt['name'] = splitext_plus(op.basename(sample))[0]
        dt['config'] = config
        dt['bed'] = op.abspath(sample)
        data.append([dt])
    return data

def detect(args):
    assert args.files, "Need a set of fastq files"

    for mem in [6, 8, 10, 12, 14, 16, 18, 20, 24]:
        safe_makedir("mem%s" % mem)
        resources = {'name': "mem%s" % mem, 'mem': mem, 'cores': 1}
        data = _prepare_samples(args)
        data = _update_algorithm(data, resources)
        data = cluster.send_job(trymemory, data, args, resources)

if __name__ == "__main__":
    parser = ArgumentParser(description="task related to allele methylation specific")
    parser = arguments.myargs(parser)

    parser.add_argument("files", nargs="*", help="Bam files.")
    args = parser.parse_args()

    detect(args)
