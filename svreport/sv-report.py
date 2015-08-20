"""
It needs as input the final folder of bcbio
python sv-report.py --run report /path/to/final/bcbio

12% chances to work :)
"""
import os
import os.path as op
import yaml
import glob

from argparse import ArgumentParser
from collections import Counter, defaultdict

import pybedtools

from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists, splitext_plus, tmpfile, safe_makedir
from bcbio.install import _get_data_dir

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import gzip

from collections import Counter
import pandas as pd
# from ggplot import *

import vcf


def _sv_dist(fn_in):
    count = Counter()
    vcf_reader = vcf.Reader(open(fn_in, 'r'))
    samples = vcf_reader.samples
    remove_sr = remove_common = 0
    for record in vcf_reader:
        if record.genotype(samples[0])['GT'] == '0/0':
            continue
        if record.genotype(samples[1])['GT'] == '0/0':
            if record.genotype(samples[0])['SR'] > 5:
                # print record.genotype(samples[0])['GT']
                # print record.genotype(samples[1])['GT']
                # print record.ALT
                count[(record.CHROM, record.var_subtype)] += 1
            else:
                remove_sr += 1
        else:
            remove_common += 1
    print "Removed common %s, removed SR %s" % (remove_sr, remove_common)
    return count

def _parse_count(count, sample):
    tab = []
    row = 0
    for k in count:
        row +=1
        tab.append([sample, k[0], k[1], count[k]])
    return  tab

def simple_report(data, args):
    summary = []
    for sample in data:
        print sample[0]['files']['lumpy-pair.vcf']
        dt = _sv_dist(sample[0]['files']['lumpy-pair.vcf'])
        dt = pd.DataFrame(_parse_count(dt, sample[0]['name']))
        # print dt.columns
        dt.columns = ['sample', 'chrom', 'sv', 'counts']
        # p = ggplot(aes(x="chrom", y="counts", fill="sv"), data=dt) + geom_bar(position = 'dodge', stat = 'identity')
        out_file = op.join(args.out, sample[0]['name'] + "_lumpy.tsv")
        dt.to_csv(out_file, sep='\t', index=False)
        summary.append(dt)
        #ggsave(p, sample[0]['name'] + "_lumpy.png")

    out_file = op.join(args.out, "lumpy.tsv")
    pd.concat(summary).to_csv(out_file, sep='\t', index=False)

def _get_samples(out_dir):
    data = defaultdict(dict)
    for fn in glob.glob(op.join(out_dir,'*/*')):
        if fn.endswith('tbi') or fn.endswith('bai'):
            continue
        rel_path = fn.split(out_dir)[1][1:]
        sample = rel_path.split("/")[0]
        fn_type, ext = splitext_plus(rel_path.split("/")[1].replace(sample + "-", ""))
        if sample.find("Tumor") > -1:
            data[sample][fn_type + ext.replace(".gz", "")] = fn
    return data

def _prepare_samples(args):
    """
    create dict for each sample having all information
    """
    # if args.galaxy:
    #    system_config = args.galaxy
    # else:
    #    system_config = os.path.join(_get_data_dir(), "galaxy", "bcbio_system.yaml")
    # config = yaml.load(open(system_config))
    # config['algorithm'] = {}
    data = []
    raw = _get_samples(args.files[0])
    for sample in raw:
        dt = {}
        dt['name'] = sample
        dt['files'] = raw[sample]
        # dt['config'] = config
        data.append([dt])
    return data


if __name__ == "__main__":
    parser = ArgumentParser(description="clean SV VCF files.")
    parser.add_argument("--region", help="bed file with regions.")
    parser.add_argument("files", nargs="*", help="Bam files.")
    parser.add_argument("--run", required=1, help="Clean VCF file.", choices=['report'])
    args = parser.parse_args()

    safe_makedir(args.out)
    if args.run == 'report':
        data = _prepare_samples(args)
        simple_report(data, args)
