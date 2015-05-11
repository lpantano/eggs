import os.path as op

from ichwrapper.log import run
from bcbio.distributed.transaction import file_transaction, tx_tmpdir

def trymemory(data, args):
    bed_file = data['bed']
    mem = data['config']['resources']['mem']
    out_file = op.join("mem%s" % mem, op.basename(bed_file))
    fails = []
    cmd = "export VAR_DICT_OPTS='-Xms750m -Xmx{mem}000m -XX:+UseSerialGC -Djava.io.tmpdir={temp_dir}' && vardict-java -G /cm/shared/apps/bcbio/20150103-devel/data/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -f 0.1 -N tumor-70_0001 -b \"/scratch/lpantano/montreal/analysis/work/align/tumor-70_0001/7_2015-05-08_analysis-sort.bam|/scratch/lpantano/montreal/analysis/work/align/norm-70_0001/3_2015-05-08_analysis-sort.bam\" -c 1 -S 2 -E 3 -g 4 {bed_file}"
    try:
        with file_transaction(out_file) as tx_out:
            with tx_tmpdir() as temp_dir:
                run(cmd.format(**locals()))

    except:
        fails = data['bed']
    return fails
