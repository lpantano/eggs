"""
Run multiple tools to measure accuracy of isomiRs genome based annotation
"""

from argparse import ArgumentParser
import os
import contextlib

from os.path import exists as is_there
from os.path import abspath as full
from bcbio.provenance import do


@contextlib.contextmanager
def ch_directory(dir):
    cur_dir = os.getcwd()
    os.chdir(dir)
    yield dir
    os.chdir(cur_dir)


def _annotate(input, mirbase):
    output = "mirbase.bed"
    cmd = ("bedtools intersect -bed -wo -s -f 0.80 -abam"
           " {input} -b {mirbase} >| {output}")
    if not is_there(output):
        do.run(cmd.format(**locals()), "")
    return output


def _star(input, index, mirbase):
    with ch_directory("star"):
        output = "star_map.bam"
        cmd = ("STAR --genomeDir {index} --readFilesIN {input}"
               " --outFilterMultimapNmax 50"
               " --outSAMattributes NH HI NM"
               " --alignIntronMax 1")
        cmd_bam = "samtools view -Sbh Aligned.out.sam >| {output}"

        if not is_there("Aligned.out.sam"):
            do.run(cmd.format(**locals()), "")
        if not is_there(output):
            do.run(cmd_bam.format(**locals()), "")

        mirbase_output = _annotate(output, mirbase)
    return mirbase_output


def _bowtie2(input, index, mirbase):
    with ch_directory("bowtie2"):
        output = "bowtie2_map.bam"
        cmd = ("bowtie2 -k 50 -L 18 -x {index}"
               " -U {input}"
               " >| hits.sam")
        cmd_bam = "samtools view -Sbh hits.sam >| {output}"

        if not is_there("hits.sam"):
            do.run(cmd.format(**locals()), "")
        if not is_there(output):
            do.run(cmd_bam.format(**locals()), "")

        mirbase_output = _annotate(output, mirbase)
    return mirbase_output


def _hisat(input, index, mirbase):
    with ch_directory("hisat"):
        output = "hisat_map.bam"
        cmd = ("hisat -k 50 -L 18 -x {index}"
               " -U {input}"
               " >| hits.sam")
        cmd_bam = "samtools view -Sbh hits.sam >| {output}"

        if not is_there("hits.sam"):
            do.run(cmd.format(**locals()), "")
        if not is_there(output):
            do.run(cmd_bam.format(**locals()), "")

        mirbase_output = _annotate(output, mirbase)
    return mirbase_output


if __name__ == "__main__":
    parser = ArgumentParser(description="Run different tools in simulated fasta file")
    parser.add_argument("--fasta", required=True, help="short reads")
    parser.add_argument("--mirbase", required=True, help="bed file with mirbase annotation")
    parser.add_argument("--star", help="star index")
    parser.add_argument("--bowtie2", help="bowtie2 index")
    parser.add_argument("--hisat", help="hisat index")
    args = parser.parse_args()

    if args.star:
        print "doing STAR"
        _star(full(args.fasta), full(args.star), full(args.mirbase))
    if args.bowtie2:
        print "doing bowtie2"
        _bowtie2(full(args.fasta), full(args.bowtie2), full(args.mirbase))
    if args.hisat:
        print "doing hisat"
        _hisat(full(args.fasta), full(args.hisat), full(args.mirbase))
