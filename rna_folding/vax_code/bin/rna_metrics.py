"""
Simple command-line utility for calculating basic RNA metrics for a single sequence.

usage: rna_metrics.py [-h] [--linear_partition_home LINEAR_PARTITION_HOME]
                      [--beam_size BEAM_SIZE]
                      fasta_file

Calculate some RNA metrics.

positional arguments:
  fasta_file            Path to a fasta file

optional arguments:
  -h, --help            show this help message and exit
  --linear_partition_home LINEAR_PARTITION_HOME
                        Path to LinearPartition. If not provided,
                        $LINEAR_PARTITION_HOME will be used.If this is also
                        missing, AUP and MFE will not be calculated.
  --beam_size BEAM_SIZE
                        The beam size to be used by LinearPartition (default
                        100). Use 0 for infinite beam.
"""

import argparse
import os

from Bio import SeqIO

from codeRNA.rna import rna
from codeRNA.rna.rna_structure import run_linear_partition_get_bpp, run_linear_partition_get_mfe

parser = argparse.ArgumentParser(description='Calculate some RNA metrics.')
parser.add_argument('fasta_file', type=str, help='Path to a fasta file')
parser.add_argument('--linear_partition_home', type=str, default=None,
                    help='Path to LinearPartition. If not provided, $LINEAR_PARTITION_HOME will be used.'
                         'If this is also missing, AUP and MFE will not be calculated.')
parser.add_argument('--beam_size', type=int, default=100,
                    help='The beam size to be used by LinearPartition (default 100). Use 0 for infinite beam.')
args = parser.parse_args()

if __name__ == "__main__":
    seq = SeqIO.read(args.fasta_file, "fasta").seq

    linear_partition_home = args.linear_partition_home
    if linear_partition_home is None:
        linear_partition_home = os.getenv("LINEAR_PARTITION_HOME")
    if linear_partition_home is None:
        print("Path to linear partition is missing. Use --linear_partition_home_argument.")
        bpp_mat = None
    else:
        print("Using the following installation of Linear Partition to calculate secondary structue: {}".format(
            linear_partition_home
        ))
        bpp_mat = run_linear_partition_get_bpp(str(seq), linear_partition_home=linear_partition_home,
                                               beam_size=args.beam_size)
        mfe = run_linear_partition_get_mfe(str(seq), linear_partition_home=linear_partition_home,
                                               beam_size=args.beam_size)
    r = rna.RNA(str(seq), bpp_matrix=bpp_mat, mfe=mfe)
    print("U content: {}".format(r.get_uridine_content()))
    print("GC content: {}".format(r.get_gc_content()))
    print("GC3 content: {}".format(r.get_gc3_content()))
    print("CAI: {}".format(r.get_cai()))
    print("%Min window: {}".format(r.get_min_window()))
    print("Hairpin count: {}".format(r.get_hairpin_count()))
    if r.bpp_matrix is not None:
        print("AUP: {}".format(r.get_aup()))
    if r.mfe is not None:
        print("MFE: {}".format(r.get_mfe()))
