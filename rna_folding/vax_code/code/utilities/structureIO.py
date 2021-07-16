"""
Helper functions for processing IO from various structure prediction tools.
"""


import numpy as np
import re


def get_bpp_matrix_from_linearpartition_output(file_name, sequence):
    """Return a BPP (base pairing probabilities) matrix from LinearPartition output.

    See: https://github.com/LinearFold/LinearPartition

    :param file_name: LinearPartition output file containing base pairing probabilities
    :param sequence: RNA sequence that was used as input fo LinearPartition run
    :return: numpy array of shape sequence length * sequence length containing base pairing probabilities
    """

    bpp_matrix = np.zeros([len(sequence), len(sequence)])

    with open(file_name,'r') as f:
        for line in f.readlines():
            if (not re.match("^\s*$", line)) and (not line.startswith(">")):
                ind_1, ind_2, prob = line.strip().split(' ')
                ind_1 = int(ind_1)-1
                ind_2 = int(ind_2)-1
                prob = float(prob)
                bpp_matrix[ind_1, ind_2] = prob
                bpp_matrix[ind_2, ind_1] = prob

    return bpp_matrix
