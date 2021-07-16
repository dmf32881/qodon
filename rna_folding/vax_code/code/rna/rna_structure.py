"""
Functions and wrappers for predicting RNA secondary structure.
"""
import os
import subprocess
import re

from codeRNA.utilities.structureIO import get_bpp_matrix_from_linearpartition_output


def run_linear_partition_get_bpp(sequence, linear_partition_home=os.getenv("LINEAR_PARTITION_HOME"), beam_size=100):
    """Run LinearPartition and return a BPP (base pairing probabilities) matrix.
    See: https://github.com/LinearFold/LinearPartition

    :param sequence: nucleotide sequence
    :param beam_size: The beam size (default 100). Use 0 for infinite beam.
    :param linear_partition_home: path to linear partition installation, if not provided, the function will look for
        $LINEAR_PARTITION_HOME
    :return: numpy array of shape sequence length * sequence length containing base pairing probabilities
    """
    tmp_file = os.path.join(os.getcwd(), "tmp_file.out")
    p1 = subprocess.Popen(["echo", sequence], stdout=subprocess.PIPE)
    p2 = subprocess.Popen([os.path.join(linear_partition_home, 'linearpartition'), "-V", "--verbose", "-o",
                           tmp_file, "-b", str(beam_size)], stdin=p1.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    p2.communicate()
    bpp_matrix = get_bpp_matrix_from_linearpartition_output(tmp_file, sequence)
    if os.path.exists(tmp_file):
        os.remove(tmp_file)
    return bpp_matrix


def run_linear_partition_get_mfe(sequence, linear_partition_home=os.getenv("LINEAR_PARTITION_HOME"), beam_size=100):
    """Run LinearPartition and return MFE (Minimum Free Energy) value.
    See: https://github.com/LinearFold/LinearPartition

    :param sequence: nucleotide sequence
    :param linear_partition_home: path to linear partition installation, if not provided, the function will look for
        $LINEAR_PARTITION_HOME
    :param beam_size: The beam size (default 100). Use 0 for infinite beam.
    :return: float MFE
    """
    p1 = subprocess.Popen(["echo", sequence], stdout=subprocess.PIPE)
    p2 = subprocess.Popen([os.path.join(linear_partition_home, 'linearpartition'), "-V", "-p", "-b", str(beam_size)],
                          stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, error = p2.communicate()
    mfe = re.search(rb'Free Energy of Ensemble: (.*) kcal/mol', error).group(1)
    return float(mfe)


if __name__ == "__main__":
    seq = "GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACC"
    bpp_matrix = run_linear_partition_get_bpp(sequence=seq, beam_size=100)
    print(bpp_matrix)
    mfe = run_linear_partition_get_mfe(sequence=seq, beam_size=100)
    print(mfe)