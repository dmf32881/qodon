"""
Objects for representing RNAs and calculating sequence and structure-based metrics.
"""
from typing import List, Any

from scipy.stats import gmean
import numpy as np

from Bio.Seq import Seq
from dnachisel import DnaOptimizationProblem, AvoidHairpins

from codeRNA.codon_data import cai_weights
from codeRNA.codon_data import codon_tables
from codeRNA.codon_data import hive_counts


class RNA(Seq):
    """A class for representing RNA sequence and (optionally) structure.

    Provides methods for calculating various structure and sequence related metrics (such as get_uridine_content,
    get_cai, etc.).

    Inherits from Biopython.Seq and therefore also provides a number of string like methods (such as count,
    find, split and strip) and some biological methods, such as complement, reverse_complement, back_transcribe,
    and translate.
    """

    def __init__(self, sequence, mfe=None, bpp_matrix=None):
        # replace T's to U's to make into RNA sequence
        sequence = sequence.upper().replace('T', 'U')
        super().__init__(sequence)
        self.sequence = sequence
        self.mfe = None
        self.bpp_matrix = None
        self.bp_count = len(self.sequence)
        # check if length is divisible by 3 (required e.g. for CAI calculations)
        self.is_coding = self.bp_count % 3 == 0
        self.set_mfe(mfe)
        self.protein_sequence = None
        if self.is_coding:
            self.protein_sequence = self.translate()
        self.set_bpp_matrix(bpp_matrix)

    def translate(self, table="Standard", stop_symbol="*", to_stop=False):
        """Return a protein sequence translated from the RNA sequence.

        :param stop_symbol: Single character string, what to use for terminators (typically “*” or "s")
        :param to_stop: Boolean, defaults to False meaning do a full translation continuing on past any stop codons
            (translated as the specified stop_symbol). If True, translation is terminated at the first in frame stop
            codon (and the stop_symbol is not appended to the returned protein sequence).
        :return: string, protein sequence
        """
        return str(Seq.translate(self.sequence, table="Standard", stop_symbol=stop_symbol, to_stop=to_stop))

    def get_uridine_content(self):
        """
        Calculate uridine content for the RNA sequence
        :return: uridine content (fraction of uridines in the sequence)
        """
        u_content = self.count("U") / self.bp_count
        return u_content

    def get_gc_content(self):
        gc_content = (self.count("G") + self.count("C")) / self.bp_count
        return gc_content

    def get_cai(self, organism="human"):
        """
        Return CAI (codon adaptation index) value based on codon usage data for a given organism.
        Raises an if the sequence is not divisible by 3 or if the selected organism is not supported.

        :param organism: target organism, defaults to "human"
        :return: codon adaptation index
        """
        if not self.is_coding:
            raise ValueError("Cannot calculate CAI for a non-coding sequence. Length must be a multiple of 3.")

        codon_weights = cai_weights.get_codon_weights(organism=organism)

        w_list = []
        for i in range(0, self.bp_count, 3):
            codon = self.sequence[i:i + 3]
            # Do not count W or M codon since there is only one that encodes them
            if codon not in ['UGG', 'AUG']:
                w_list.append(codon_weights[codon])

        return gmean(w_list)

    def get_gc3_content(self):
        """Get GC3 content metric (average percent GC content in the third (wobble) codon position).

        :return: float GC3 content
        """
        if not self.is_coding:
            raise ValueError("Cannot calculate GC3 for a non-coding sequence. Length must be a multiple of 3.")
        wobble_nucleotides = [self.sequence[i] for i in range(2, self.bp_count, 3)]
        return (wobble_nucleotides.count("G") + wobble_nucleotides.count("C")) / len(wobble_nucleotides)

    def get_median_codon_frequency(self, organism="human"):
        """
        Returns the median codon frequency for the RNA and sepecified organism.
        Raises a ValueError if RNA is non-coding (read: not divisible by 3).
        Metric originally defined in: https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-00968-8

        :param organism: Organism defining codon frequency table.
        :return: Median codon frequency.
        """
        if not self.is_coding:
            raise ValueError(
                "Cannot calculate median codon frequency for a non-coding sequence. Length must be a multiple of 3.")

        codon_frequency_table = hive_counts.get_codon_frequencies(organism)
        codon_frequencies = [codon_frequency_table[self.sequence[i:i + 3]] for i in range(0, self.bp_count, 3)]

        median_codon_frequency = np.median(codon_frequencies)

        return median_codon_frequency

    def get_mfe(self):
        return self.mfe

    def get_min_max(self, window=18, organism="human"):
        """
        Returns the %MinMax for a sliding window along the RNA and a specified organism.
        Raises a ValueError if RNA is non-coding (not divisible by 3).
        Metric originally defined in: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003412

        :param window: Range over which to take averages for min/max.
        :param organism: Organism defining codon frequency table.
        :return: List of %MinMax
        """
        # maps amino acid to list of all codons
        codon_table = codon_tables.get_codon_table()
        # frequency of each codon vs all others
        codon_frequency_table = hive_counts.get_overall_codon_frequencies(organism)
        # codon_frequency_table = hive_counts.get_codon_frequencies(organism)

        # most/least common and avg. codon frequencies coding for the same AA sequence
        max_frequencies = []
        min_frequencies = []
        avg_frequencies = []

        # self.translate() overwrite self.sequence in python3.6, so avoid it for now
        for aa in Seq.translate(self.sequence):
            # could store as a table just look up after populating once if this is slow
            frequencies_for_aa = [codon_frequency_table[x] for x in codon_table[aa]]
            max_frequencies += [max(frequencies_for_aa)]
            min_frequencies += [min(frequencies_for_aa)]
            avg_frequencies += [np.mean(frequencies_for_aa)]

        codon_frequencies = [codon_frequency_table[self.sequence[i:i + 3]] for i in range(0, self.bp_count, 3)]

        minmax = []

        for i in range(len(codon_frequencies) - window + 1):
            if sum(codon_frequencies[i:i + window]) > sum(avg_frequencies[i:i + window]):
                minmax_i = sum([codon_frequencies[i + j] - avg_frequencies[i + j] for j in range(window)])
                minmax_i /= sum([max_frequencies[i + j] - avg_frequencies[i + j] for j in range(window)])
                minmax.append(100 * minmax_i)
            else:
                minmax_i = sum([avg_frequencies[i + j] - codon_frequencies[i + j] for j in range(window)])
                minmax_i /= sum([avg_frequencies[i + j] - min_frequencies[i + j] for j in range(window)])
                minmax.append(-100 * minmax_i)  # convention from paper

        return minmax

    def get_min_window(self, window=18, organism="human"):
        """Return %Min window, percent of the nucleotide sequence that exists with a minimum below zero in the MinMax
        profile (clusters of rare codons).
        Metric originally defined in: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003412

        :param window: Range over which to take averages for min/max.
        :param organism: Organism defining codon frequency table.
        :return: float, %Min window
        """
        minmax = self.get_min_max(window=window, organism=organism)
        below_0 = [el for el in minmax if el < 0]
        return len(below_0) / len(minmax)

    def get_aup(self):
        """Return the average unpaired probability (AUP) metric.

        Calculates AUP as introduced by Wayment-Steele et al.
        (https://www.biorxiv.org/content/10.1101/2020.08.22.262931v2)

        :return: float aup value
        """
        if self.bpp_matrix is None:
            raise ValueError("Cannot calculate AUP for an RNA molecule with unknown structure. Run set_bpp_matrix().")
        prob_unp_vector = 1 - np.sum(self.bpp_matrix, axis=0)
        aup = np.mean(prob_unp_vector)
        return aup

    def get_hairpin_count(self, stem_size=20, hairpin_window=200):
        """Return number of hairpin patterns in the sequence.
        A hairpin is defined by a sequence segment which has a reverse complement “nearby” in a given window.
        See: https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html and
        github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/master/dnachisel/builtin_specifications/AvoidHairpins.py

        :param stem_size: Size of the stem of a hairpin, i.e. the length of the sequence which should have a reverse
        complement nearby to be considered a hairpin.
        :param hairpin_window: The window in which the stem’s reverse complement should be searched for.
        :return: int number of hairpin patterns in the sequence
        """
        evaluation = self._get_hairpin_info(stem_size=stem_size, hairpin_window=hairpin_window)
        return -evaluation.score  # DNAChisel assigns hairpin score as -number_of_hairpins

    def get_hairpin_locations(self, stem_size=20, hairpin_window=200):
        """Return the locations of hairpin patterns in the sequence. Returns a list of tuples representing a segment
        of the sequence (start, end).

        Warning: function is based on dnachisel, which uses Python's splicing notation, so location (5, 10) represents
        sequence[5, 6, 7, 8, 9] which corresponds to nucleotides number 6, 7, 8, 9, 10.

        A hairpin is defined by a sequence segment which has a reverse complement “nearby” in a given window.
        See: https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html and
        github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/master/dnachisel/builtin_specifications/AvoidHairpins.py

        :param stem_size: Size of the stem of a hairpin, i.e. the length of the sequence which should have a reverse
        complement nearby to be considered a hairpin.
        :param hairpin_window: The window in which the stem’s reverse complement should be searched for.
        :return: list of (int:start, int:end) tuples representing  segments of the sequence
        """
        evaluation = self._get_hairpin_info(stem_size=stem_size, hairpin_window=hairpin_window)
        locations = [(loc.start, loc.end) for loc in evaluation.locations]
        return locations

    def _get_hairpin_info(self, stem_size=20, hairpin_window=200):
        problem = DnaOptimizationProblem(
            sequence=self.sequence.replace("U", "T"),  # works on DNA sequences
            constraints=[AvoidHairpins(stem_size=stem_size, hairpin_window=hairpin_window)]
        )
        evaluation = problem.constraints_evaluations().evaluations[0]
        return evaluation

    def set_bpp_matrix(self, bpp_matrix):
        """Set the bpp_matrix (base pairing probabilities).
        Raises a Value error if the matrix is not square of if the dimensions do not match the length of the sequence

        :param bpp_matrix: numpy array with base pair probabilities
        """
        if bpp_matrix is not None:
            if bpp_matrix.shape[0] != bpp_matrix.shape[1]:
                raise ValueError("Incorrect dimensions of the input bpp_matrix - must be a square matrix")
            if bpp_matrix.shape[0] != self.bp_count:
                raise ValueError("Dimensions of bpp_matrix don't match the sequence length - must be {0}x{0}".format(
                    self.bp_count))
        self.bpp_matrix = bpp_matrix

    def set_mfe(self, mfe):
        """Set mfe value (Minimum Free Energy).
        Raises a Value error if mfe is larger than 0

        :param mfe: float Minimum Free Energy value
        """
        if (mfe is not None) and (mfe > 0):
            raise ValueError("Incorrect value of MFE: must be < 0")
        self.mfe = mfe


if __name__ == "__main__":
    seq = RNA("UTCAGCAUCACG")
    print("U content: {}".format(seq.get_uridine_content()))
    print("GC content: {}".format(seq.get_gc_content()))
    print("GC3 content: {}".format(seq.get_gc3_content()))
    print("CAI: {}".format(seq.get_cai()))
    print(f"Median codon frequency: {seq.get_median_codon_frequency()}")
    print(f"Minmax: {seq.get_min_max(3)}")
    print(f"Hairpin count: {seq.get_hairpin_count(stem_size=3, hairpin_window=8)}")
    print(f"Hairpin locations: {seq.get_hairpin_locations(stem_size=3, hairpin_window=8)}")
