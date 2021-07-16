"""
Unit tests for codeRNA/rna/rna.py
The general rule of thumb is to test individual function.
Try to test for a typical use case and an edge case. 
As lots of data/code comes from publications, maybe try to replicate a published result.

Integration testing to be determined.
"""
from codeRNA.rna.rna import RNA
import numpy as np
import pytest

# sequences from https://github.com/JeffreyJLi/codon_optimization_analysis
# other tables/data from https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-00968-8
@pytest.fixture
def idt_seq_from_paper_as_rna():
    idt = RNA(
        "ATGACTGAGTATAAACTTGTTGTCGTCGGTGCGGGTGGCGTCGGCAAGTCCGCACTTACAATCCAATTGA"
        "TCCAAAACCACTTTGTCGACGAGTACGATCCAACTATTGAGGATAGTTACCGCAAGCAGGTTGTTATCGA"
        "TGGAGAAACGTGTTTGTTGGACATTTTGGATACAGCAGGCCAGGAGGAGTACAGTGCCATGCGCGATCAA"
        "TATATGCGCACCGGGGAAGGTTTTTTGTGTGTGTTTGCCATTAACAACACGAAGAGTTTTGAAGACATCC"
        "ACCACTATCGTGAACAAATCAAGCGTGTCAAAGATAGCGAAGATGTCCCCATGGTATTAGTAGGGAATAA"
        "ATGTGACTTACCGAGCCGTACCGTTGACACAAAGCAGGCACAAGACCTGGCACGCTCATACGGTATTCCC"
        "TTCATCGAAACTTCCGCCAAAACACGTCAAGGAGTAGACGACGCTTTCTACACGCTTGTTCGCGAGATTC"
        "GTAAACATAAGGAGAAGATGAGCAAAGATGGGAAGAAAAAGAAAAAAAAGTCAAAAACAAAGTGCGTCAT"
        "TATGTAA")
    return idt

@pytest.fixture
def twist_seq_from_paper_as_rna():
    twist = RNA(
        "ATGACAGAGTACAAGCTCGTCGTTGTGGGTGCAGGCGGTGTGGGGAAATCTGCGCTCACCATTCAACTTA"
        "TCCAAAACCACTTCGTAGATGAGTACGACCCCACCATTGAAGACAGTTATCGCAAACAGGTGGTGATCGA"
        "CGGCGAGACGTGCTTGCTCGACATCCTGGATACGGCGGGACAGGAAGAATATAGCGCGATGCGTGATCAA"
        "TATATGCGTACCGGCGAAGGATTCTTGTGCGTCTTCGCAATTAACAACACCAAGAGCTTCGAGGACATCC"
        "ATCACTACCGCGAGCAGATCAAGCGCGTGAAAGATAGCGAGGACGTGCCGATGGTATTGGTCGGCAACAA"
        "GTGCGACCTGCCATCACGCACCGTTGATACGAAGCAAGCCCAAGATCTTGCCCGCAGCTACGGTATCCCA"
        "TTCATCGAGACCTCTGCCAAAACCCGTCAAGGCGTGGACGACGCATTTTACACCCTGGTGCGCGAGATCC"
        "GTAAGCACAAGGAGAAAATGTCGAAGGATGGAAAGAAGAAGAAGAAGAAATCCAAAACTAAATGCGTCAT"
        "CATGTGA")
    return twist

@pytest.fixture
def native_seq_from_paper_as_rna():
    native = RNA(
        "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTC"
        "AGAATCATTTTGTGGACGAATATGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGA"
        "AACCTGTCTCTTGGATATTCTCGACACAGCAGGTCAAGAGGAGTACAGTGCAATGAGGGACCAGTACATGAGG"
        "ACTGGGGAGGGCTTTCTTTGTGTATTTGCCATAAATAATACTAAATCATTTGAAGATATTCACCATTATAGAG"
        "AACAAATTAAAAGAGTTAAGGACTCTGAAGATGTACCTATGGTCCTAGTAGGAAATAAATGTGATTTGCCTTC"
        "TAGAACAGTAGACACAAAACAGGCTCAGGACTTAGCAAGAAGTTATGGAATTCCTTTTATTGAAACATCAGCA"
        "AAGACAAGACAGGGTGTTGATGATGCCTTCTATACATTAGTTCGAGAAATTCGAAAACATAAAGAAAAGATGA"
        "GCAAAGACGGTAAAAAGAAGAAAAAGAAGTCAAAGACAAAGTGTGTAATTATGTAA")
    return native


def test_RNA_class():
    # init
    rna = RNA("ATGCACTGTTGA")
    # rna class should replace Ts with Us
    assert rna.sequence == "AUGCACUGUUGA"
    # rna class should translate
    assert rna.translate() == "MHC*"
    # add other tests as necessary


def test_get_uridine_content():
    # init
    rna = RNA("ATGCACTGTTGA")
    assert rna.get_uridine_content() == 4 / 12


def test_get_gc_content(idt_seq_from_paper_as_rna, twist_seq_from_paper_as_rna):
    # init
    rna = RNA("ATGCACTGTTGA")
    assert rna.get_gc_content() == 5 / 12
    # paper tests
    idt = idt_seq_from_paper_as_rna
    twist = twist_seq_from_paper_as_rna
    assert np.isclose(idt.get_gc_content(), 0.451, atol=0.001)
    assert np.isclose(twist.get_gc_content(), 0.515, atol=0.001)

def test_get_gc3_content(idt_seq_from_paper_as_rna, twist_seq_from_paper_as_rna):
    # init
    rna = RNA("ATGCACTGTTGA")
    assert rna.get_gc3_content() == 2 / 4
    # paper tests
    idt = idt_seq_from_paper_as_rna
    twist = twist_seq_from_paper_as_rna
    assert np.isclose(idt.get_gc3_content(), 0.50, atol=0.02)
    assert np.isclose(twist.get_gc3_content(), 0.67, atol=0.02)


def test_get_cai():
    return


def test_get_median_codon_frequency(idt_seq_from_paper_as_rna, twist_seq_from_paper_as_rna):
    # only paper tests
    idt = idt_seq_from_paper_as_rna
    twist = twist_seq_from_paper_as_rna
    assert np.isclose(idt.get_median_codon_frequency("ecoli_pub"), 0.373, atol=0.001)
    assert np.isclose(twist.get_median_codon_frequency("ecoli_pub"), 0.373, atol=0.001)

def test_get_mfe():
    return


def test_get_min_max():
    """
    Cannot test against the paper as the paper's codon usage table does not seem to exist anymore.
    Test are against manually calculated values.
    """
    minirna = RNA("AUGAAGUCGAGGACC")
    # values with window of 5 averages
    # min : [27.44, 10.84, 7.82, 1.70, 7.86] = 11.13
    # avg : [27.44, 22.31, 9.91, 9.34, 13.57] = 16.51
    # max : [27.44, 33.78, 15.95, 21.39, 22.75] = 22.26
    # act : [27.44, 10.84, 8.87, 1.70, 22.75] = 14.32
    # minmax% = -100 * (16.51 - 14.32) / (16.51 - 11.13) = -40.7
    results = minirna.get_min_max(window=5, organism="ecoli")
    assert np.isclose(results, -40.7, atol=0.1)


def test_get_aup():
    return


def test_translate():
    rna = RNA("AUGAAGUCGAGGACC")
    assert rna.translate() == "MKSRT"
    rna = RNA("ACAGAAAGCACAACAAGAGCAAACAGCCUAGCAACAGAA")
    assert rna.translate() == "TESTTRANSLATE"
    # translate until stop
    rna = RNA("ACAGAAAGCACAUAAACAAGAGCAAACAGCCUAGCAACAGAA")
    assert rna.translate(to_stop=True) == "TEST"


def test_get_min_window(idt_seq_from_paper_as_rna, twist_seq_from_paper_as_rna):
    # only paper tests
    idt = idt_seq_from_paper_as_rna
    twist = twist_seq_from_paper_as_rna
    assert np.isclose(idt.get_min_window(window=18, organism="ecoli"), 0.27, atol=0.05)
    assert np.isclose(twist.get_min_window(window=18, organism="ecoli"), 0.17, atol=0.05)


def test_get_hairpin_count():
    # test on a sequence with hairpins on extremities
    # see https://github.com/Edinburgh-Genome-Foundry/DnaChisel/issues/37
    minirna = RNA("AUUCAAUGGGGGGGGGGGGGGGGGGGGGGGGGUAGCCUA")
    assert minirna.get_hairpin_count(stem_size=3, hairpin_window=8) == 2


def test_get_hairpin_locations():
    # test on a sequence with hairpins on extremities
    # see https://github.com/Edinburgh-Genome-Foundry/DnaChisel/issues/37
    minirna = RNA("AUUCAAUGGGGGGGGGGGGGGGGGGGGGGGGGUAGCCUA")
    locations = minirna.get_hairpin_locations(stem_size=3, hairpin_window=8)
    assert len(locations) == 2
    assert locations[0] == (0, 6)
    assert locations[1] == (32, 39)