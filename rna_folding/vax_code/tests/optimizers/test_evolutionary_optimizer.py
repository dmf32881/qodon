"""
Unit tests for codeRNA/optimizers/evolutionary_optimizer.py
"""
import pytest

from codeRNA.optimizers.evolutionary_optimizer import EvolutionaryRNAOptimizer
from codeRNA.rna.rna import RNA


@pytest.fixture
def rna_evo():
    pytest.aa_seq = "MEDAKNIKKGPAPFYPLEDGTAGEQLHKAMKRYALVPGTIAFTDAHIEVNITYAEYFEMSVRLAEAMKRYGLNTNHRIVVCSENSLQFFMPVLGALFIGVAVAPA"
    pytest.get_u_content = lambda seq: seq.count("U") / len(seq)
    pytest.population_size = 50
    pytest.rna_evo = EvolutionaryRNAOptimizer(protein_sequence=pytest.aa_seq, scoring_function=pytest.get_u_content,
                                              population_size=pytest.population_size)


def test_get_initial_sequences(rna_evo):
    # test that the initial population is of right size
    initial_population = pytest.rna_evo.get_initial_sequences()
    assert len(initial_population) == pytest.population_size

    # test that all sequences in initial population code for the same protein
    initial_population_protein_seqs = [RNA(seq).protein_sequence for seq, score in initial_population]
    assert all([seq == pytest.aa_seq for seq in initial_population_protein_seqs])

    # test that scores are correctly sorted
    scores = [score for seq, score in initial_population]
    assert all(scores[i] <= scores[i + 1] for i in range(len(scores) - 1))


def test_get_codon_map(rna_evo):
    # test that codon map has the right number of amino acids
    codon_map = pytest.rna_evo._get_codon_map()
    assert len(codon_map) == 21  # Including stop

    # test that codon map has the right number of codons
    codons = [x for el in codon_map.values() for x in el["codons"]]
    assert len(codons) == 64
    assert len(set(codons)) == 64


def test_get_optimized_sequences(rna_evo):
    # test that function raises ValueError if called too early
    with pytest.raises(ValueError):
        pytest.rna_evo.get_optimized_sequences()

    # test that number of optimized sequences is equal to the population size
    pytest.rna_evo.optimize(max_iterations=10)
    optimized_seqs = pytest.rna_evo.get_optimized_sequences()
    assert len(optimized_seqs) == pytest.population_size

    # test that all sequences in initial population code for the same protein
    optimized_protein_seqs = [RNA(seq).protein_sequence for seq, score in optimized_seqs]
    assert all([seq == pytest.aa_seq for seq in optimized_protein_seqs])

    # test that scores are correctly sorted
    scores = [score for seq, score in optimized_seqs]
    assert all(scores[i] <= scores[i + 1] for i in range(len(scores) - 1))


def test_optimize(rna_evo):
    # test that scores are improved after optimization run
    initial_scores = [score for seq, score in pytest.rna_evo.get_initial_sequences()]
    initial_score_min = min(initial_scores)

    pytest.rna_evo.optimize(max_iterations=20)
    optimized_scores = [score for seq, score in pytest.rna_evo.get_optimized_sequences()]
    optimized_score_min = min(optimized_scores)

    assert initial_score_min > optimized_score_min
