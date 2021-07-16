"""
Evolutionary algorithm for optimizing RNA sequences.
The algorithm iteratively improves the initial population of candidate sequences by applying genetic operators of
mutation and recombination.

See: https://www.tensorflow.org/probability/api_docs/python/tfp/optimizer/differential_evolution_minimize

"""

import random
import timeit
from datetime import timedelta
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import tensorflow_probability as tfp
import tensorflow as tf
import os, re, tempfile, subprocess

from codeRNA.codon_data.hive_counts import get_codon_frequencies
from codeRNA.codon_data.codon_tables import get_codon_table
from codeRNA.rna.rna import RNA


class EvolutionaryRNAOptimizer(object):
    def __init__(self, protein_sequence, scoring_function, population_size=100,
                 organism="human", codon_table_name="standard"):
        """Create an EvolutionaryRNAOptimizer object.

        Example scoring function (minimizing uridine content):
        lambda seq: RNA(seq).get_uridine_content()

        :param protein_sequence: str, protein sequence
        :param scoring_function: function used to calculate scores for individual sequences (to be minimized by the
            evolutionary algorithm)
        :param population_size: int, the size of the population to evolve, must be larger than 4
        :param organism: str, target organism, defaults to "human"
        :param codon_table_name: str, name of the codon table, defaults to "standard"
        """
        self.protein_sequence = protein_sequence
        self.scoring_function = scoring_function
        self.organism = organism
        self.codon_table_name = codon_table_name
        self.par = True
        self.population_size = population_size
        self.codon_map = self._get_codon_map()
        self.initial_population_tensor = self._generate_initial_population()
        self.optimization_results = None

    def optimize(self, max_iterations=500, differential_weight=0.01, crossover_prob=0.025, seed=None):
        """Run evolutionary algorithm.
        The results can be checked by running get_ranked_optimized_sequences().

        :param max_iterations: int, maximum  number of iterations to evolve the population for
        :param differential_weight: float, parameter controlling strength of mutation. Must be positive and less than 2
        :param crossover_prob: float, the probability of recombination per site; must be between 0 and 1
        :param seed: int, random seed. If None, no seed is applied
        """

        start = timeit.default_timer()
        if self.scoring_function != 'linearpartition':
            optimization_results = tfp.optimizer.differential_evolution_minimize(
                self._objective,
                initial_population=self.initial_population_tensor,
                max_iterations=max_iterations,
                differential_weight=differential_weight,
                crossover_prob=crossover_prob,
                seed=seed
            )
        else:
            optimization_results = tfp.optimizer.differential_evolution_minimize(
                self._linearpartition_objective,
                initial_population=self.initial_population_tensor,
                max_iterations=max_iterations,
                differential_weight=differential_weight,
                crossover_prob=crossover_prob,
                seed=seed
            )
        elapsed_time = timeit.default_timer() - start
        print(f"Elapsed time: {timedelta(seconds=elapsed_time)}")
        self.optimization_results = optimization_results

        top_seq, top_score = self.get_optimized_sequences()[0]
        print(f"Best score after {max_iterations}: {top_score:.2f}")

    def get_optimized_sequences(self):
        """ Return a ranked list of optimized sequences and their associated scores (ascending order).
        Raises ValueError if called before optimize()

        :return: list of tuples (sequence, score), sorted by score (ascending)
        """
        if self.optimization_results is None:
            raise ValueError("No optimized sequences available. Run optimize() first")
        nseqs = self._convert_to_nseqs(self.optimization_results.final_population)
        scores = list(self.optimization_results.final_objective_values.numpy())
        sorted_nseqs = sorted(zip(nseqs, scores), key=lambda x: x[1])
        return sorted_nseqs

    def get_initial_sequences(self):
        """ Return a ranked list of initial population sequences (not the optimized ones!) and their associated scores
        (ascending order).

        :return: list of tuples (sequence, score), sorted by score (ascending)
        """
        nseqs = self._convert_to_nseqs(self.initial_population_tensor)
        scores = list(self._objective(self.initial_population_tensor).numpy())
        sorted_nseqs = sorted(zip(nseqs, scores), key=lambda x: x[1])
        return sorted_nseqs

    def _generate_initial_population(self):
        """
        Assembles initial population of random RNA sequences coding for given protein sequence.
        Codon usage reflects the fequencies from the chosen species.

        :return: initial population (Tensor objects to be used by the alg, not actual sequences)
        """
        initial_population = []
        for i in range(self.population_size):
            member = []
            for aa in self.protein_sequence:
                #  random.choices function: https://pynative.com/python-weighted-random-choices-with-probability/
                index_options = range(len(self.codon_map[aa]["codons"]))
                index = random.choices(index_options, weights=self.codon_map[aa]["weights"], k=1)
                member.extend(index)
            initial_population.append(member)

        # convert to TF object
        initial_population = tf.convert_to_tensor((initial_population), float)
        return initial_population

    def _convert_to_nseqs(self, population):
        """ Helper function to convert Tensor objects returned by the optimizer to RNA sequences.
        Aka Dillion's magic hack.

        :param population: population of sequences (tensor objects)
        :return: list of str, RNA sequences
        """
        # This is a hack. TF deals with continuous valued functions. We need discrete and finite.
        # Whatever values are assigned, make them ints and take the absolute value.
        population = np.absolute(np.array(population).astype(int))

        # Now we want to do something with the values. It's possible that some values exceed the
        # number of codons for the given position, so take the modulus. This is effectively a hashing
        # function. It's not mathematically rigorous, but it's good enough.
        # Finally, convert list of indices to the RNA sequence.
        get_seq = lambda se: ''.join([self.codon_map[res]['codons'][se[i] % self._get_nc(res)]
                                      for i, res in enumerate(self.protein_sequence)])
        n_seqs = [get_seq(se) for se in population]
        return n_seqs

    def _get_nc(self, aa):
        """ Helper function to get the number of codons for the given aminoacid

        :param aa: str, single-letter amino acid symbol
        :return: number of possible codons
        """
        return len(self.codon_map[aa]['codons'])

    def _objective(self, population):
        """Objective function for TF to minimize.

        NOTE: TF uses gradient descent to minimize continuous valued functions.
        The approach used here is not mathematically sound. It's a hack. But it gets the job done.

        :param population: population of sequences (tensor objects)
        :return: tf.float32, list of scores calculated for the population
        """
        # Map continuous valued tensor to RNA sequence
        n_seqs = self._convert_to_nseqs(population)
        # Use the imported scoring function to score all sequences.
        scores = [self.scoring_function(s) for s in n_seqs]
        # Return TF object
        return tf.cast(scores, float)

    def _linearpartition_objective(self, population, beam_size=100, linear_partition_home="/hpc/projects/upt/vx_codeRNA/software/LinearPartition/"):
        """Objective function for TF to minimize.
        NOTE: TF uses gradient descent to minimize continuous valued functions.
        The approach used here is not mathematically sound. It's a hack. But it gets the job done.
    
        Run LinearPartition and return MFE (Minimum Free Energy) value.
        See: https://github.com/LinearFold/LinearPartition

        :param sequence: nucleotide sequence
        :param beam_size: The beam size (default 100). Use 0 for infinite beam.
        :param linear_partition_home: path to linear partition installation, if not provided, the function will look for
            $LINEAR_PARTITION_HOME
        :return: tf.float32, list of scores calculated for the population
        """
    
        # Submit all processes
        processes = []
        for s in self._convert_to_nseqs(population):
            # Open temp file handle to stash output
            f_out = tempfile.TemporaryFile()
            f_err = tempfile.TemporaryFile()
            # Run executable -- will spawn child process on new thread
            seq = subprocess.Popen(["echo", s], stdout=subprocess.PIPE)
            p = subprocess.Popen([os.path.join(linear_partition_home, 'linearpartition'), "-V", "-p", "-b", str(beam_size)],
                          stdin=seq.stdout, stdout=f_out, stderr=f_err)
            # Stash output handles and read it when the process completes
            processes.append((p, f_err))
        
        scores = []
        # Handle output (synchronously) when processes complete
        for p, f in processes:
            # Wait for process to complete before handling output
            p.wait()
            # Rewind to the beginning of the file so it can be dumped
            f.seek(0)
            # Extract MFE results
            mfe = re.search(rb'Free Energy of Ensemble: (.*) kcal/mol', f.read()).group(1)
            # Convert to float and store in list -- convert to tf object later
            scores.append(float(mfe))
            # Close file handle
            f.close()
    
        # Return TF object
        return tf.cast(scores, float)

    def _get_codon_map(self):
        """ Helper function the get conveniently formatted dict with codon and codon frequency information per amino
        acid.

        :return: dict of dicts, maps amino acids to codons and codon frequencies in the given organism
        """
        freqs = get_codon_frequencies(self.organism)
        codon_table = get_codon_table(self.codon_table_name)
        codon_map = dict()
        for aa, codons in codon_table.items():
            codon_map[aa] = {"codons": codons,
                             "weights": [freqs[codon] for codon in codons]}
        return codon_map


if __name__ == "__main__":
    aa_seq = "MEDAKNIKKGPAPFYPLEDGTAGEQLHKAMKRYALVPGTIAFTDAHIEVNITYAEYFEMSVRLAEAMKRYGLNTNHRIVVCSENSLQFFMPVLGALFIGVAVAPA"

    get_u_content = lambda seq: RNA(seq).get_uridine_content()
    revo = EvolutionaryRNAOptimizer(protein_sequence=aa_seq, organism="human", scoring_function=get_u_content)
    revo = EvolutionaryRNAOptimizer(protein_sequence=aa_seq, organism="human", scoring_function="linearpartition", population_size=48)

    revo.optimize(max_iterations=100)
