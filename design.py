# Available here: https://mygithub.gsk.com/dmf32881/quantum_codon_opt
from qodon.src.classical_ga import CodonOptimization
from rna_folding.rna_fold import RNAFold
import tensorflow as tf, numpy as np
import tensorflow_probability as tfp
from Bio.Seq import Seq
from Bio import SeqIO


class QuDesign(object):
    
    def __init__(self,seq):
        self.seq = seq
        self.initial_members = self._get_initial_pop()
        self.execute()
    
    def execute(self):
        # More tricks. 
        # Differential_weight: controls strength of mutations. We basically want to turn this off.
        # Crossover_prob: set this low. Need to think more about why this helps.
        optim_results = tfp.optimizer.differential_evolution_minimize(
            self._objective,
            initial_population=self.initial_members,
            max_iterations=10,
            differential_weight=0.01,
            crossover_prob=0.1,
        )
        # Translate "best" result back to protein sequence to verify it is valid
        nseq = self._convert_to_nseqs(optim_results.final_population)[np.argmin(optim_results.final_objective_values)]
        print('TF with 100 iterations:',np.min(optim_results.final_objective_values))
    
    def _get_initial_pop(self):

        # Run all of the preprocessing from the CodonOptimization class, but
        # no need to run the Genetic Algorithm (GA)
        co = CodonOptimization(seq,lazy=True)
        
        self.code_map = co.code_map
        
        # Pull out initial population generation by previous command and
        # convert to TF object
        initial_members = tf.convert_to_tensor(([_[1] for _ in co.population]),np.float32)
        
        return initial_members
    
    def _objective(self,members):
        '''
        Objective function for TF to minimize
        
        NOTE: TF uses gradient descent to minimize continuous valued functions.
        The approach used here is not mathematically sound. It's a hack. But
        it gets the job done. 
        
        '''
        
        # Map continuous valued tensor to RNA sequence
        n_seqs = self._convert_to_nseqs(members)
        
        # Use the imported scoring function to score all sequences.
        scores = [self._tf_fold(s) for s in n_seqs]
        
        # Return TF object
        return tf.cast(scores, np.float32)
    
    # Helper function to get number of possible codons for an amino acid
    def _get_nc(self,res):
        return len(self.code_map[res]['codons'])

    def _tf_fold(self, nseq):
        rna_ss = RNAFold(nseq, min_stem_len=4, min_loop_len=4)
        results = rna_ss.compute_dwave_sa()
        return results.first.energy
    
    def _convert_to_nseqs(self, members):
        # This is a hack. TF deals with continuous valued functions. We need discrete and finite.
        # So let's cheat. Whatever values are assigned, make them ints and take the absolute value.
        members = np.absolute(np.array(members).astype(int))
        
        # Now we want to do something with the values. It's possible that some values exceed the
        # number of codons for the given position, so take the modulus. This is effectively a hashing
        # function. It's not mathematically rigorous, but it's good enough.
        # Finally, convert list of indices to the RNA sequence.
        get_seq = lambda se: ''.join([self.code_map[res]['codons'][se[i] % self._get_nc(res)] for i, res in enumerate(seq)])
        n_seqs = [get_seq(se) for se in members]
        return n_seqs


if __name__ == "__main__":
    
    seq = str(SeqIO.read('spike_trim.fasta','fasta').seq)
    exe = QuDesign(seq)