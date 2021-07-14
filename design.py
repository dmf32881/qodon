from qodon.src.classical_ga import CodonOptimization
from rna_folding.rna_fold import RNAFold
import tensorflow as tf, numpy as np
import tensorflow_probability as tfp
from Bio.Seq import Seq
from Bio import SeqIO
import warnings


class QuDesign(object):
    
    def __init__(self,seq,codon_opt_it=100,rna_fold_it=10000):
        self.seq = seq
        self.codon_opt_it = codon_opt_it
        self.rna_fold_it = rna_fold_it
        self.initial_members = self._get_initial_pop()
        
        self._validate()
        self.execute()
    
    def execute(self):
        '''
        Main execution. Run tensorflow optimizer for codon optimization. Objective
        function computes RNA structure with D-Wave's SA algorithm.
        
        ''' 
        
        # Differential_weight: controls strength of mutations. We basically want to turn this off.
        # Crossover_prob: set this low. Need to think more about why this helps.
        optim_results = tfp.optimizer.differential_evolution_minimize(
            self._objective,
            initial_population=self.initial_members,
            max_iterations=self.codon_opt_it,
            differential_weight=0.01,
            crossover_prob=0.1,
        )
        
        # Assign results as class attributes
        self. nseq = self._convert_to_nseqs(optim_results.final_population)[np.argmin(optim_results.final_objective_values)]
        self.mfe = np.min(optim_results.final_objective_values)
    
    def _get_initial_pop(self):
        '''
        Re-use code from the qodon package to compute an initial population.
        
        '''
        # Run all of the preprocessing from the CodonOptimization class, but
        # no need to run the Genetic Algorithm (GA)
        co = CodonOptimization(seq,lazy=True)
        
        # Store the "codon map". This contains information about Codon Usage Bias.
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
        '''
        Extract number of possible codons for each amino acid
        
        '''
        return len(self.code_map[res]['codons'])

    def _tf_fold(self, nseq):
        '''
        Compute Minimum Free Energy (MFE) of RNA fold.
        
        '''
        rna_ss = RNAFold(nseq, min_stem_len=4, min_loop_len=4)
        results = rna_ss.compute_dwave_sa(sweeps=self.rna_fold_it)
        return results.first.energy
    
    def _convert_to_nseqs(self, members):
        '''
        Continuous --> discrete transformation
        
        Doesn't make mathematical sense but it works.
        
        '''
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
    
    def _validate(self):
        '''
        Validate user input.
        
        '''
        
        if not isinstance(self.seq,str):
            
            raise TypeError('''
            Input protein sequence must be a string! User provided
            input with type {}
            
            '''.format(type(self.seq)))
        
        self.seq = self.seq.upper()
        
        aas = 'ACDEFGHIKLMNPQRSTVWY'
        
        if not any(_ not in aas for _ in self.seq):
            print('Not a valid input sequence!')
        
        if set(self.seq).issubset(set('GCAU')):
            warnings.warn("Input protein sequence looks like an RNA sequence!")
        
        if set(self.seq).issubset(set('GCAT')):
            warnings.warn("Input protein sequence looks like an DNA sequence!")
        
        if not isinstance(self.codon_opt_it, int):
            raise TypeError('''
            codon_opt_it must be a positive integer!
            
            ''')
        
        if self.codon_opt_it < 1:
            raise ValueError('''
            codon_opt_it must be at least 1!
            
            ''')
        
        if not isinstance(self.rna_fold_it, int):
            raise TypeError('''
            rna_fold_it must be a positive integer!
            
            ''')
        
        if self.rna_fold_it < 1:
            raise ValueError('''
            rna_fold_it must be at least 1!
            
            ''')


if __name__ == "__main__":
    
    seq = str(SeqIO.read('examples/spike_trim.fasta','fasta').seq)
    exe = QuDesign(seq)