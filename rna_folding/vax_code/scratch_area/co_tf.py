# Available here: https://mygithub.gsk.com/dmf32881/quantum_codon_opt
from quantum_codon_opt.src.classical_ga import CodonOptimization
from quantum_codon_opt.src.scoring import SeqScorer
import tensorflow as tf, numpy as np
import tensorflow_probability as tfp
from Bio.Seq import Seq

# Spike protein
seq = '''MFVFLVLLPL VSSQCVNLTT RTQLPPAYTN SFTRGVYYPD KVFRSSVLHS TQDLFLPFFS
     NVTWFHAIHV SGTNGTKRFD NPVLPFNDGV YFASTEKSNI IRGWIFGTTL DSKTQSLLIV
     NNATNVVIKV CEFQFCNDPF LGVYYHKNNK SWMESEFRVY SSANNCTFEY VSQPFLMDLE
     GKQGNFKNLR EFVFKNIDGY FKIYSKHTPI NLVRDLPQGF SALEPLVDLP IGINITRFQT
     LLALHRSYLT PGDSSSGWTA GAAAYYVGYL QPRTFLLKYN ENGTITDAVD CALDPLSETK
     CTLKSFTVEK GIYQTSNFRV QPTESIVRFP NITNLCPFGE VFNATRFASV YAWNRKRISN
     CVADYSVLYN SASFSTFKCY GVSPTKLNDL CFTNVYADSF VIRGDEVRQI APGQTGKIAD
     YNYKLPDDFT GCVIAWNSNN LDSKVGGNYN YLYRLFRKSN LKPFERDIST EIYQAGSTPC
     NGVEGFNCYF PLQSYGFQPT NGVGYQPYRV VVLSFELLHA PATVCGPKKS TNLVKNKCVN
     FNFNGLTGTG VLTESNKKFL PFQQFGRDIA DTTDAVRDPQ TLEILDITPC SFGGVSVITP
     GTNTSNQVAV LYQDVNCTEV PVAIHADQLT PTWRVYSTGS NVFQTRAGCL IGAEHVNNSY
     ECDIPIGAGI CASYQTQTNS PRRARSVASQ SIIAYTMSLG AENSVAYSNN SIAIPTNFTI
     SVTTEILPVS MTKTSVDCTM YICGDSTECS NLLLQYGSFC TQLNRALTGI AVEQDKNTQE
     VFAQVKQIYK TPPIKDFGGF NFSQILPDPS KPSKRSFIED LLFNKVTLAD AGFIKQYGDC
     LGDIAARDLI CAQKFNGLTV LPPLLTDEMI AQYTSALLAG TITSGWTFGA GAALQIPFAM
     QMAYRFNGIG VTQNVLYENQ KLIANQFNSA IGKIQDSLSS TASALGKLQD VVNQNAQALN
     TLVKQLSSNF GAISSVLNDI LSRLDKVEAE VQIDRLITGR LQSLQTYVTQ QLIRAAEIRA
     SANLAATKMS ECVLGQSKRV DFCGKGYHLM SFPQSAPHGV VFLHVTYVPA QEKNFTTAPA
     ICHDGKAHFP REGVFVSNGT HWFVTQRNFY EPQIITTDNT FVSGNCDVVI GIVNNTVYDP
     LQPELDSFKE ELDKYFKNHT SPDVDLGDIS GINASVVNIQ KEIDRLNEVA KNLNESLIDL
     QELGKYEQYI KWPWYIWLGF IAGLIAIVMV TIMLCCMTSC CSCLKGCCSC GSCCKFDEDD
     SEPVLKGVKL HYT'''.replace(' ','').replace('\n','')

# Detect GPU's
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

# Run all of the preprocessing from the CodonOptimization class, but
# no need to run the Genetic Algorithm (GA)
co = CodonOptimization(seq,ntrials=5000,auto_exec=False)

# Pull out initial population generation by previous command and
# convert to TF object
initial_members = tf.convert_to_tensor(([_[1] for _ in co._get_initial_members()]),np.float32)

# Helper function to get number of possible codons for an amino acid
get_nc = lambda res: len(co.code_map[res]['codons'])

def convert_to_nseqs(members):
    # This is a hack. TF deals with continuous valued functions. We need discrete and finite.
    # So let's cheat. Whatever values are assigned, make them ints and take the absolute value.
    members = np.absolute(np.array(members).astype(int))
    
    # Now we want to do something with the values. It's possible that some values exceed the
    # number of codons for the given position, so take the modulus. This is effectively a hashing
    # function. It's not mathematically rigorous, but it's good enough.
    # Finally, convert list of indices to the RNA sequence.
    get_seq = lambda se: ''.join([co.code_map[res]['codons'][se[i] % get_nc(res)] for i, res in enumerate(seq)])
    n_seqs = [get_seq(se) for se in members]
    return n_seqs

def objective(members):
    '''
    Objective function for TF to minimize
    
    NOTE: TF uses gradient descent to minimize continuous valued functions.
    The approach used here is not mathematically sound. It's a hack. But
    it gets the job done. 
    
    FOR DEMONSTRATION PURPOSES ONLY.
    
    '''
    
    # Map continuous valued tensor to RNA sequence
    n_seqs = convert_to_nseqs(members)
    
    # Use the imported scoring function to score all sequences.
    scores = [SeqScorer(s).score for s in n_seqs]
    
    # Return TF object
    return tf.cast(scores, np.float32)

# More tricks. 
# Differential_weight: controls strength of mutations. We basically want to turn this off.
# Crossover_prob: set this low. Need to think more about why this helps.
optim_results = tfp.optimizer.differential_evolution_minimize(
    objective,
    initial_population=initial_members,
    max_iterations=500,
    differential_weight=0.01,
    crossover_prob=0.1,
)
# Translate "best" result back to protein sequence to verify it is valid
nseq = convert_to_nseqs(optim_results.final_population)[np.argmin(optim_results.final_objective_values)]
if seq != str(Seq(nseq).transcribe().translate()):
    print('craaaap')
print('TF with 500 iterations:',np.min(optim_results.final_objective_values))

if False:
    # Compare to GA with defaults
    co.execute()
    print('GA with 500 iterations:',co.score)


