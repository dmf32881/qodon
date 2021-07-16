CODON_WEIGHTS = {
    "human": dict(UUU=0.85185, UUC=1.0, UUA=0.2, UUG=0.325, UCU=0.79167, UCC=0.91667, UCA=0.625, UCG=0.20833,
                  UAU=0.78571, UAC=1.0, UAA=0.6383, UAG=0.51064, UGU=0.85185, UGC=1.0, UGA=1.0, UGG=1.0, CUU=0.325,
                  CUC=0.5, CUA=0.175, CUG=1.0, CCU=0.90625, CCC=1.0, CCA=0.875, CCG=0.34375, CAU=0.72414, CAC=1.0,
                  CAA=0.36986, CAG=1.0, CGU=0.38095, CGC=0.85714, CGA=0.52381, CGG=0.95238, AUU=0.76596, AUC=1.0,
                  AUA=0.3617, AUG=1.0, ACU=0.69444, ACC=1.0, ACA=0.77778, ACG=0.30556, AAU=0.88679, AAC=1.0,
                  AAA=0.75439, AAG=1.0, AGU=0.625, AGC=1.0, AGA=1.0, AGG=1.0, GUU=0.3913, GUC=0.52174, GUA=0.26087,
                  GUG=1.0, GCU=0.675, GCC=1.0, GCA=0.575, GCG=0.275, GAU=0.85185, GAC=1.0, GAA=0.72414, GAG=1.0,
                  GGU=0.47059, GGC=1.0, GGA=0.73529, GGG=0.73529),
    "hamster": dict(UUU=0.693, GUU=0.311, UCG=0.212, GCC=1.0, ACA=0.72, GGA=0.784, UUG=0.317, AGG=1.0, GCA=0.533,
                    UAC=1.0, AAU=0.688, AAG=1.0, ACC=1.0, CAG=1.0, CUG=1.0, UCC=0.963, UUC=1.0, GAG=1.0, ACU=0.619,
                    CUA=0.165, GGC=1.0, CUC=0.487, AGA=0.991, CGC=0.856, ACG=0.266, CGG=0.883, AGC=1.0, CGU=0.441,
                    CAA=0.307, CCU=0.908, AGU=0.619, GAC=1.0, UUA=0.129, UAA=0.583, GUG=1.0, AUU=0.647, GUA=0.219,
                    CAU=0.662, UAG=0.583, CAC=1.0, AUA=0.261, UGG=1.0, GGU=0.545, GUC=0.536, UAU=0.718, UGA=1.0,
                    GCG=0.213, CCC=1.0, AAA=0.562, CGA=0.568, GGG=0.751, GAU=0.691, CUU=0.302, UCU=0.799, UGU=0.746,
                    AUG=1.0, AUC=1.0, CCG=0.281, GAA=0.638, AAC=1.0, CCA=0.872, UCA=0.556, GCU=0.71, UGC=1.0)
}


def get_codon_weights(organism="human"):
    if organism not in CODON_WEIGHTS:
        raise AttributeError('Invalid organism name for CAI calculations. Valid names: "{}"'.format(
            '", "'.join(CODON_WEIGHTS.keys())
        ))
    codon_weights = CODON_WEIGHTS[organism]
    return codon_weights
