CODON_TABLES = {"standard":
                    {'A': ['GCA', 'GCC', 'GCG', 'GCU'],
                     'C': ['UGC', 'UGU'],
                     'D': ['GAC', 'GAU'],
                     'E': ['GAA', 'GAG'],
                     'F': ['UUC', 'UUU'],
                     'G': ['GGA', 'GGC', 'GGG', 'GGU'],
                     'H': ['CAC', 'CAU'],
                     'I': ['AUA', 'AUC', 'AUU'],
                     'K': ['AAA', 'AAG'],
                     'L': ['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'],
                     'M': ['AUG'],
                     'N': ['AAC', 'AAU'],
                     'P': ['CCA', 'CCC', 'CCG', 'CCU'],
                     'Q': ['CAA', 'CAG'],
                     'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU'],
                     '*': ['UAA', 'UAG', 'UGA'],
                     'S': ['AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU'],
                     'T': ['ACA', 'ACC', 'ACG', 'ACU'],
                     'V': ['GUA', 'GUC', 'GUG', 'GUU'],
                     'W': ['UGG'],
                     'Y': ['UAC', 'UAU']}
                }


def get_codon_table(codon_table_name="standard"):
    if codon_table_name not in CODON_TABLES:
        raise AttributeError('Invalid codon table name. Valid names: "{}"'.format(
            '", "'.join(CODON_TABLES.keys())
        ))
    codon_table = CODON_TABLES[codon_table_name]
    return codon_table
