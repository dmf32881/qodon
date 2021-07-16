"""
Codon usage counts/frequencies for scoring RNA sequences.
Codon usage frequencies are calculated for degenerate codons and as a whole.


Inspiration: https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-00968-8
And: https://github.com/JeffreyJLi/codon_optimization_analysis/

Currently only the tables for humans and E. coli is implemented.

Other tables can be found below by searching the full species name:
https://hive.biochemistry.gwu.edu/dna.cgi?cmd=tissue_codon_usage&id=586358&mode=cocoputs
"""
from codeRNA.codon_data.codon_tables import get_codon_table

HIVE_COUNTS = {
    "ecoli":{   "UUU":624098963, "UUC":451745454, "UUA":387996616, "UUG":374025846,
                "CUU":320096412, "CUC":304601370, "CUA":111020619, "CUG":1435663922,
                "AUU":836514697, "AUC":682532156, "AUA":143061192, "AUG":766478081,
                "GUU":512084887, "GUC":419381351, "GUA":306079034, "GUG":718217392,
                "UCU":243515357, "UCC":248216212, "UCA":218378482, "UCG":247653765,
                "CCU":204217118, "CCC":158427504, "CCA":238116553, "CCG":628802910, 
                "ACU":253980275, "ACC":635377411, "ACA":219450527, "ACG":407126465, 
                "GCU":434095669, "GCC":707498214, "GCA":575393774, "GCG":906462962, 
                "UAU":457776532, "UAC":340385225, "UAA":61281837, "UAG":9032292,
                "CAU":358918732, "CAC":265759641, "CAA":420439371, "CAG":811838472, 
                "AAU":512409788, "AAC":599800325, "AAA":943480826, "AAG":302630272, 
                "GAU":897451396, "GAC":532087064, "GAA":1092099062, "GAG":506595622, 
                "UGU":150316650, "UGC":183384315, "UGA":36929633, "UGG":428492249,
                "CGU":571197974, "CGC":597341472, "CGA":107634828, "CGG":167489061, 
                "AGU":256466509, "AGC":445436616, "AGA":74512491, "AGG":47440847, 
                "GGU":677829496, "GGC":795200802, "GGA":240152498, "GGG":317247601 
            },
    "human":{   "UUU":1345484, "UCU":1325110, "UAU":949353, "UGU":818211,
                "UUC":1367806, "UCC":1350473, "UAC":1052006, "UGC":846046,
                "UUA":689145, "UCA":1108355, "UAA":34502, "UGA":62824,
                "UUG":1052279, "UCG":316763, "UAG":27354, "UGG":910905,
                "CUU":1106851, "CCU":1496300, "CAU":929989, "CGU":356134,
                "CUC":1391429, "CCC": 1473428, "CAC":1142446, "CGC":680327,
                "CUA":584767, "CCA": 1467896, "CAA":1106752, "CGA":501278,
                "CUG":2813643, "CCG": 480590, "CAG":2756517, "CGG":829391,
                "AUU":1295762, "ACU": 1117630, "AAU":1446611, "AGU":1098600,
                "AUC":1459807, "ACC": 1388473, "AAC":1427415, "AGC":1536075,
                "AUA":637556, "ACA": 1294206, "AAA":2169279, "AGA":1046642,
                "AUG":1677888, "ACG": 436790, "AAG":2483951, "AGG":951166,
                "GUU":922234, "GCU": 1475932, "GAU":1882177, "GGU":845041,
                "GUC":1050009, "GCC": 2003500, "GAC":1894516, "GGC":1541147,
                "GUA":602084, "GCA": 1335414, "GAA":2644966, "GGA":1343257,
                "GUG":2011761, "GCG": 464292, "GAG":3078049, "GGG":1192278
            },
    }

# these are the degenerate codon frequencies.
# e.g. how frequently is each threonine codon used
HIVE_FREQUENCIES = {x:{} for x in HIVE_COUNTS.keys()}

# these are overall codon frequencies,
HIVE_FREQUENCIES_OVERALL = {x:{} for x in HIVE_COUNTS.keys()}

# copied from supplement of above paper
HIVE_FREQUENCIES['ecoli_pub'] = {'UUU': 58.0, 'UUC': 42.0, 'UUA': 13.2, 'UUG': 12.7,
                                 'CUU': 10.9, 'CUC': 10.4, 'CUA': 3.7, 'CUG': 49.2,
                                 'AUU': 50.6, 'AUC': 41.2, 'AUA': 8.2, 'AUG': 100.0,
                                 'GUU': 26.2, 'GUC': 21.4, 'GUA': 15.6, 'GUG': 36.8,
                                 'UCU': 14.7, 'UCC': 15.0, 'UCA': 13.0, 'UCG': 14.9,
                                 'CCU': 16.5, 'CCC': 12.7, 'CCA': 19.3, 'CCG': 51.5,
                                 'ACU': 16.7, 'ACC': 42.3, 'ACA': 14.2, 'ACG': 26.8,
                                 'GCU': 16.5, 'GCC': 27.0, 'GCA': 21.8, 'GCG': 34.7,
                                 'UAU': 57.4, 'UAC': 42.6, 'UAA': 61.5, 'UAG': 7.5,
                                 'CAU': 57.6, 'CAC': 42.4, 'CAA': 34.1, 'CAG': 65.9,
                                 'AAU': 46.0, 'AAC': 54.0, 'AAA': 76.0, 'AAG': 24.0,
                                 'GAU': 62.8, 'GAC': 37.2, 'GAA': 68.4, 'GAG': 31.6,
                                 'UGU': 45.1, 'UGC': 54.9, 'UGA': 31.0, 'UGG': 100.0,
                                 'CGU': 37.2, 'CGC': 38.7, 'CGA': 6.7, 'CGG': 10.3,
                                 'AGU': 15.4, 'AGC': 27.0, 'AGA': 4.5, 'AGG': 2.7,
                                 'GGU': 33.6, 'GGC': 39.3, 'GGA': 11.6, 'GGG': 15.5}
# divide by 100 to get frequencies
for k in HIVE_FREQUENCIES['ecoli_pub'].keys():
    HIVE_FREQUENCIES['ecoli_pub'][k] = HIVE_FREQUENCIES['ecoli_pub'][k]/100 

def get_codon_counts(organism="human"):

    valid_organisms = HIVE_COUNTS.keys()
  
    if organism not in valid_organisms:
        raise AttributeError(f"Invalid organism name for HIVE counts. Valid names: {', '.join(valid_organisms)}")
  
    return HIVE_COUNTS[organism]

def get_codon_frequencies(organism="human"):
    # these frequencies are the usage rates for each codon per amino acid
    # e.g. AUG is used 100% for methionine
    valid_organisms = HIVE_FREQUENCIES.keys()
  
    if organism not in valid_organisms:
        raise AttributeError(f"Invalid organism name for HIVE frequencies. Valid names: {', '.join(valid_organisms)}")

    if len(HIVE_FREQUENCIES[organism]) < 1:
        # if not calculated get frequency for each codon for each AA
        codon_table = get_codon_table("standard")
        for aa in codon_table.keys():
            aa_total_count = sum(HIVE_COUNTS[organism][x] for x in codon_table[aa])
            for codon in codon_table[aa]:
                HIVE_FREQUENCIES[organism][codon] = HIVE_COUNTS[organism][codon]/aa_total_count

    return HIVE_FREQUENCIES[organism]

def get_overall_codon_frequencies(organism="human"):
    # these are the overall frequencies
    # e.g. for ecoli methionine is use 26.2/1000
    if len(HIVE_FREQUENCIES_OVERALL[organism]) < 1:
    # populate
        for k in HIVE_FREQUENCIES_OVERALL.keys():
            total_counts = sum(HIVE_COUNTS[k].values())
            HIVE_FREQUENCIES_OVERALL[k] = { k2:v/total_counts for k2, v in HIVE_COUNTS[k].items() }

    return HIVE_FREQUENCIES_OVERALL[organism]