from codeRNA.codon_data import cai_weights
from codeRNA.codon_data import codon_tables


def get_rare_codons_from_cai_weights(organism="human", codon_table_name="standard", cai_weight_threshold=0.5):

    codon_weights = cai_weights.get_codon_weights(organism=organism)

    codon_table = codon_tables.get_codon_table(codon_table_name=codon_table_name)

    rare_codons = []
    for aa, codon_options in codon_table.items():
        weighted_aa_codons = [(codon, codon_weights[codon]) for codon in codon_options]
        # sort by weights
        weighted_aa_codons = sorted(weighted_aa_codons, key=lambda x: x[1], reverse=True)
        # select only codons with weights below the specified threshold
        rare_aa_codons = [codon for codon, weight in weighted_aa_codons if weight <= cai_weight_threshold]
        # if all codons for the given amino acid have CAI weights >= cai_weight_threshold,
        # exclude the top codon (index 0 since the list has been sorted) from the rare codon list
        if len(rare_aa_codons) == len(codon_options):
            rare_aa_codons = rare_aa_codons[1:]
        rare_codons.extend(rare_aa_codons)

    return rare_codons


def get_complement_codons(codon_list):
    all_codons = [el for sublist in codon_tables.CODON_TABLES["standard"].values() for el in sublist]
    complement_codons = [el for el in all_codons if el not in codon_list]
    return complement_codons


def check_if_each_aa_represented_in_codon_list(codon_list, codon_table_name="standard"):
    codon_table = codon_tables.get_codon_table(codon_table_name=codon_table_name)

    for aa, aa_codon_options in codon_table.items():
        if not any([el in aa_codon_options for el in codon_list]):
            print("No codons for amino acid {} in the provided list of codons".format(aa))
            return False
    return True


if __name__ == "__main__":
    rare_codons = get_rare_codons_from_cai_weights(cai_weight_threshold=0.5)
    complement_codons = get_complement_codons(rare_codons)
    assert check_if_each_aa_represented_in_codon_list(complement_codons), "Not all amino acids represented"
    print(rare_codons)

