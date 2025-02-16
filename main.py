# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    def get_species_name(header: str) -> str:
        parts = header[1:].split('_')  # Remove '>' and split by underscores
        return f"{parts[0]} {parts[1]}"  # Join first two parts (Genus species)
    
    species = [
        (get_species_name(gg_header), gg_seq),
        (get_species_name(mm_header), mm_seq),
        (get_species_name(br_header), br_seq),
        (get_species_name(tt_header), tt_seq)
    ]

    aligner = NeedlemanWunsch(
        sub_matrix_file="substitution_matrices/BLOSUM62.mat",
        gap_open=-10,
        gap_extend=-1
    )

    scores = []
    for name, seq in species:
        score, _, _ = aligner.align(hs_seq, seq)
        scores.append((name, score))

    sorted_species = sorted(scores, key=lambda x: x[1], reverse=True)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
     # Print species in order of similarity
    print("Species ordered by similarity to human BRD2:")
    for name, score in sorted_species:
        print(name)

    # Print alignment scores
    print("\nAlignment scores:")
    for name, score in sorted_species:
        print(f"{name}: {score}")

if __name__ == "__main__":
    main()
