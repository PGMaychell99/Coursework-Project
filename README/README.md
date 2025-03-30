DNA IDENTIFICATION SERVICE
This code is a DNA identification service that utilises the Biopython library to compute the following:
    - Compare a target sequence to multiple sequences
    - Identify the closest sequence to the target and any differences between them
    - Construct a phylogenetic tree

FEATURES
    - Parses DNA sequences from FASTA files.
    - Aligns sequences using pairwise alignment.
    - Calculates distances between DNA sequences.
    - Constructs and visualises phylogenetic trees to identify the mystery DNA's closest match.

REQUIREMENTS
    - Python 3.8+
    - Biopython
    - matplotlib

USAGE
Ensure the FASTA files are located in the ./data/project_dog_dna/ directory. Then run the script
    'python dna_sequence_comparison.py'

FUNCTIONS OVERVIEW
dna_seq_id(multi_seqfn: str) -> dict: Parses a FASTA file and returns a dictionary of sequence IDs and sequences.

find_best_match(mystery_seq: str, sequences: dict) -> tuple: Compares each sequence in the multi-sequence file to the mystery sequence and returns the ID and sequence of the best sequence match, the alignment score and the alignment between the two closest sequences.

calculate_differences(aligned_seq1: str, aligned_seq2: str) -> tuple: Calculates and returns the base differences between the mystery sequence and the closest match sequence.

build_phylogenetic_tree(sequences: dict) -> Phylo.BaseTree.Tree and plot_phylogenetic_tree(tree: Phylo.BaseTree.Tree): aligns all sequences, creates, plots and shows a phylogenetic tree based on the distance between each sequence.

OUTPUT
Printed results showing the closest breed match to the mystery DNA. A visual phylogenetic tree showing sequence relationships


