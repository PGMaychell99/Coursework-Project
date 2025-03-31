from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio import Phylo
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

multi_seqfn = "./data/project_dog_dna/dog_breeds.fa"
mystery_seqfn = "./data/project_dog_dna/mystery.fa"

def dna_seq_id(multi_seqfn: str) -> dict:
    """ Reads in FASTA file containing multiple sequences and returns a dictionary of sequences.

        This function reads in the FASTA file containing multiple sequences (multi_seqfn), iterates over each item in the file, assigns the ID to key and the sequence as the value.

        Parameters - str(FASTA file containing sequence)
        Returns - dict(dictionary of key(ID)-value(sequence) pairs)
    """
    sequences = {}

    with open(multi_seqfn, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences[record.id] = str(record.seq) #uses biopython to iterate over each line in the file and store the sequences in a dictionary
    return sequences

def find_best_match(mystery_seq: str, sequences: dict) -> tuple:
    """ Find the closest sequence in the multi-sequence file to the mystery sequence.

        Takes the dictionary from FASTA file conversion dictionary in dna_seq_id, 

        Parameters - str(mystery sequence FASTA file containing the target dog breed sequence) and dict(multiple sequences stored in dictionary from previous function)
        Returns - tuple(returns the following elements:
                        1. best_seq_id(str): The sequence ID of the best match
                        2. best_match(str): The sequence that is the best match
                        3. best_score(float): The alignment score of the best match
                        4. best_alignment(Alignment object): The alignment object that represents the best alignment between the mystery sequence and the best match)

        Example:
        >>> sequences = {
                    "seq1": "AGCTAGCTA",
                    "seq2": "AGCTCGTAA",
                    "seq3": "CGATCGTA"}
        >>> 'seq1', 'AGCTAGCTA', 90.0, <alignment object>
    """
    aligner = PairwiseAligner()
    best_match = None
    best_score = -1
    best_alignment = None
    best_seq_id = None #initiate the variables for storing the best matches from the sequences
    
    mystery_seq_object = Seq(mystery_seq) #convert the mystery sequence string into a seq object

    for i, (seq_id, seq) in enumerate(sequences.items(), start=1): #iterate over the sequences in the multi sequence file starting from pos 1
        seq_obj = Seq(seq) #convert the current sequence into a sequence object
        alignments = aligner.align(mystery_seq_object, seq_obj) #align mystery seq and current seq
        score = alignments[0].score #store the alignment score

        if score > best_score: #if the current score is better than the previous score, store the information of the current sequence overwriting the previous
            best_score = score
            best_match = seq
            best_seq_id = seq_id
            best_alignment = alignments[0]

    return best_seq_id, best_match, best_score, best_alignment

def calculate_differences(aligned_seq1: str, aligned_seq2: str) -> tuple:
    """ Calculate the differneces between two aligned DNA sequences.

        This function aligns two sequences using pairwise sequence alignment and compares their aligned versions.
        It returns the differences between the sequences, including the positions and the bases that differ.

        Parameters - str(two sequences for comparison)

        Returns - tuple(returns the following elements:
                        1. aligned_seq1(str): The aligned version of seq1
                        2. aligned_seq2(str): The aligned version of seq2
                        3. differences(list): A list of tuples, each representing a difference between the two aligned sequences containing:
                            - position(int): The position in the aligned sequence where the difference occurs
                            - base1(str): The base in seq1 at the difference position
                            - base2(str): The base in seq2 at the difference position)

        Example"
        >>> seq1 = "AGCTAGCTA"
            seq2 = "AGCTCGTAA"

            'AGCTAGCTA', 'AGCTCGTAA', [(5,'A','C'), (7,'T','G'),(8,'A','T')]
    """

    differences = [] #initialise empty list to store the differences in
    for i, (base1, base2) in enumerate(zip(aligned_seq1, aligned_seq2)): #check if the base pairs at the current index are matched, if not then store in the differences list
        if base1 != base2 and base1 != '-' and base2 != '-': #check there are no gaps interferring ('-')
            differences.append((i, base1, base2))
    
    print(f"Differences (Position, Mystery Base, Matched Base):")
    for diff in differences: #print each difference between the sequences
        print(diff)

    return aligned_seq1, aligned_seq2, differences

def build_phylogenetic_tree(sequences: dict) -> Phylo.BaseTree.Tree:
    """ Build a phylogenetic tree from a dictionary of sequences.

        Computes a distance matrix between the sequences (based on identity) and the constructs a phylogenetic tree using the Neighnour-Joining method.
        It uses PairwiseAligner for pairwise sequence alignment and DistanceCalculator to generate the distance matrix.

        Parameters - dict(dictionary containing multiple sequences and their identities)
        Returns - Phylo.BaseTree.Tree()
    """
    aligner = PairwiseAligner() #initialise the PairwiseAligner
    aligned_sequences = [] #empty list to store all the pairwise alignments
    sequence_ids = list(sequences.keys()) #list for storing the IDs
    for i in range(len(sequences)): #perform pairwise alignmnet between all sequences
        for j in range(i+1, len(sequences)):
            seq1, seq2 = sequences[sequence_ids[i]], sequences[sequence_ids[j]]
            alignments = aligner.align(Seq(seq1), Seq(seq2))
            aligned_sequences.append(alignments[0])
        
    calculator = DistanceCalculator('identity') #create a DistanceCalculator object and multi-seq align to get sequence distances
    distance_matrix = calculator.get_distance(aligned_sequences)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix) #Generate tree using distance constructor and Neighbour-Joining

    return tree

def plot_phylogenetic_tree(tree: Phylo.BaseTree.Tree):
    """ Plot the phylogenetic tree using matplotlib.

        Parameters - tree(A Phylo.BaseTree.Tree object representing the phylogenetic tree)
    """
    fig = plt.figure(figsize=(10, 10))#define the size of the plot
    ax = fig.add_subplot(111)
    Phylo.draw(tree, do_show=False, axes=ax) #draw the tree on subplot (ax) so that it shows in plt.show()
    plt.show()

if __name__ == '__main__':
    multi_seqfn = "./data/project_dog_dna/dog_breeds.fa"
    mystery_seqfn = "./data/project_dog_dna/mystery.fa"

    sequences = dna_seq_id(multi_seqfn)
    mystery_seq_dict = dna_seq_id(mystery_seqfn)

    if mystery_seq_dict:
        mystery_seq = list(mystery_seq_dict.values())[0]
    else:
        raise ValueError("Mystery sequence file is empty or incorrectly formatted")


    best_seq_id, best_match, best_score, best_alignment = find_best_match(mystery_seq, sequences)

    aligned_seq1, aligned_seq2, differences = calculate_differences(best_match, mystery_seq)

    print(f"Best Match: {best_seq_id}")
    print(f"Alignment Score: {best_score}")
    print(f"Aligned Mystery Sequence:\n{best_match}")
    print(f"Aligned Best Match Sequence:\n{best_alignment}")
    
    tree = build_phylogenetic_tree(sequences)
    plot_phylogenetic_tree(tree)
