from Bio import SeqIO
from Bio import pairwise2

def DNA_Seq_ID(multi_seqfn: str = "dog_breeds.fa") -> sequence [str]:
    """ Reads in FASTA file and returns a dictionary of sequences.


        Example:"""
    sequences = {}
    with open(multi_seqfn, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences[record.id] = str(record.seq) #uses biopython to iterate over each line in the file and store the sequences in a dictionary
    return sequences

def find_best_match(mystery_seqfn: str = "mystery.fa") -> sequence [str]:
    """ Find the closest sequence in the multi-sequence file to the mystery sequence.

        Takes the dictionary from FASTA file conversion dictionary
    """
    best_match = None
    best_score = -1
    best_alignment = None
    best_seq_if = None

    for seq_id, seq in sequences.items():
        alignments = pairwise2.align.globalxx(mystery_seqfn, seq, score_only = True)
        score = alignments

        if score > best_score:
            best_score = score
            best_match = seq
            best_seq_id = seq_id

    return best_seq_id, best_match, best_score

def calculate_differences(seq1[str], seq2[str]) -> str:
    alignments = pairwise2.align.globalxx(seq1, seq2)
    aligned_seq1, aligned_seq2, _, _, _ = alignments[0]

    differences = []
    for i, (base1, base2) in enumerate(zip(aligned_seq1, aligned_seq2)):
        if base1 != base2:
            differences.append((i, base1, base2))
    
    return aligned_seq1, aligned_seq2, differences

sequences = read_fasta()
mystery_seq_dict = read_fasta()
mystery_seq = list(mystery_seq_dict.values())[0]

best_seq_id, best_match, best_score = find_best_match()

print(f"Best Match: {best_seq_id}")
print(f"Alignment Score: {best_score}")
print(f"Aligned Mystery Sequence:\n{aligned_mystery}")
print(f"Aligned Best Match Sequence:\n{aligned_best}")
print(f"Differences (Position, Mystery Base, Matched Base):")
for diff in differences:
    print(diff)