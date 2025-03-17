from Bio import SeqIO
from Bio import pairwise2
multi_seqfn = "./data/project_dog_dna/dog_breeds.fa"
mystery_seqfn = "./data/project_dog_dna/mystery.fa"

def DNA_Seq_ID(multi_seqfn: str) -> dict:
    """ Reads in FASTA file and returns a dictionary of sequences.

        Description

        Parameters -
        Returns - 

        Example:
    """
    dog_sequences = {}

    with open(multi_seqfn, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            dog_sequences[record.id] = str(record.seq) #uses biopython to iterate over each line in the file and store the sequences in a dictionary
    return dog_sequences

def find_best_match(mystery_seq: str, sequences: dict) -> tuple:
    """ Find the closest sequence in the multi-sequence file to the mystery sequence.

        Takes the dictionary from FASTA file conversion dictionary

        Parameters - 
        Returns - 

        Example:
    """
    best_match = None
    best_score = -1
    best_alignment = None
    best_seq_id = None

    for seq_id, seq in sequences.items():
        score = pairwise2.align.globalxx(mystery_seq, seq, score_only = True)
        

        if score > best_score:
            best_score = score
            best_match = seq
            best_seq_id = seq_id
            best_alignment = pairwise2.align.globalxx(mystery_seq, seq)[0]

    return best_seq_id, best_match, best_score, best_alignment

def calculate_differences(seq1: str, seq2: str) -> list:
    """
    
    
    """
    alignments = pairwise2.align.globalxx(seq1, seq2)

    if not alignments:
        return None, None, []
    aligned_seq1, aligned_seq2, _, _, _ = alignments[0]

    differences = []
    for i, (base1, base2) in enumerate(zip(aligned_seq1, aligned_seq2)):
        if base1 != base2:
            differences.append((i, base1, base2))
    
    return aligned_seq1, aligned_seq2, differences

sequences = DNA_Seq_ID(multi_seqfn)
mystery_seq_dict = DNA_Seq_ID(mystery_seqfn)

if mystery_seq_dict:
    mystery_seq = list(mystery_seq_dict.values())[0]
else:
    raise ValueError("Mystery sequence file is empty or incorrectly formatted")

best_seq_id, best_match, best_score, best_alignment = find_best_match(mystery_seq, sequences)

print(f"Best Match: {best_seq_id}")
print(f"Alignment Score: {best_score}")
print(f"Aligned Mystery Sequence:\n{best_match}")
print(f"Aligned Best Match Sequence:\n{best_alignment}")
#print(f"Differences (Position, Mystery Base, Matched Base):")
#for diff in differences: #need to be defined within the for loop - not outside of it
    #print(diff)