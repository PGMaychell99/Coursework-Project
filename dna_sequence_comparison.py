from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import os
import pytest

multi_seqfn = "./data/project_dog_dna/dog_breeds.fa"
mystery_seqfn = "./data/project_dog_dna/mystery.fa"

def dna_seq_id(multi_seqfn: str) -> dict:
    """ Reads in FASTA file and returns a dictionary of sequences.

        Description

        Parameters -
        Returns - 

        Example:
    """
    sequences = {}

    with open(multi_seqfn, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences[record.id] = str(record.seq) #uses biopython to iterate over each line in the file and store the sequences in a dictionary
    return sequences

def find_best_match(mystery_seq: str, sequences: dict) -> tuple:
    """ Find the closest sequence in the multi-sequence file to the mystery sequence.

        Takes the dictionary from FASTA file conversion dictionary

        Parameters - 
        Returns - 

        Example:
    """
    aligner = PairwiseAligner()
    best_match = None
    best_score = -1
    best_alignment = None
    best_seq_id = None
    
    mystery_seq_object = Seq(mystery_seq)

    for i, (seq_id, seq) in enumerate(sequences.items(), start=1):
        seq_obj = Seq(seq)
        alignments = aligner.align(mystery_seq_object, seq_obj)
        score = alignments[0].score

        if score > best_score:
            best_score = score
            best_match = seq
            best_seq_id = seq_id
            best_alignment = alignments[0]

    return best_seq_id, best_match, best_score, best_alignment

def calculate_differences(seq1: str, seq2: str) -> tuple:
    """
    
    
    """
    aligner = PairwiseAligner()
    alignments = aligner.align(seq1, seq2)

    if not alignments:
        return None, None, []
    
    aligned_seq1, aligned_seq2, _, _, _ = alignments[0]

    differences = []
    for i, (base1, base2) in enumerate(zip(aligned_seq1, aligned_seq2)):
        if base1 != base2:
            differences.append((i, base1, base2))
    
    print(f"Differences (Position, Mystery Base, Matched Base):")
    for diff in differences: #need to be defined within the for loop - not outside of it
        print(diff)

    return aligned_seq1, aligned_seq2, differences

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

    print(f"Best Match: {best_seq_id}")
    print(f"Alignment Score: {best_score}")
    print(f"Aligned Mystery Sequence:\n{best_match}")
    print(f"Aligned Best Match Sequence:\n{best_alignment}")
