from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import os
import sys
import pytest

from dna_sequence_comparison import dna_seq_id, find_best_match, calculate_differences

def test_dna_seq_id(mock_fasta_files):
    multi_seqfn, _ = mock_fasta_files
    sequences = dna_seq_id(multi_seqfn)

    assert len(sequences) == 3 #check if the dictionary contains 3 sequences

    assert "seq1" in sequences
    assert sequences["seq1"] == "AGCTAGCTAG" #check if sequence "seq1" is in the dictionary and whether it matches the expected value

def test_find_best_match(mock_align, mock_fasta_files):
    multi_seqfn, mystery_seqfn = mock_fasta_files #get the multi-seq and mystery sequence files
    sequences = dna_seq_id(multi_seqfn) #store the multi-seq files in dictionary

    with open(mystery_seqfn, "r") as f: 
        mystery_seq = f.readlines()[1].strip() #take the sequence from the test mystery seq file

    best_seq_id, best_match, best_score, best_aligment = find_best_match(mystery_seq, sequences)

    assert best_seq_id == "seq1" #check best matching seq is "seq1"
    assert best_score == 10 #check that the alignment score is 10 (matches the mocked return value)

def test_calculate_differences():
    seq1 = "AGCTAGCTAG"
    seq2 = "TCGAAGCTAG"

    aligned_seq1, aligned_seq2, differences = calculate_differences(seq1, seq2) #call the function to find the aligned sequences and their differences

    assert len(differences) > 0 #check there is at least one difference
    assert differences[0] == (0, 'A', 'T') #check first differences occurs at position 0, where the bases are A(seq1) and T(seq2)

if __name__ == '__main__':
    pytest.main()