import pytest
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

from "DNA Sequence Comparison Code" import dna_seq_id, find_best_match, calculate_differences

def mock_fasta_files():
    mock_fasta = ">seq1\nAGCTAGCTAG\n>seq2\nTCGAAGCTAG\n>seq3\nAGCTAGCTAA/n"
    multi_seqfn = "test_dog_breeds.fa"
    mystery_seqfn = "test_mystery.fa"

    with open(multi_seqfn, "w") as f:
        f.write(mock_fasta)

    mock_mystery = "AGCTAGCTAG"
    with open(mystery_seqfn, "W") as f:
        f.write(f">mystery\n{mock_mystery}")

    yield multi_seqfn, mystery_seqfn

    os.remove(multi_seqfn)
    os.remove(mystery_seqfn)

def test_dna_seq_id(mock_fasta_files):
    multi_seqfn, _ = mock_fasta_files
    sequences = dna_seq_id(multi_seqfn)

    assert len(sequences) == 3 #check if the dictionary contains 3 sequences

    assert "seq1" in sequences
    assert sequences["seq1"] == "AGCTAGCTAG" #check if sequence "seq1" is in the dictionary and whether it matches the expected value

def test_find_best_match(mock_align, mock_fasta_files):
    