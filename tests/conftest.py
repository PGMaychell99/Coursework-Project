import os 
import pytest

def mock_fasta_files():
    mock_fasta = ">seq1\nAGCTAGCTAG\n>seq2\nTCGAAGCTAG\n>seq3\nAGCTAGCTAA/n"
    multi_seqfn = "test_dog_breeds.fa"
    mystery_seqfn = "test_mystery.fa"

    with open(multi_seqfn, "w") as f:
        f.write(mock_fasta)

    mock_mystery = "AGCTAGCTAG"
    with open(mystery_seqfn, "w") as f:
        f.write(f">mystery\n{mock_mystery}")

    yield multi_seqfn, mystery_seqfn

    os.remove(multi_seqfn)
    os.remove(mystery_seqfn)