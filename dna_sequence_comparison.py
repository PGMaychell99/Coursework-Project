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
    
    return aligned_seq1, aligned_seq2, differences

#Test code using pytest

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
#print(f"Differences (Position, Mystery Base, Matched Base):")
#for diff in differences: #need to be defined within the for loop - not outside of it
    #print(diff)