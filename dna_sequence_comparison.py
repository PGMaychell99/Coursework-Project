__all__ = [
    'calculate_differences',
    'dna_sequences_by_id',
    'find_best_match',
]


# -- IMPORTS --

# -- Standard libraries --
from pathlib import Path

# -- 3rd party libraries --
import Bio

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

# -- Internal libraries --


def dna_sequences_by_id(multi_seqfn: str) -> dict[str, str]:
    """:py:class:`str` : Reads in FASTA file and returns a dictionary of sequences.

        Reads in FASTA file and returns a dictionary of sequences.

        Parameters
        ----------
        multi_seqfn : str
            A multi-sequence FASTA string

        Returns
        -------
        dict
            A dictionary of sequences keyed by sequence ID.

        Examples
        --------
        >>> sequences = dna_sequences_by_id("./data/project_dog_dna/dog_breeds.fa")
        >>> assert isinstance(sequences, dict)
    """
    sequences = {}

    with open(multi_seqfn, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences[record.id] = str(record.seq) #uses biopython to iterate over each line in the file and store the sequences in a dictionary

    return sequences


def find_best_match(mystery_seq: str, sequences: dict[str, str]) -> tuple[str, str, float, Bio.Align.Alignment]:
    """:py:class:`tuple` : Find the closest sequence in the multi-sequence file to the mystery sequence.

        Takes the dictionary from FASTA file conversion dictionary.

        Parameters
        ----------
        mystery_seq : str
            The mystery sequence (target) to compare individual dog DNA
            sequences against.

        Returns
        -------
        tuple
            A tuple of information for the best sequence match, consisting
            of the best score, best match/seq., best sequence ID, and best
            alignment.

        Examples
        --------
    """
    aligner = PairwiseAligner()
    best_match = None
    best_score = -1
    best_alignment = None
    best_seq_id = None
    
    mystery_seq_object = Seq(mystery_seq)
    print(f'Mystery seq (target): {mystery_seq_object[:20]}...')
    print(f'Num. sequences to compare: {len(sequences)}')

    for i, (seq_id, seq) in enumerate(sequences.items(), start=1):
        seq_obj = Seq(seq)
        alignments = aligner.align(mystery_seq_object, seq_obj)
        score = alignments[0].score
        print(f'Comparing seq #{i}: {seq_id}:{seq[:20]}... (score={score}, best_score={best_score}, best_seq_id={best_seq_id})')

        if score > best_score:
            best_score = score
            best_match = seq
            best_seq_id = seq_id
            best_alignment = alignments[0]

    return best_seq_id, best_match, best_score, best_alignment

def calculate_differences(seq1: str, seq2: str) -> tuple[str, str, list]:
    """:py:class:`list` : Calculate differences between two sequences.
    
    Parameters
    ----------
    seq1 : str
        The first sequence

    seq2 : str
        The second sequence

    Returns
    -------
    tuplke
        A tuple of objects including alignments for sequences #1 and #2, and
        their differences.
    
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


if __name__ == '__main__':
    multi_seqfn = Path("./data/project_dog_dna/dog_breeds.fa").resolve()
    mystery_seqfn = Path("./data/project_dog_dna/mystery.fa").resolve()

    sequences = dna_sequences_by_id(multi_seqfn)
    mystery_seq_dict = dna_sequences_by_id(mystery_seqfn)

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