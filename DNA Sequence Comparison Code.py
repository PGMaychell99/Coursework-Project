from Bio import SeqIO
from Bio import pairwise2

def DNA_Seq_ID(multi_seqfn: str = "dog_breeds.fa") -> sequence [str]:
    """ Reads in FASTA file and returns a dictionary of sequences.


        Example:"""
    sequences = {}
    with open(multi_seqfn, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences[record.id] = str(record.seq)
    return sequences

def find_best_match()

    