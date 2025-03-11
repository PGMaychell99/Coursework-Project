def DNA_Seq_ID(multi_seqfn: str = "dog_breeds.fa", mystery_seqfn: str = "mystery.fa") -> sequence [str]:
    from Bio import SeqIO
    from Bio import pairwise2

    