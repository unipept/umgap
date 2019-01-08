#!/bin/python3
"""
Experiments
"""

from collections import Counter
import matplotlib.pyplot as plt

DATA_DIR = "/home/rien/Development/Thesis/data/metacyc/data/"


PROTSEQ_FILE = 'protseq.fsa'
FASTA_HEADER_PREFIX = ">gnl|META|"
KMER_LENGTH = 9

def read_protseqs(filename):
    """
    Yields tuples of (protein, sequence)
    """
    with open(filename, 'r') as fasta:
        header = fasta.readline()
        while header:
            assert header.startswith(FASTA_HEADER_PREFIX)
            seq = fasta.readline()
            yield header[len(FASTA_HEADER_PREFIX):header.index(" ")], seq
            header = fasta.readline()

def protseq_kmers(protseq_generator, k):
    """
    Yields tuples of (kmer, protein)
    """
    for prot, seq in protseq_generator:
        yield from ((seq[i:i+k], prot) for i in range(len(seq) - k - 1))

def kmer_multiplicity(k):
    kmers_prots = set(protseq_kmers(read_protseqs(f"{DATA_DIR}/{PROTSEQ_FILE}"), k))
    kmers = [k[0] for k in kmers_prots]
    kmer_proteins_counts = [c[1] for c in Counter(kmers).most_common()]
    return Counter(kmer_proteins_counts).most_common()

def kmer_proteins_plot(k):
    countplot_data = kmer_multiplicity(k)
    figure, axis = plt.subplots()
    axis.set_xscale("log")
    plt.plot([c[1] for c in countplot_data], [c[0] for c in countplot_data], "o")
    plt.show()
    figure.savefig("kmer_proteins.pdf", bbox_inches="tight")

#kmer_proteins_plot(KMER_LENGTH)





