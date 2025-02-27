from collections import Counter
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import math
import time
import sys
import os


def calculate_sequential(sequence):
    counts = Counter(sequence)
    gc_content = (counts['G'] + counts['C']) / len(sequence) * 100
    return gc_content


def calculate_gc_content(chunk):
    count = Counter(chunk)
    return count['G'] + count['C']


def calculate_parallel(sequence, num_workers=5, use_multiprocessing=True):
    chunk_size = math.ceil(len(sequence) / num_workers)
    chunks = [sequence[i:i + chunk_size]
              for i in range(0, len(sequence), chunk_size)]

    executor_class = ProcessPoolExecutor if use_multiprocessing else ThreadPoolExecutor

    with executor_class(max_workers=num_workers) as executor:
        gc_counts = list(executor.map(calculate_gc_content, chunks))

    total_gc = sum(gc_counts)
    total_length = len(sequence)
    gc_content = (total_gc / total_length) * 100

    return gc_content


def read_sequence_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            sequence = file.read().replace('\n', '').replace(
                ' ', '').upper()  # formats in case it isn't done properly
            return sequence
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        sys.exit(1)


# measures the runtime
def measure_runtime(function, *args, **kwargs):
    start_time = time.time()
    result = function(*args, **kwargs)
    elapsed = time.time() - start_time
    return result, elapsed

# returns percent of each nucleotide (A, T, G, C) in the DNA sequence
def nucleotide_composition(sequence):
    dna_len = len(sequence)
    if dna_len == 0:
        return {"A": 0, "T": 0, "G": 0, "C":0}
    
    from collections import Counter

    counts = Counter(sequence)

    composition = {
        "A": (counts['A'] / dna_len) * 100,
        "T": (counts['T'] / dna_len) * 100,
        "G": (counts['G'] / dna_len) * 100,
        "C": (counts['C'] / dna_len) * 100
    }

    return composition

# translates dna sequence into amino acid sequence based on codon table.
def dna_protein_translation(sequence):

    # standard gemnetic code 
    codon_table = {
        # Phenylalanine
        "TTT": "F", "TTC": "F",

        # Leucine
        "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",

        # Isoleucine
        "ATT": "I", "ATC": "I", "ATA": "I",

        # Methionine (start codon)
        "ATG": "M",

        # Valine
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",

        # Serine
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",

        # Proline
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",

        # Threonine
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",

        # Alanine
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",

        # Tyrosine
        "TAT": "Y", "TAC": "Y",

        # Histidine
        "CAT": "H", "CAC": "H",

        # Glutamine
        "CAA": "Q", "CAG": "Q",

        # Asparagine
        "AAT": "N", "AAC": "N",

        # Lysine
        "AAA": "K", "AAG": "K",

        # Aspartic Acid
        "GAT": "D", "GAC": "D",

        # Glutamic Acid
        "GAA": "E", "GAG": "E",

        # Cysteine
        "TGT": "C", "TGC": "C",

        # Tryptophan
        "TGG": "W",

        # Arginine
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",

        # Glycine
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",

        # stop codons
        "TAA": "*", "TAG": "*", "TGA": "*"
    }

    protein = []

    # Translate in steps of 3 (codons)
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        if len(codon) < 3:
            break  # ignore incomplete codon at the end
        amino_acid = codon_table.get(codon, 'X')  # 'X' for unknown codon if any
        protein.append(amino_acid)

    return "".join(protein)




def main():
    file_path = "dna_sequence.txt"

    num_processes = 5

    # reads/validates DNA sequence
    dna_sequence = read_sequence_from_file(file_path)

    # calculates GC content sequentially
    gc_seq, time_seq = measure_runtime(calculate_sequential, dna_sequence)

    # calculates GC content using ProcessPool Executor (parallelism)
    gc_parallel_processes, time_parallel_processes = measure_runtime(
        calculate_parallel, dna_sequence, num_workers=num_processes, use_multiprocessing=True
    )

    #finds composition of each nucleotide in sequence
    composition = nucleotide_composition(dna_sequence)

    # translates dna sequence to protein (amino acid) sequence
    proteins = dna_protein_translation(dna_sequence)

    # Display Results
    print("\nGC Content Calculation Results")
    print("--------------------------------")
    print(f"GC Content (Sequential): {gc_seq:.2f}%")
    print(f"Sequential Time: {time_seq:.6f} seconds\n")

    print(f"GC Content (Parallel): {gc_parallel_processes:.2f}%")
    print(
        f"Parallel Time: {time_parallel_processes:.6f} seconds\n")
        
    print(composition)

    print(proteins)


if __name__ == "__main__":
    main()
