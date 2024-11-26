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


def calculate_parallel(sequence, num_workers=5, use_multiprocessing=False):
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
                ' ', '').upper()  # formats just in case its messed up
            return sequence
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        sys.exit(1)
    except ValueError as ve:
        print(f"Error: {ve}")
        sys.exit(1)


# measures the runtime
def measure_runtime(function, *args, **kwargs):
    start_time = time.time()
    result = function(*args, **kwargs)
    elapsed = time.time() - start_time
    return result, elapsed


def main():
    file_path = "dna_sequence.txt"

    num_processes = 5

    # reads/validates DNA sequence
    dna_sequence = read_sequence_from_file(file_path)

    # calculates GC content sequentially
    gc_seq, time_seq = measure_runtime(calculate_sequential, dna_sequence)

    # Parallel GC Content Calculation using ProcessPoolExecutor
    # calculates GC content using ProcessPool Executor (parallelism)
    gc_parallel_processes, time_parallel_processes = measure_runtime(
        calculate_parallel, dna_sequence, num_workers=num_processes, use_multiprocessing=True
    )

    # Display Results
    print("\nGC Content Calculation Results")
    print("--------------------------------")
    print(f"GC Content (Sequential): {gc_seq:.2f}%")
    print(f"Sequential Time: {time_seq:.6f} seconds\n")

    print(f"GC Content (Parallel): {gc_parallel_processes:.2f}%")
    print(
        f"Parallel Time: {time_parallel_processes:.6f} seconds\n")


if __name__ == "__main__":
    main()
    