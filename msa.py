import sys
from itertools import combinations
from multiprocessing import Pool

import numpy as np
from typing import List, Tuple, Dict
import textwrap

from nw2 import scoring_path, traceback_alignment

def load_sequence_manual():
    """
    Load a sequence from manual input.

    Returns:
        tuple: (sequence, name)
    """
    name = input("Enter sequence name: ")
    seq = input(f"Enter sequence for {name}: ").strip().upper()
    return seq, name


def load_sequences_from_fasta(file_path):
    """
    Load multiple sequences from a FASTA file.

    Args:
        file_path (str): Path to the FASTA file

    Returns:
        tuple: (sequences, seq_names)
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()

        sequences = []
        seq_names = []
        current_seq = ""
        current_name = ""

        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                # Save the previous sequence if there was one
                if current_seq:
                    sequences.append(current_seq.upper())
                    seq_names.append(current_name)

                # Start a new sequence
                current_name = line[1:]  # Remove the '>' character
                current_seq = ""
            elif line:  # Only add non-empty lines
                current_seq += line

        # Add the last sequence
        if current_seq:
            sequences.append(current_seq.upper())
            seq_names.append(current_name)

        return sequences, seq_names
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return [], []
    except Exception as e:
        print(f"Error loading FASTA file: {str(e)}")
        return [], []

def worker(params):
    """
    Worker function for parallel computation of pairwise alignments

    Args:
        params: Tuple of ((i, seq_i), (j, seq_j), (match, mismatch, gap))

    Returns:
        Tuple: ((i, j), (aligned_i, aligned_j, score))
    """
    (i, seq_i), (j, seq_j), scoring_params = params
    match, mismatch, gap = scoring_params

    score_table, direction = scoring_path(seq_i, seq_j, match, mismatch, gap)
    aligned_i, aligned_j, score = traceback_alignment(seq_i, seq_j, direction, match, mismatch, gap)

    return ((i, j), (aligned_i, aligned_j, score))

def compute_all_pairwise_alignments(sequences, match=1, mismatch=-1, gap=-2):
    """
    Compute all pairwise alignments between sequences and their scores
    Args:
        sequences (list): List of sequences to align
        match (int): Score for matching
        mismatch (int): Score for mismatching
        gap (int): Score for gap penalty

    Returns:
        tuple: (score_matrix, alignment_pairs)
            - score_matrix: Matrix of alignment scores between sequences
            - alignment_pairs: Dictionary of pairwise alignments
    """
    n = len(sequences)
    score_matrix = np.zeros((n, n))
    alignment_pairs = {}

    # Create tasks for parallelization
    tasks = []
    for i, j in combinations(range(n), 2):
        tasks.append(((i, sequences[i]), (j, sequences[j]), (match, mismatch, gap)))

    # Parallel process
    with Pool() as pool:
        results = pool.map(worker, tasks)

    # Result process
    for (i, j), (aligned_i, aligned_j, score) in results:
        alignment_pairs[(i, j)] = (aligned_i, aligned_j)
        score_matrix[i, j] = score
        score_matrix[j, i] = score

    return score_matrix, alignment_pairs


def find_center_sequence(score_matrix):
    """
    Find the center sequence (most similar to all others)

    Args:
        score_matrix (np.array): Matrix of alignment scores
    Returns:
        int: Index of center sequence
    """
    row_sums = np.sum(score_matrix, axis=1)
    center_idx = np.argmax(row_sums)
    return int(center_idx)


def align_similar(s1, s2):
    """
    Find positions where gaps need to be inserted to make two sequences identical

    Args:
        s1 (str): First sequence
        s2 (str): Second sequence

    Returns:
        tuple: (changes_for_s1, changes_for_s2) - indices where gaps need insertion
    """
    change1, change2 = [], []
    i = 0

    # Copy sequences to avoid modifying originals
    s1_copy, s2_copy = s1, s2

    while s1_copy != s2_copy:
        if i >= len(s1_copy):
            # Need to extend s1 with remainder of s2
            s1_copy += s2_copy[i:]
            change1.extend(range(i, i + len(s2_copy[i:])))
            break

        if i >= len(s2_copy):
            # Need to extend s2 with remainder of s1
            s2_copy += s1_copy[i:]
            change2.extend(range(i, i + len(s1_copy[i:])))
            break

        if s1_copy[i] != s2_copy[i]:
            if s1_copy[i] == '-':
                # Insert gap in s2
                s2_copy = s2_copy[:i] + '-' + s2_copy[i:]
                change2.append(i)
            else:
                # Insert gap in s1
                s1_copy = s1_copy[:i] + '-' + s1_copy[i:]
                change1.append(i)

        i += 1

    return sorted(change1), sorted(change2)


def adjust(string_list, indices):
    """
    Insert gaps at specified positions in all strings in a list

    Args:
        string_list (list): List of strings to modify
        indices (list): Positions where gaps should be inserted
    """
    for idx in range(len(string_list)):
        s = string_list[idx]
        for pos in indices:
            if pos <= len(s):
                s = s[:pos] + '-' + s[pos:]
        string_list[idx] = s


def merge_alignments(center_idx, alignments, sequences, match=1, mismatch=-1, gap=-2):
    """
    Merge pairwise alignments into a multiple sequence alignment

    Args:
        center_idx (int): Index of center sequence
        alignments (dict): Dictionary of pairwise alignments
        sequences (list): Original sequences
        match, mismatch, gap: Scoring parameters

    Returns:
        list: Multiple sequence alignment
    """
    n = len(sequences)
    other_indices = [i for i in range(n) if i != center_idx]

    # Start with the center sequence
    msa_rows = [sequences[center_idx]]

    # Process other sequences in order
    for idx in other_indices:
        # Get pairwise alignment between center and current sequence
        if (center_idx, idx) in alignments:
            center_aligned, seq_aligned = alignments[(center_idx, idx)]
        else:
            center_aligned, seq_aligned = alignments[(idx, center_idx)]

        # If this is the first sequence added to the MSA
        if len(msa_rows) == 1:
            msa_rows = [center_aligned, seq_aligned]
            continue

        # Find positions where gaps need to be inserted
        ch_index1, ch_index2 = align_similar(msa_rows[0], center_aligned)

        # Create copy of current MSA and the new sequence to add
        new_msa = msa_rows.copy()
        new_seq = [seq_aligned]

        # Insert gaps in MSA and new sequence
        adjust(new_msa, ch_index1)
        adjust(new_seq, ch_index2)

        # Add new sequence to MSA
        msa_rows = new_msa + new_seq

    return msa_rows

def center_star_alignment(sequences, match=1, mismatch=-1, gap=-2):
    """
    Perform multiple sequence alignment using the center star method.

    Args:
        sequences (list): List of sequences to align
        match (int): Score for matching characters
        mismatch (int): Score for mismatching characters
        gap (int): Score for gap penalty

    Returns:
        tuple: (msa, score_matrix, center_idx)
            - msa: List of aligned sequences in the MSA
            - score_matrix: Matrix of alignment scores
            - center_idx: Index of the center sequence
    """
    print("Computing pairwise alignments...")
    score_matrix, alignments = compute_all_pairwise_alignments(sequences, match, mismatch, gap)
    print(f"Score matrix:\n{score_matrix}")

    print("Finding center sequence...")
    center_idx = find_center_sequence(score_matrix)
    print(f"Center sequence is sequence #{center_idx + 1} ({sequences[center_idx]})")

    print("Merging alignments...")
    msa = merge_alignments(center_idx, alignments, sequences, match, mismatch, gap)

    return msa, score_matrix, center_idx


def compute_alignment_statistics(msa):
    """
    Compute statistics for the multiple sequence alignment.

    Args:
        msa (list): List of aligned sequences

    Returns:
        dict: Dictionary with statistics
    """
    n = len(msa)
    seq_length = len(msa[0])

    stats = {
        "total_columns": seq_length,
        "identity_percentage": 0,
        "match_count": 0,
        "mismatch_count": 0,
        "gap_count": 0,
        "conserved_columns": 0,
        "column_stats": []
    }

    # Per-column analysis
    for col in range(seq_length):
        column = [seq[col] for seq in msa]
        gaps = column.count('-')
        non_gaps = [c for c in column if c != '-']

        if len(non_gaps) == 0:
            continue

        most_common = max(set(non_gaps), key=non_gaps.count)
        matches = sum(1 for c in non_gaps if c == most_common)
        mismatches = len(non_gaps) - matches

        stats['gap_count'] += gaps
        stats['match_count'] += matches
        stats['mismatch_count'] += mismatches

        if matches == len(non_gaps) and gaps == 0:
            stats['conserved_columns'] += 1

        stats['column_stats'].append({
            'position': col + 1,
            'gaps': gaps,
            'matches': matches,
            'mismatches': mismatches
        })

    total_possible = n * (n - 1) / 2 * seq_length
    stats['identity_percentage'] = (stats['match_count'] / (stats['match_count'] + stats['mismatch_count'])) * 100 \
        if (stats['match_count'] + stats['mismatch_count']) > 0 else 0
    return stats

def print_msa(msa, seq_names=None):
    """
    Print the multiple sequence alignment in a readable format

    Args:
        msa (list): List of aligned sequences
        seq_names (list): Optional list of sequence names
    """
    if not seq_names:
        seq_names = [f"Seq {i + 1}" for i in range(len(msa))]

    # Find max length of sequence names for proper formatting
    max_name_len = max(len(name) for name in seq_names)

    print()

    # Print the MSA
    for i, seq in enumerate(msa):
        print(f"{seq_names[i]:{max_name_len}} ", end="")
        for j, char in enumerate(seq):
            print(char, end="")
        print()


def compute_identity_matrix(msa):
    """
    Compute the percentage identity matrix from an MSA

    Args:
        msa (list): List of aligned sequences

    Returns:
        np.array: Matrix of percent identities between sequences
    """
    n = len(msa)
    identity_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            # Count matches and total non-gap positions
            matches = 0
            total = 0

            for k in range(len(msa[0])):
                if msa[i][k] != '-' and msa[j][k] != '-':  # Neither is a gap
                    total += 1
                    if msa[i][k] == msa[j][k]:  # Match
                        matches += 1

            identity = matches / total if total > 0 else 0
            identity_matrix[i, j] = identity
            identity_matrix[j, i] = identity

    return identity_matrix


def save_alignment_to_file(msa: List[str], params: Dict, stats: Dict, filename: str, seq_names: List[str] = None):
    """
    Save alignment with parameters and statistics to file
    """
    with open(filename, 'w') as f:
        f.write("PROGRAM PARAMETERS:\n")
        f.write(f"Match score: {params['match']}\n")
        f.write(f"Mismatch score: {params['mismatch']}\n")
        f.write(f"Gap penalty: {params['gap']}\n")

        f.write("MULTIPLE SEQUENCE ALIGNMENT:\n")
        max_name_len = max(len(name) for name in seq_names) if seq_names else 6
        for i, seq in enumerate(msa):
            name = seq_names[i] if seq_names else f"Seq {i + 1}"
            wrapped_seq = textwrap.wrap(seq, 80)
            for j, line in enumerate(wrapped_seq):
                prefix = f"{name:{max_name_len}} " if j == 0 else " " * (max_name_len + 1)
                f.write(f"{prefix}{line}\n")

        f.write("\nSTATISTICS:\n")
        f.write(f"Identity percentage: {stats['identity_percentage']:.2f}%\n")
        f.write(f"Total matches: {stats['match_count']}\n")
        f.write(f"Total mismatches: {stats['mismatch_count']}\n")
        f.write(f"Total gaps: {stats['gap_count']}\n")
        f.write(f"Alignment length: {stats['total_columns']}\n")