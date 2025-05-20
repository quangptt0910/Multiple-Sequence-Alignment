import sys
import numpy as np
from typing import List, Tuple, Dict
import textwrap

from nw2 import scoring_path2, traceback_alignment

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

    for i in range(n):
        for j in range(i + 1, n):
            # Compute pairwise alignment
            score_table, direction = scoring_path2(sequences[i], sequences[j], match=match, mismatch=mismatch, gap=gap)

            # Get the alignment
            aligned_i, aligned_j, score = traceback_alignment(sequences[i], sequences[j], direction)

            score_matrix[i, j] = score
            score_matrix[j, i] = score

            alignment_pairs[(i, j)] = (aligned_i, aligned_j)
            alignment_pairs[(j, i)] = (aligned_j, aligned_i)

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
    return np.argmax(row_sums)


def merge_alignments(center_seq_idx, alignments, sequences):
    """
    Correctly merges alignments with proper gap propagation
    """
    n = len(sequences)
    center_seq = sequences[center_seq_idx]
    aligned_seqs = [list(center_seq)]  # Start with center sequence

    # Initialize all aligned sequences with gaps from the center's pairwise alignments
    for i in range(n):
        if i == center_seq_idx:
            continue

        # Get pairwise alignment between center and sequence i
        center_aligned, other_aligned = alignments[(center_seq_idx, i)]

        # Track gap positions to insert in all sequences
        gap_positions = []
        current_pos = 0
        for c in center_aligned:
            if c == '-':
                gap_positions.append(current_pos)
            else:
                current_pos += 1

        # Insert gaps into all existing sequences
        for pos in reversed(gap_positions):
            for seq in aligned_seqs:
                seq.insert(pos, '-')

        # Add the newly aligned sequence
        aligned_seqs.append(list(other_aligned))

    # Ensure equal length
    max_len = max(len(seq) for seq in aligned_seqs)
    for seq in aligned_seqs:
        while len(seq) < max_len:
            seq.append('-')

    return [''.join(seq) for seq in aligned_seqs]

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

    print("Finding center sequence...")
    center_idx = find_center_sequence(score_matrix)
    print(f"Center sequence is sequence #{center_idx + 1}")

    print("Merging alignments...")
    msa = merge_alignments(center_idx, alignments, sequences)

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