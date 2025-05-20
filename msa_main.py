import argparse

from msa import *

def display_alignment(msa, seq_names, params, stats):
    """Enhanced graphical display"""
    print("\n" + "=" * 50)
    print(" Multiple Sequence Alignment Results ".center(50))
    print("=" * 50)

    print("\nScoring Parameters:")
    print(f" Match: {params['match']}")
    print(f" Mismatch: {params['mismatch']}")
    print(f" Gap: {params['gap']}")

    print("\nAlignment:")
    print_msa(msa, seq_names)

    print("\nStatistics:")
    print(f" Identity: {stats['identity_percentage']:.2f}%")
    print(f" Matches: {stats['match_count']} | Mismatches: {stats['mismatch_count']}")
    print(f" Gaps: {stats['gap_count']} | Alignment length: {stats['total_columns']}")

def handle_output(msa, params, stats, seq_names, default_path=""):
    """Handle output saving logic"""
    save_choice = input("\nSave results to file? (y/n): ").lower()
    if save_choice == 'y':
        output_path = input(f"Enter output path ({default_path}): ").strip() or default_path
        if output_path:
            save_alignment_to_file(msa, params, stats, output_path, seq_names)
            print(f"Results saved to {output_path}")
        else:
            print("Invalid path. Results not saved.")

def run_center_star_alignment():
    """
    Run the center star alignment algorithm from the console
    """
    print("Multiple Sequence Alignment - Center Star Method")
    print("==============================================")

    # Ask for input method
    print("\nHow would you like to input sequences?")
    print("1. Manual input")
    print("2. Load from FASTA files")

    choice = input("Choose an option (1-2): ")

    sequences = []
    seq_names = []

    if choice == '1':
        num_seqs = int(input("Enter number of sequences: "))

        for i in range(num_seqs):
            seq, name = load_sequence_manual()
            sequences.append(seq)
            seq_names.append(name)

    elif choice == '2':
        # Load from a FASTA file (handles both single and multiple sequences)
        path = input("Enter path to FASTA file: ")
        loaded_sequences, loaded_names = load_sequences_from_fasta(path)

        if loaded_sequences:
            sequences = loaded_sequences
            seq_names = loaded_names
            print(f"Loaded {len(sequences)} sequences from the file.")
        else:
            print("No sequences were loaded.")
            return

    else:
        print("Invalid choice")
        return

    # Check if we have sequences to align
    if not sequences:
        print("No sequences to align.")
        return

    # Scoring scheme
    scoring = get_scoring_scheme()

    # Run the center star alignment
    msa, score_matrix, center_idx = center_star_alignment(
        sequences,
        scoring['match'],
        scoring['mismatch'],
        scoring['gap'],
    )

    stats = compute_alignment_statistics(msa)

    # Display results
    display_alignment(msa, seq_names, scoring, stats)

    # Handle output
    handle_output(msa, scoring, stats, seq_names, "alignment.txt")


def example_star_alignment():
    """
    Run an example of the center star alignment from the lecture example :D
    """
    print("\n" + "=" * 50)
    print(" Example Mode ".center(50))
    print("=" * 50)

    sequences = [
        "ATTGCCATT",  # S1
        "ATGGCCATT",  # S2
        "ATCCATTTTT",  # S3
        "ATCTTCTT",  # S4
        "ACTGACCT"  # S5
    ]
    seq_names = ["S1", "S2", "S3", "S4", "S5"]

    # Get scoring parameters
    print("Default example sequences:")
    for name, seq in zip(seq_names, sequences):
        print(f"{name}: {seq}")

    scoring = get_scoring_scheme()

    # Perform alignment
    msa, score_matrix, center_idx = center_star_alignment(
        sequences,
        match=scoring['match'],
        mismatch=scoring['mismatch'],
        gap=scoring['gap']
    )
    stats = compute_alignment_statistics(msa)

    # Display results
    display_alignment(msa, seq_names, scoring, stats)

    # Handle output
    handle_output(msa, scoring, stats, seq_names, "example_alignment.txt")
def get_scoring_scheme(match=1, mismatch=-1, gap=-2):
    print("\nScoring Scheme Configuration:")
    return {
        'match': int(input(f"Match score (default {match}): ") or match),
        'mismatch': int(input(f"Mismatch penalty (default {mismatch}): ") or mismatch),
        'gap': int(input(f"Gap open penalty (default {gap}): ") or gap),
    }

def main():
    """command-line interface"""
    parser = argparse.ArgumentParser(description='Center Star MSA')
    parser.add_argument('-i', '--input', help='Input FASTA file')
    parser.add_argument('-o', '--output', help='Output file name')
    parser.add_argument('--match', type=int, default=1, help='Score for matches (default: 1)')
    parser.add_argument('--mismatch', type=int, default=-1, help='Score for mismatches (default: -1)')
    parser.add_argument('--gap', type=int, default=-2, help='Penalty for gaps (default: -2)')
    args = parser.parse_args()

    sequences, seq_names = load_sequences_from_fasta(args.input)
    if not sequences:
        print(f"Error loading sequences from {args.input}")
        return

    # Set parameters
    scoring = {
            'match': args.match,
            'mismatch': args.mismatch,
            'gap': args.gap
        }

    # Perform alignment
    msa, score_matrix, center_idx = center_star_alignment(
        sequences,
        match=scoring['match'],
        mismatch=scoring['mismatch'],
        gap=scoring['gap']
    )
    stats = compute_alignment_statistics(msa)

    # Display results
    display_alignment(msa, seq_names, scoring, stats)

    # Handle output
    if args.output:
        save_alignment_to_file(msa, scoring, stats, args.output, seq_names)
        print(f"Results saved to {args.output}")
    else:
        handle_output(msa, scoring, stats, seq_names, "cmd_alignment.txt")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main()
    else:
        print("Multiple Sequence Alignment - Center Star Method")
        print("===============================================")
        print("1. Interactive alignment")
        print("2. Example alignment")
        print("3. Exit")

        choice = input("Choose an option (1-3): ")

        if choice == '1':
            run_center_star_alignment()
        elif choice == '2':
            example_star_alignment()
        else:
            print("Exiting program.")