import sys
import matplotlib.pyplot as plt
import numpy as np

# import seaborn as sns # Not used in the relevant parts, can be kept if used elsewhere
# import matplotlib.patches as mpatches # Not explicitly used for this task, can be kept

"""
Quantitative matching of sequence pairs
Prepare a properly functioning program for quantitative comparison of pairs of protein or DNA
coding sequences using the Needleman-Wunsch global matching algorithm
The program should control the correctness of the input data, generate correct results, as well as be
concise, legible (with comments) and consistent with the style rules of the programming language
used.
IMPORTANT NOTE:
sequence
matching point
SHOW how the match between2 sequences
"""


def load_sequence_manual():
    """
    Prompt the user to enter two sequences
    :return: sequence1, sequence2
    """
    seq1 = input("Please enter the first sequence: ").strip().upper()
    seq2 = input("Please enter the second sequence: ").strip().upper()
    # Basic validation: ensure sequences are not empty and contain valid characters (optional, add as needed)
    if not seq1 or not seq2:
        print("Error: Both sequences must be provided.")
        sys.exit(1)
    # Example: check for allowed characters (e.g., for DNA)
    # allowed_chars = "ATGC"
    # if not all(c in allowed_chars for c in seq1) or not all(c in allowed_chars for c in seq2):
    #     print(f"Error: Sequences can only contain characters from {allowed_chars}.")
    #     sys.exit(1)
    return seq1, seq2


def load_sequence_fasta(path1, path2):
    """
    Load one sequence from each of two FASTA files.

    Args:
        path1 (str): path to first FASTA file (one record)
        path2 (str): path to second FASTA file (one record)
    Returns: seq1, seq2 (str)
    """

    def _read_single_fasta(path):
        seq_lines = []
        try:
            with open(path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('>'):
                        continue
                    seq_lines.append(line.upper())  # Convert to upper for consistency
        except FileNotFoundError:
            print(f"Error: File not found at {path}")
            sys.exit(1)
        except Exception as e:
            print(f"Error reading FASTA file {path}: {e}")
            sys.exit(1)

        if not seq_lines:
            print(f"No sequence data found in {path}")
            sys.exit(1)
        return ''.join(seq_lines)

    seq1 = _read_single_fasta(path1)
    seq2 = _read_single_fasta(path2)
    return seq1, seq2


def dotplot(seq1, seq2, window_size=1, threshold=1):
    """
     Dot plot table from sequence 1 and sequence 2
    :param seq1: first sequence
    :param seq2: second sequence
    :param window_size: size of window. Default is 1
    :param threshold: threshold for dot plot. Default is 1
    :return: the dot plot table - 1 when they have the same character, otherwise 0
    """
    len1 = len(seq1)
    len2 = len(seq2)
    matrix = np.zeros((len1, len2))

    if window_size == 1:
        for i in range(len1):
            for j in range(len2):
                if seq1[i] == seq2[j]:
                    matrix[i][j] = 1
    else:
        for i in range(len1 - window_size + 1):
            for j in range(len2 - window_size + 1):
                matches = 0
                for k in range(window_size):
                    # Check bounds: i+k must be < len1 and j+k < len2
                    if seq1[i + k] == seq2[j + k]:
                        matches += 1
                if matches >= threshold:
                    # Mark the start of the window match.
                    # Or, more commonly, mark all cells in the window diagonal.
                    # For this implementation, we'll mark the start [i,j]
                    # Or, to fill the diagonal segment:
                    # for l_k in range(window_size):
                    #    matrix[i+l_k][j+l_k] = 1
                    # Sticking to original logic of marking matrix[i][j]
                    matrix[i][j] = 1
    return matrix


def dotplot2Ascii(dp, seq1, seq2, heading, filename):
    """
    Print the dot plot table from sequence 1 and sequence 2 to the text file
    :param dp: dot plot table
    :param seq1: sequence 1 (to take the first sequence name)
    :param seq2: sequence 2 (to take the second sequence name)
    :param heading: Heading name of the result file
    :param filename: filename to save the result
    :return: a (text) file of dot plot table with the heading and 2 sequences
    """
    try:
        with open(filename, 'w') as f:
            f.write(heading + "\n\n")
            f.write("      " + " ".join(list(seq2)) + "\n")
            f.write("    +" + "-" * (2 * len(seq2) - 1) + "\n")
            for i in range(len(seq1)):
                f.write(f"{seq1[i]}   | ")
                for j in range(len(seq2)):
                    if dp[i][j] == 1:
                        f.write("* ")
                    else:
                        f.write(". ")
                f.write("\n")
        print(f"Dot plot table saved to {filename}")
    except Exception as e:
        print(f"Error saving dot plot: {e}")


def dotplotGraphic(dp, labelA, labelB, heading, filename="dotplot.png"):
    """
    Generate a graphical dot plot.
    :param dp: dot plot matrix (numpy array)
    :param labelA: sequence for Y-axis
    :param labelB: sequence for X-axis
    :param heading: title for the plot
    :param filename: name of the file to save the plot
    """
    fig, ax = plt.subplots(figsize=(max(6, len(labelB) * 0.3), max(6, len(labelA) * 0.3)))

    # Use imshow for a cleaner plot if dp is dense, or scatter if sparse
    # For typical dot plots, scatter is fine.
    rows, cols = np.where(dp == 1)
    ax.scatter(cols, rows, marker='o', color='black', s=max(100 / max(len(labelA), len(labelB)), 1))

    ax.set_xlabel(f"Sequence B: {labelB if len(labelB) < 30 else labelB[:27] + '...'}")
    ax.set_ylabel(f"Sequence A: {labelA if len(labelA) < 30 else labelA[:27] + '...'}")
    ax.set_title(heading)

    # Set ticks to be in the center of the cells, corresponding to characters
    ax.set_xticks(np.arange(len(labelB)))
    ax.set_yticks(np.arange(len(labelA)))

    # Set tick labels to sequence characters
    ax.set_xticklabels(list(labelB))
    ax.set_yticklabels(list(labelA))

    # Invert Y axis to have origin at top-left, common for matrices
    ax.invert_yaxis()

    # Add grid lines
    ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    try:
        plt.savefig(filename)
        print(f"Dot plot graphic saved to {filename}")
    except Exception as e:
        print(f"Error saving dot plot graphic: {e}")
    plt.show()


def scoring_matrix_multiple_paths(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Calculate the scoring matrix and a direction matrix that stores ALL optimal paths
    for sequence alignment using Needleman-Wunsch.
    Directions: 0 for diagonal, 1 for up, 2 for left.
    """
    len1, len2 = len(seq1), len(seq2)
    score_m = np.zeros((len1 + 1, len2 + 1), dtype=int)

    # direction_m stores lists of optimal directions from each cell
    direction_m = [[[] for _ in range(len2 + 1)] for _ in range(len1 + 1)]

    # Initialize first row and column for scores and directions
    for i in range(1, len1 + 1):
        score_m[i, 0] = i * gap
        direction_m[i][0] = [1]  # From Up
    for j in range(1, len2 + 1):
        score_m[0, j] = j * gap
        direction_m[0][j] = [2]  # From Left
    # direction_m[0][0] remains [] (it's the starting point for traceback)

    # Fill the matrices
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            char_s1 = seq1[i - 1]
            char_s2 = seq2[j - 1]

            match_score_val = match if char_s1 == char_s2 else mismatch

            diag_val = score_m[i - 1, j - 1] + match_score_val
            up_val = score_m[i - 1, j] + gap
            left_val = score_m[i, j - 1] + gap

            current_max_score = max(diag_val, up_val, left_val)
            score_m[i, j] = current_max_score

            # Store all directions that lead to the max score
            if diag_val == current_max_score:
                direction_m[i][j].append(0)  # Diagonal
            if up_val == current_max_score:
                direction_m[i][j].append(1)  # Up
            if left_val == current_max_score:
                direction_m[i][j].append(2)  # Left

            # Optional: sort directions for deterministic behavior if multiple paths start identically
            # direction_m[i][j].sort()

    return score_m, direction_m


def find_all_optimal_alignments(seq1, seq2, score_matrix, direction_matrix):
    """
    Traceback through the direction matrix to find all optimal alignments.
    Returns a list of tuples: (aligned_seq1, aligned_seq2, path_coordinates_list)
    """
    all_found_paths = []

    # Recursive helper function for traceback
    def backtrack_recursive(r_i, r_j, current_aligned_s1_chars, current_aligned_s2_chars, current_path_coords_list):
        # Add current cell (score matrix index) to this path's coordinates
        current_path_coords_list.append((r_i, r_j))

        # Base case: reached the top-left corner (0,0)
        if r_i == 0 and r_j == 0:
            all_found_paths.append(
                ("".join(reversed(current_aligned_s1_chars)),
                 "".join(reversed(current_aligned_s2_chars)),
                 list(reversed(current_path_coords_list)))  # Store a reversed copy of coords
            )
            current_path_coords_list.pop()  # Backtrack: remove current cell for other explorations
            return

        # Get possible moves from the direction matrix for cell (r_i, r_j)
        possible_moves = direction_matrix[r_i][r_j]

        for move_type in possible_moves:
            if move_type == 0:  # Diagonal move
                if r_i > 0 and r_j > 0:  # Ensure we can move diagonally
                    backtrack_recursive(r_i - 1, r_j - 1,
                                        current_aligned_s1_chars + [seq1[r_i - 1]],
                                        current_aligned_s2_chars + [seq2[r_j - 1]],
                                        current_path_coords_list)
            elif move_type == 1:  # Up move
                if r_i > 0:  # Ensure we can move up
                    backtrack_recursive(r_i - 1, r_j,
                                        current_aligned_s1_chars + [seq1[r_i - 1]],
                                        current_aligned_s2_chars + ['-'],
                                        current_path_coords_list)
            elif move_type == 2:  # Left move
                if r_j > 0:  # Ensure we can move left
                    backtrack_recursive(r_i, r_j - 1,
                                        current_aligned_s1_chars + ['-'],
                                        current_aligned_s2_chars + [seq2[r_j - 1]],
                                        current_path_coords_list)

        current_path_coords_list.pop()  # Backtrack: remove current cell after exploring all its children

    # Start traceback from the bottom-right corner of the score matrix
    initial_i, initial_j = len(seq1), len(seq2)
    backtrack_recursive(initial_i, initial_j, [], [], [])

    if not all_found_paths and (len(seq1) > 0 or len(seq2) > 0):  # Check if sequences are not both empty
        # This case should ideally not be reached if scoring/direction matrix is correct,
        # as NW always finds a path.
        print("Warning: No optimal alignments were found during traceback. This might indicate an issue.")
        # You could add a default path based on the old single-path traceback if needed as a fallback,
        # but it's better to ensure the multi-path logic is sound.

    return all_found_paths


def visualize_alignment_multiple_paths(seq1, seq2, score_m, all_paths_coords, max_len=100,
                                       filename="alignment_matrix.png"):
    """
    Visualize the alignment scoring matrix and highlight all optimal traceback paths.
    :param seq1: First sequence (string)
    :param seq2: Second sequence (string)
    :param score_m: The scoring matrix (NumPy array)
    :param all_paths_coords: A list of path coordinate lists (each list contains (row,col) tuples for one path)
    :param max_len: Maximum number of rows/columns to display for very large sequences
    :param filename: Name of the file to save the plot
    """
    rows, cols = score_m.shape

    # Adjust for visualization grid (adding one row/col for sequence labels)
    # Determine visual dimensions, respecting max_len
    viz_display_rows = min(rows + 1, max_len + 1 if max_len < rows else rows + 1)
    viz_display_cols = min(cols + 1, max_len + 1 if max_len < cols else cols + 1)

    fig_width = max(8, 0.3 * viz_display_cols)
    fig_height = max(8, 0.3 * viz_display_rows)
    plt.figure(figsize=(fig_width, fig_height))
    ax = plt.gca()

    # Collect all unique cells that are part of any optimal path for highlighting
    # Path coordinates are 0-indexed score matrix indices (r, c)
    # Visualization grid is 1-indexed, so a cell (r,c) in score_m is at (r+1, c+1) in viz grid
    highlight_cells_viz_indices = set()
    if all_paths_coords:
        for single_path_coords in all_paths_coords:
            for r_idx, c_idx in single_path_coords:
                highlight_cells_viz_indices.add((r_idx + 1, c_idx + 1))

    default_font_size = 8
    # Adjust font size dynamically based on matrix size to prevent overlap
    font_size = max(2, min(default_font_size, default_font_size * 30 / max(1, rows, cols)))

    # Draw the grid and fill cells for the visualized portion
    for viz_r in range(viz_display_rows):  # viz_r is 0-indexed visual row
        for viz_c in range(viz_display_cols):  # viz_c is 0-indexed visual col

            # Cell border
            rect = plt.Rectangle((viz_c, viz_r), 1, 1, fill=False, edgecolor='gray', linewidth=0.5)
            ax.add_patch(rect)

            text_x, text_y = viz_c + 0.5, viz_r + 0.5  # Center of the cell

            # Top-left corner cell (0,0 in viz grid)
            if viz_r == 0 and viz_c == 0:
                ax.text(text_x, text_y, " ø ", ha='center', va='center', fontweight='bold', fontsize=font_size)

            # Sequence 2 characters on the first visual row (viz_r = 0)
            # score_m column j corresponds to viz_c = j+1
            # seq2 char at index j-1 corresponds to score_m column j, so viz_c = j+1
            elif viz_r == 0 and viz_c > 0:  # Skips (0,0)
                if viz_c - 1 < len(seq2) and viz_c - 1 < cols:  # viz_c-1 is effectively seq2 index
                    ax.text(text_x, text_y, seq2[viz_c - 1], ha='center', va='center', fontsize=font_size,
                            fontweight='bold')

            # Sequence 1 characters on the first visual col (viz_c = 0)
            # score_m row i corresponds to viz_r = i+1
            # seq1 char at index i-1 corresponds to score_m row i, so viz_r = i+1
            elif viz_c == 0 and viz_r > 0:  # Skips (0,0)
                if viz_r - 1 < len(seq1) and viz_r - 1 < rows:  # viz_r-1 is effectively seq1 index
                    ax.text(text_x, text_y, seq1[viz_r - 1], ha='center', va='center', fontsize=font_size,
                            fontweight='bold')

            # Score values from the score_m
            # score_m[r_score, c_score] is at viz_r = r_score+1, viz_c = c_score+1
            # So, r_score = viz_r-1, c_score = viz_c-1
            elif viz_r > 0 and viz_c > 0:  # These are the score cells
                r_score, c_score = viz_r - 1, viz_c - 1
                if r_score < rows and c_score < cols:  # Check if within score_m bounds
                    # Highlight if part of any optimal path
                    if (viz_r, viz_c) in highlight_cells_viz_indices:  # Check using viz_grid indices
                        fill_color = 'lightblue' if (r_score, c_score) != (
                        rows - 1, cols - 1) else 'mediumseagreen'  # Different color for final score
                        rect_fill = plt.Rectangle((viz_c, viz_r), 1, 1, fill=True, facecolor=fill_color, alpha=0.6)
                        ax.add_patch(rect_fill)

                    ax.text(text_x, text_y, str(score_m[r_score, c_score]),
                            ha='center', va='center', fontsize=font_size,
                            color='black' if not (viz_r,
                                                  viz_c) in highlight_cells_viz_indices else 'black')  # Ensure text is readable

    # Display overall alignment score
    alignment_score = score_m[rows - 1, cols - 1]  # Score is at the bottom-right of score_m
    # Position score text outside the main grid, adjust as needed
    # ax.text(viz_display_cols / 2, -0.5, f"Optimal Score: {alignment_score}", ha='center', va='bottom', fontsize=font_size + 2, fontweight='bold')
    plt.figtext(0.5, 0.01, f"Optimal Alignment Score: {alignment_score}", ha="center", fontsize=12,
                bbox={"facecolor": "lightgray", "alpha": 0.5, "pad": 5})

    ax.set_xlim(-0.5, viz_display_cols - 0.5)
    ax.set_ylim(viz_display_rows - 0.5, -0.5)  # Inverted y-axis for matrix style (0 at top)

    ax.set_xticks([])  # No numeric ticks on axes
    ax.set_yticks([])

    for spine in ax.spines.values():  # Remove frame spines
        spine.set_visible(False)

    plt.title(f"Needleman-Wunsch Alignment Matrix ({seq1} vs {seq2})", fontsize=14, y=1.05)
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])  # Adjust layout to make space for figtext

    try:
        plt.savefig(filename)
        print(f"Alignment matrix visualization saved to {filename}")
    except Exception as e:
        print(f"Error saving alignment matrix visualization: {e}")
    plt.show()


def run_needleman_wunsch_paths_visualization(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Run the Needleman-Wunsch algorithm to find ALL optimal paths and visualize the results.
    """
    print(f"\nAligning Sequence 1: {seq1}")
    print(f"Aligning Sequence 2: {seq2}")
    print(f"Parameters: Match={match}, Mismatch={mismatch}, Gap={gap}\n")

    score_matrix, direction_m = scoring_matrix_multiple_paths(seq1, seq2, match, mismatch, gap)

    all_alignments_data = find_all_optimal_alignments(seq1, seq2, score_matrix, direction_m)

    if not all_alignments_data:
        print("No optimal alignments found. This may indicate an issue or very unusual sequences/parameters.")
        # Optionally, you can still try to visualize the score matrix if it's generated
        # visualize_alignment_multiple_paths(seq1, seq2, score_matrix, [], filename="alignment_matrix_no_path.png")
        return None, None, None  # Or handle as an error

    print(f"Found {len(all_alignments_data)} optimal alignment(s).")
    print(f"Optimal Score: {score_matrix[len(seq1), len(seq2)]}\n")

    all_path_coords_lists = []
    for idx, (al1, al2, p_coords) in enumerate(all_alignments_data):
        print(f"Alignment #{idx + 1}:")
        # Create a match string: | for match, . for mismatch, space for gap
        match_str = "".join(
            "|" if al1[k] == al2[k] and al1[k] != '-' else
            ("." if al1[k] != '-' and al2[k] != '-' else " ")
            for k in range(len(al1))
        )
        print(f"  {al1}")
        print(f"  {match_str}")
        print(f"  {al2}\n")
        all_path_coords_lists.append(p_coords)

    visualize_alignment_multiple_paths(seq1, seq2, score_matrix, all_path_coords_lists,
                                       filename=f"NW_{seq1}_vs_{seq2}_matrix.png")

    return all_alignments_data, score_matrix


    # --- The following are older/alternative implementations from the original code ---
    # --- They are not directly used by the new run_needleman_wunsch_M_paths_visualization ---
def scoring_path(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Calculate the scoring matrix for sequence alignment (original single path logic)
    """
    len1, len2 = len(seq1), len(seq2)
    path = np.zeros((len1 + 1, len2 + 1), dtype=int)
    direction = np.zeros_like(path, dtype=int)  # 0 diagonal, 1 Up, 2 left

    for i in range(1, len1 + 1):
        path[i, 0] = i * gap
        direction[i, 0] = 1
    for j in range(1, len2 + 1):
        path[0, j] = j * gap
        direction[0, j] = 2

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            diag_score = path[i - 1, j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            up_score = path[i - 1, j] + gap
            left_score = path[i, j - 1] + gap
            scores = [(diag_score, 0), (up_score, 1), (left_score, 2)]
            best_score, best_dir = max(scores, key=lambda x: (x[0], -x[1]))  # Tie-breaking
            path[i, j] = best_score
            direction[i, j] = best_dir
    return path, direction


def traceback_alignment(seq1, seq2, direction):
    """
    Trace back through the direction matrix to find the optimal alignment.

    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        direction (numpy.ndarray): Direction matrix

    Returns:
        tuple: (aligned_seq1, aligned_seq2, alignment_score)
    """
    i, j = len(seq1), len(seq2)
    aligned_seq1, aligned_seq2 = "", ""
    score = 0

    while i > 0 or j > 0:
        if i > 0 and j > 0 and direction[i, j] == 0:  # Diagonal
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            if seq1[i - 1] == seq2[j - 1]:
                score += 1  # Match
            else:
                score -= 1  # Mismatch
            i -= 1
            j -= 1
        elif i > 0 and direction[i, j] == 1:  # Up
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            score -= 2  # Gap
            i -= 1
        elif j > 0:  # Left
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            score -= 2  # Gap
            j -= 1

    return aligned_seq1, aligned_seq2, score

def scoring_path2(seq1, seq2, match=1, mismatch=-1, gap=-2):  # Original alternative
    len1, len2 = len(seq1), len(seq2)
    score = np.zeros((len1 + 1, len2 + 1), dtype=int)
    direction_m = np.zeros_like(score, dtype=int)
    for i in range(1, len1 + 1): score[i, 0] = i * gap; direction_m[i, 0] = 1
    for j in range(1, len2 + 1): score[0, j] = j * gap; direction_m[0, j] = 2
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            diag = score[i - 1, j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            up = score[i - 1, j] + gap
            left = score[i, j - 1] + gap
            best, dir_idx = max((diag, 0), (up, 1), (left, 2), key=lambda x: (x[0], -x[1]))
            score[i, j], direction_m[i, j] = best, dir_idx
    return score, direction_m


def traceback(seq1, seq2, direction_m):  # Original alternative
    i, j = len(seq1), len(seq2)
    a1, a2 = [], []
    while i > 0 or j > 0:
        if i == 0:
            a1.append('-'); a2.append(seq2[j - 1]); j -= 1
        elif j == 0:
            a1.append(seq1[i - 1]); a2.append('-'); i -= 1
        else:
            d = direction_m[i, j]
            if d == 0:
                a1.append(seq1[i - 1]); a2.append(seq2[j - 1]); i -= 1; j -= 1
            elif d == 1:
                a1.append(seq1[i - 1]); a2.append('-'); i -= 1
            else:
                a1.append('-'); a2.append(seq2[j - 1]); j -= 1
    return ''.join(reversed(a1)), ''.join(reversed(a2))


def align_pair(seq1, seq2, match=1, mismatch=-1, gap=-2):  # Uses original alternatives
    score_matrix, direction_matrix = scoring_path2(seq1, seq2, match, mismatch, gap)
    return traceback(seq1, seq2, direction_matrix)


def format_and_print_alignment(al1, al2):  # Original
    L = max(len(al1), len(al2))
    al1 = al1.ljust(L, '-')
    al2 = al2.ljust(L, '-')
    print(" ".join(list(al1)))  # Added list() for clarity if al1/al2 are strings
    print(" ".join(list(al2)))


if __name__ == "__main__":
    # Example
    seq_C = "AG"
    seq_D = "AAG"  # A-G vs AAG, or AG- vs AAG
    # run_needleman_wunsch_M_paths_visualization(seq_C, seq_D, match=1, mismatch=-1, gap=-1)

    seq_E = "ATTGCCATT"
    seq_F = "ATCCAATTTT"
    print(f"\n--- Running Needleman-Wunsch for {seq_E} and {seq_F} ---")
    run_needleman_wunsch_paths_visualization(seq_E, seq_F, match=1, mismatch=-1, gap=-2)
