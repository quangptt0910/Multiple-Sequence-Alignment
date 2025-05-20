# Multiple Sequence Alignment (MSA) - Center Star Method

## ![Bioinformatics](https://img.shields.io/badge/Bioinformatics-MSA-blue)

This repository provides an implementation of the Center Star Multiple Sequence Alignment (MSA) algorithm in Python. The Center Star method is a heuristic approach that constructs a multiple alignment by first computing all pairwise alignments, selecting the "center" sequence (the one with the highest sum of pairwise scores), and then merging the pairwise alignments into a final multiple alignment.

## Repository Contents

* `msa.py` and `nw2.py`: Core modules for pairwise alignment, scoring, traceback, and the Center Star MSA.
* `msa_main.py`: Command-line interface and interactive/exemplary modes.
* Utility functions for FASTA I/O, statistics computation, and result display.

## Algorithm Scheme

1. **Load sequences**
   Sequences can be provided manually or loaded from a FASTA file.

2. **Compute pairwise alignments**
   Use a Needleman–Wunsch variant to score and trace back every pair of sequences.

3. **Build score matrix**
   Store each pairwise alignment score in an $n \times n$ matrix.

4. **Select center sequence**
   Sum each row of the score matrix; choose the sequence with maximal total score as the center.

5. **Merge alignments**
   Starting from the center sequence, iteratively merge each pairwise alignment by propagating gaps across all previously aligned sequences to maintain column consistency.

6. **Compute MSA statistics**
   Calculate identity percentage, conserved columns, gaps, matches, and mismatches per column and overall.

## Requirements

* Python 3.7+
* `numpy` 
* `matplotlib`(in case for visualize the nw2.py)
* No external alignment libraries are required; all algorithms are implemented in the repository.

Install dependencies via pip:

```bash
pip install numpy matplotlib
```

## Time and Space Complexity Analysis

### Time Complexity

* **Pairwise alignments**: Needleman–Wunsch on two sequences of length $L$ takes $O(L^2)$ time. For $n$ sequences, there are $\binom{n}{2} = O(n^2)$ pairs, giving $O(n^2 \cdot L^2)$.
* **Center selection**: Summation over each row of an $n\times n$ matrix is $O(n^2)$.
* **Merging alignments**: Each merge may require propagating up to $O(L)$ gaps over up to $O(n)$ sequences, for $O(nL)$ per merge, repeated $n-1$ times: $O(n^2 L)$.

Overall dominant cost: **$O(n^2 L^2)$**.

### Space Complexity

* **Score matrix**: $O(n^2)$ entries.
* **Dynamic programming table**: for each pairwise alignment, uses $O(L^2)$ space.
* **MSA result**: stores $n$ sequences of length up to $O(L+g)$ where $g$ is total gaps, worst-case $O(nL)$.

Overall auxiliary space: **$O(n^2 + L^2)$** per alignment, plus final alignment storage $O(nL + ng)$.

## Example Comparisons of Sequence Groups

We illustrate performance on 2 example fasta file, each have 4 evolutionarily related/unrelated:
*With scoring scheme match = 1, mismatch = -1 and gap = -2*
## **Related** 
### related_msa.fa
```text
>Homolog1
ATGGCGTACTGAT
>Homolog2
ATGGGGTATCAA
>Homolog3
ATGGCGACTATG
>Homolog4
ATGGCCGACTGA
```
```bash
# Result
Alignment:

Homolog1 ATGGCGTACTGAT-
Homolog2 ATGGGGTA-TCAA-
Homolog3 ATGGCG-ACT-ATG
Homolog4 ATGGCCGACTGA--

Statistics:
 Identity: 89.80%
 Matches: 44 | Mismatches: 5
 Gaps: 7 | Alignment length: 14
```
## **Unrelated** 
```text
>ProteinA
MVLSPADKTN
>DNAB
GATCCCTA
>EnzymeX
LLYPTAER
>Random
ACGTACGT
```
```bash
# Result
Alignment:

ProteinA --ACGTACG-T-
DNAB     -MVLSPADK-TN
EnzymeX  --GATCCCTA--
Random   LLYPTA-ER---

Statistics:
 Identity: 47.06%
 Matches: 16 | Mismatches: 18
 Gaps: 14 | Alignment length: 12
```
## Usage

```bash
# Interactive mode:
python msa_main.py

# Command-line FASTA alignment:
python msa_main.py -i input.fasta -o output.txt --match 1 --mismatch -1 --gap -2
```

Follow on-screen prompts to configure scoring and save results.

## Output

* Aligned sequences printed to console with identity statistics.
* Optional file output includes parameters, wrapped alignments, and statistics.

---
