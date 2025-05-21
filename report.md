# Multiple Sequence Alignment (MSA) - Center Star Method

#### *https://github.com/quangptt0910/Multiple-Sequence-Alignment*

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
*https://www.uniprot.org/uniref/UniRef90_A0A009L5F2*
![img.png](related_result.png)
### related_msa.fa
```text
>tr|A0A009L5F2|A0A009L5F2_9GAMM Metallopeptidase M24 family protein OS=Acinetobacter sp. 1475718 OX=1310652 GN=J551_4367 PE=4 SV=1
MIAGYYTALERTQFFDHCDDASLRIWEANVKVHEAGLKLIRPGARCSDIAKELNEIFYEE
GLLQYRTFGYGHSFGVLSHYYGREAGLELREDIDTVIEPGMVISMEPMIMIPQGMAGAGG
YREHDILVVTETGNENITKFPYGPEHNIIK
>UPI00366FA0D9 status=active
YNVVYTDWQRDNYFKAVTQEIPNKGRVGIEFDHLPLERLNKLKAVLPNVQFVDISAACMK
LRMIKSAEEIAHITQGAAVCDIGGFAVTAAIREGVPEHEVALASTQSMVREIAKRFPNSE
LMDTWTWFQSGINTDGAHNPVTTRQVRKGDILSLNCFSMIAGYYTALERTQFFDHCDDAS
LRIWEANVKVHEAGLKLIRPGARCSDIAKELNEIFYEEGLLQYRTFGYGHSFGVLSHYYG
REAGLELREDIDTVIEPGMVISMEPMIMIPQGMAGAGGYREHDILVVTETGNENITKFPY
GPEHNIIK
>UPI0004F540D4 status=active
SLNCFSMIAGYYTALERTQFFDHCDDASLRIWEANVKVHEAGLKLIRPGARCSDIAKELN
EIFYEEGLLQYRTFGYGHSFGVLSHYYGREAGLELREDIDTVIEPGMVISMEPMIMIPQG
MAGAGGYREHDILVVTETGNENITKFPYGPEHNIIK
>UPI0012502E68 status=active
NINYYADFIYCSFGRFYGLVVSPEKVVTISANIDAGQPWRRTIGDYNVVYTDWQRDNYFK
AVAQEIPNKGRVGIEFDHLPLERLNKLKAVLPNVQFVDISSACMKLRMIKSAEEIAHITQ
GAAVCDIGGFAVTAAIREGVPEHEVALASTQSMVREIAKRFPNSELMDTWTWFQSAINTD
GAHNPVTTRQVQKGDILSLNCFSMIAGYYTALERTQFFDHCDDASLRIWEANVKVHEAGL
KLIRPGARCSDIAKELNEIFYEEGLLQYRTFGYGHSFGVLSHYYGREAGLELREDIDTVI
EPGMVISMEPMIMIPQGMAGAGGYREHDL
>UPI0004F57FE2 status=active
QLMEELEIDHVLFSSIHNINYYADFIYCSFGRFYGLVVSPEKVVTISANIDAGQPWRRTV
GDYNVVYTDWQRDNYCKAVAQEIPNRGRVGIEFDHLPLERLNKLKTVLPNVQFVDISAAC
MKLRMIKSAEEIAHITQGAAVCDIGGFAVTAAIREGVPEHEVALASTQSMVREIAKRFPN
SELMDTWTWFQSGINTDGAHNPVTTRQVQKGDILSLNCFSMIAGYYTALERTQFFDHCDD
ASLRIWEANVKVHEAGLKLIRPGARCSDIAKELNEIFYAEGLLQ
```
```bash
# Result
...
Alignment:

tr|A0A009L5F2|A0A009L5F2_9GAMM MetallopeptidaseM24... ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------MIAGYYTALERTQFFDHCDDASLRIWEANVKVHEAGLKLIRPGARCSDIAKELNEIFYEEGLLQYRTFGYGHSFGVLSHYYGREAGLELREDIDTVIEPGMVISMEPMIMIPQGMAGAGGYREHDILVVTETGNENITKFPYGPEHNIIK
UPI00366FA0D9 status=active                           --------------------------------------------------------------YNVVYTDWQRDNYFKAVTQEIPNKGRVGIEFDHLPLERLNKLKAVLPNVQFVDISAACMKLRMIKSAEEIAHITQGAAVCDIGGFAVTAAIREGVPEHEVALASTQSMVREIAKRFPNSELMDTWTWFQSGINTDGAHNPVTTRQVRKGDILSLNCFSMIAGYYTALERTQFFDHCDDASLRIWEANVKVHEAGLKLIRPGARCSDIAKELNEIFYEEGLLQYRTFGYGHSFGVLSHYYGREAGLELREDIDTVIEPGMVISMEPMIMIPQGMAGAGGYREHDILVVTETGNENITKFPYGPEHNIIK
UPI0004F540D4 status=active                           ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------SLNCFSMIAGYYTALERTQFFDHCDDASLRIWEANVKVHEAGLKLIRPGARCSDIAKELNEIFYEEGLLQYRTFGYGHSFGVLSHYYGREAGLELREDIDTVIEPGMVISMEPMIMIPQGMAGAGGYREHDILVVTETGNENITKFPYGPEHNIIK
UPI0012502E68 status=active                           NINYYADFIYCSFGRFYGLVVSPEKVVTISANIDAGQPWRRTIGD-----------------YNVVYTDWQRDNYFKAVAQEIPNKGRVGIEFDHLPLERLNKLKAVLPNVQFVDISSACMKLRMIKSAEEIAHITQGAAVCDIGGFAVTAAIREGVPEHEVALASTQSMVREIAKRFPNSELMDTWTWFQSAINTDGAHNPVTTRQVQKGDILSLNCFSMIAGYYTALERTQFFDHCDDASLRIWEANVKVHEAGLKLIRPGARCSDIAKELNEIFYEEGLLQYRTFGYGHSFGVLSHYYGREAGLELREDIDTVIEPGMVISMEPMIMIPQGMAGAGGYREHD-L-----------------------
UPI0004F57FE2 status=active                           QLMEELEIDHVLFSSIHNINYYADFIYCSFGRFYGLVVSPEKVVTISANIDAGQPWRRTVGDYNVVYTDWQRDNYCKAVAQEIPNRGRVGIEFDHLPLERLNKLKTVLPNVQFVDISAACMKLRMIKSAEEIAHITQGAAVCDIGGFAVTAAIREGVPEHEVALASTQSMVREIAKRFPNSELMDTWTWFQSGINTDGAHNPVTTRQVQKGDILSLNCFSMIAGYYTALERTQFFDHCDDASLRIWEANVKVHEAGLKLIRPGARCSDIAKELNEI-----------------F------Y-AE-GL-L-----------------------Q-------------------------------------

Statistics:
 Identity: 95.76%
 Matches: 1175 | Mismatches: 52
 Gaps: 623 | Alignment length: 370
 ...
```
![related.png](related.png)

## **Unrelated** 
*https://www.uniprot.org/uniref/UniRef50_A0A009MDJ4*
![related_result.png](related_result.png)
```text
>tr|A0A009MDJ4|A0A009MDJ4_9GAMM Major Facilitator Superfamily protein OS=Acinetobacter sp. 1475718 OX=1310652 GN=J551_3466 PE=4 SV=1
MSTEYSQPSAIERETMRKVAWRLLPFLTLCYVIAIIDRGNIGMASLQMNEDLGLTATAFG
FASSLFFFAYFLFEVPSNLAMQKFGAKVWIARIMVTWGGNIIRNCIR
>tr|A0A829S2Y2|A0A829S2Y2_9GAMM Major Facilitator Superfamily protein OS=Acinetobacter sp. 225588 OX=1310788 GN=J687_3779 PE=4 SV=1
MSTEYSQPSAIERETMRKVAWRLLPFLTLCYVIAIIDRGNIGMASLQMNEDLGLTATAFG
FASSLFFFAYFLFEVPSNLAMQKFGAKVWIARIMVTWGGNIIRNCIR
>tr|A0A951B0D1|A0A951B0D1_9PROT MFS transporter (Fragment) OS=Acetobacteraceae bacterium OX=1909293 GN=JOZ42_14295 PE=4 SV=1
MDNAAVERSAMRRVGWRLVPFLIFAYFVAFLDRVNVGFAAIQMNKDVGLTATIYGWGAGI
FFLGYFVLEVPSNLALERFGARIWIA
>UPI001C5A93E3 status=active
ALPIFGMASLQMNEDLGLTATAFGFASSLFFFAYFLFEVPSNLAMQKFGAKVWIARIMVT
WGVISSATAFVENSTSLYIMRFLLGAAEAGFFPGVILYLTYWL
```
```bash
# Result
...
Alignment:

tr|A0A009HVK9|A0A009HVK9_ACIB9 Poly-beta-1,6-N-acetyl-D-glucosamine synthase OS=Acinetobacter baumannii (strain 1295743) OX=1310613 GN=J512_0835 PE=3 SV=1 M--L-V------S---KY--EY--L---V----CIDGDALLHPHAVLWLMQPFLNFPRIGAVTGNPRILNRSSILGKLQVGEFSSIIGLIKRAQRTYGRIFTVSGVIAAFRKTALVRVGFWSDDKITEDIDISWKLQMDHWDIQYIPQALCYIYMPETFKGLWKQRLRWAQGGVEVLLEYIPKMFKLRLRRMWPVMLEALVSIIW--SYVMIMIFILFFVGLFVDLPQQFQINSLMPQWYGVILGGTCLVQFLVSLWIDHRYDRGRLFRNYFWVIWYPLFFWLLTLFTSVVAVPKTIF-NTKKRARWVSPDRGFRGDHS-
tr|A0A378B747|A0A378B747_KLEPO Biofilm PGA synthesis N-glycosyltransferase PgaC OS=Klebsiella pneumoniae subsp. ozaenae OX=574 GN=pgaC PE=3 SV=1           M-------------I-----HY--------------------P--------H------VGAVTGNPRIRTRSTLIGRIQVGEFSSIIGLIKRTQRIYGRVFTVSGVIAAFRRQALADVGYWSPDMITEDIDISWKLQLRHWDIFFEPRALCWILMPETLKGLWKQRLRWAQGGAEVFLVNLRRLFHWQHHRMWPLFLEYACSTLWAFAYAITML-LFILSRLMPLPDNLAVHSLFPPEFTGMLLGLMCLLQFVASLYIERRYEKK-VAKSLFWVIWFPMVYWMIGLLTTLVAFPKVMVKRQRSRARWVSPNRG-K--GSL
tr|A0A6N0ZAD0|A0A6N0ZAD0_9PROT Glycosyltransferase OS=Accumulibacter sp OX=2053492 GN=HT579_08560 PE=3 SV=1                                                MTPLPVLLGDFWSAILAYVFFYPVLMSVVWMIGGIDGDARLDPCALNWLVRHF-DEPDVAAVTGNPRVRNRSTLLGRIQVGEFSSIVGLIKRAQRSVGLIFTVAGVIAAFRRSALHEVGYWSPEKLTEDVDISWKLQLAGWEVRFEPRALCWILMPETLRGLWKQRLRWSMGGTQVLLDYTPQLLSPKNAAMWPLCVEYGMSITW--AYLFTLLAAYRVIDLLVYRIDLQSAPILLLGWAGLLIATTCLLQIMPSLVLDRPYDRG-LLKNYFWMIWYPIIYWAITAATSVVAVPKILFRDSEKRARWTSPDRGIKPDASR
UPI0035ADA118 status=active                                                                                                                                MTPLPVLLGDFWSAILAYVFFYPVLMSVVWMIGGIDGDARLDPCALNWLVRHF-DEPDVAAVTGNPRVRNRSPLLGRIQVGEFSSIVGLIKRAQRSVGLIFTVAWVIAAFRHSALHEVGYWSPEKLTEDVDISWKLQLAGWEVRFEPRALCWILMPETLRGLWKQRLRWSMGGTQVLLDYTPQLLSPKNAAMWPLCVEYGMSITS--AYLFALLATYRVIDLLVYRIDLQSAPILLLWWTGLLIATTCLLQIMLSLVLDRPYDRE-PLKNYFWMIWYPIIYWAITAATSVVAVPKILFRDSEKRARWTSPDRGIKPDASR

Statistics:
 Identity: 78.28%
 Matches: 930 | Mismatches: 258
 Gaps: 92 | Alignment length: 320
 ...
```
![unrelated.png](unrelated.png)

## Usage

```bash
# Interactive mode:
python msa_main.py

# Command-line FASTA alignment:
python msa_main.py -i example.fa -o output.txt --match 1 --mismatch -1 --gap -2
```

Follow on-screen prompts to configure scoring and save results.

## Output

* Aligned sequences printed to console with identity statistics.
* Optional file output includes parameters, wrapped alignments, and statistics.

---
