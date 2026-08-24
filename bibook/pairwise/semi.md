---
kernelspec:
  name: python3
  display_name: Python 3
---


# Semi-global Alignment

## Background

We also have a third type of alignments, semi-global alignments, which are also known as semi-local or glocal alignments.

## Principle

Here we cross breed the Needleman-Wunsch and the Smith-Waterman algorithm, by initiating the dynamic programming matrix as for a local alignment, but using the same recursion as in the global alignment. 

## Implementation

The implementation of the Semi-global alignment algorithm involves initializing an alignment matrix with zeros, and then iteratively updating its cells based on the scores of adjacent cells and the scoring schema (match, mismatch, and gap penalties).

## Applications

It is particularly useful in situations where we strive to align a short sequence to a longer sequence.

## The definitions of the problem, and the solution

Given two sequences $a_1,\ldots,a_N$ and $b_1,\ldots,b_M$, a scoring function $d(x,y)$, find a *semi-global* alignment that gives an optimal (maximal) score.

The solution can be found by studying the dynamic programming matrix, $S$, of size $(N+1,M+1)$, by using the recursions defined in equations {eq}`sg-init` and {eq}`sg-recursion`.

```{math}
:label: sg-init

\begin{align*}
S_{0,0} = 0, & \\
S_{i,0} =& 0, & \textrm{for all}\ i \\
S_{0,j} =& 0, & \textrm{for all}\ j 
\end{align*}
```

```{math}
:label: sg-recursion

S_{i,j}=\max
\begin{cases}
   S_{i-1,j-1} & + d(a_i,b_j)\\
   S_{i-1,j} & + d(a_i,-)\\
   S_{i,j-1} & + d(-,b_j).
\end{cases}
```

The optimal alignment is found by backtracing from the maximal bottommost or rightmost element, $\max(\max_i S_{i,M},\max_j S_{N,j})$, to the first encountered leftmost or topmost (0) element.  

## Example

The typical use of a semi-global alignment is to place a *short* sequence inside a
much *longer* one: a sequencing read against a reference, or a gene against a whole
chromosome. A global alignment would insist on aligning the read to the whole
reference and pay a gap penalty for every position of the reference the read does not
cover, which for a long reference completely swamps the score.

Here we align the short sequence $b=$`GCTTA` against the longer sequence
$a=$`ACGTTGCATTACG`, using the scoring function
$d(x,y)= \begin{cases}3 & \textrm{if } x=y\\-2 & \textrm{if } x= \textrm{- or } y=\textrm{-}\\ -1 & \textrm{otherwise}\end{cases}$.

We begin by initiating the whole first row and the whole first column with zeros,
exactly as for a local alignment, as shown in {numref}`fig-sg-init`. This is what lets
the short sequence start anywhere along the long one, free of charge.

```{code-cell} python
:tags: [hide-input]
:label: fig-sg-init
:caption: Initialization of the dynamic programming matrix. As for a local alignment the first row and the first column are set to zero, so that the short sequence may start anywhere along the long one without cost.

from alignviz import Alignment, Style

# mode="semi" initiates with zeros but keeps the global recursion.
sg = Alignment("ACGTTGCATTACG", "GCTTA", mode="semi")

# Smaller cells, so that the tall matrix fits the page.
compact = Style(cell=44, value_size=16, letter_size=18)

sg.figure(upto=sg.n_init, style=compact)
```

The rest of the matrix is then filled in with the *global* recursion of
Equation {eq}`sg-recursion`. Note that there is no reset to zero here; unlike a local
alignment, a semi-global alignment is free to run through a negative score on its way
to a better one. The filled in matrix is shown in {numref}`fig-sg-fill`.

```{code-cell} python
:tags: [hide-input]
:label: fig-sg-fill
:caption: The filled in matrix. The recursion is the same as for a global alignment, so negative values are kept rather than reset to zero.

sg.figure(style=compact)
```

Finally we trace the alignment backwards, starting from the largest element of the
last row or the last column, and stopping as soon as we reach the first row or the
first column, as shown in {numref}`fig-sg-bt`. That is what lets the short sequence
stop anywhere along the long one, again free of charge.

```{code-cell} python
:tags: [hide-input]
:label: fig-sg-bt
:caption: We trace backwards from the largest element of the last row or column, and stop at the first row or column. The parts of the long sequence outside the aligned region are left unaligned and cost nothing.

sg.figure(path=sg.traceback(), style=compact)
```

The traceback starts in the last column rather than in the bottom-right corner, and
stops in the first column rather than in the top-left corner. It places `GCTTA` over
the stretch `GCATTA` of the long sequence,

```
GCATTA
GC-TTA
```

for a score of $5 \times 3 - 2 = 13$: five matches, and one gap for the `A` of the long
sequence that the short one has nothing to pair with. That internal gap *is* charged.
What is not charged is the leading `ACGTT` and the trailing `CG` of the long sequence,
which simply hang outside the alignment. A global alignment would have had to pay a gap
penalty for each of those seven positions, turning a good hit into a bad score.

## Summary: Comparing Global, Local, and Semi-global Alignment

Sequence alignment algorithms can be broadly categorized into three types: global (Needleman-Wunsch), local (Smith-Waterman), and semi-global alignment. Each serves a distinct purpose depending on the biological question and the nature of the sequences being compared.

- **Global alignment (Needleman-Wunsch)** aligns two sequences from end to end, optimizing the alignment across their entire lengths. It is best suited for sequences of similar length and overall similarity.
- **Local alignment (Smith-Waterman)** finds the highest-scoring matching region (subsequence) between two sequences, making it ideal for identifying conserved domains or motifs within otherwise divergent sequences.
- **Semi-global alignment** is a hybrid approach, allowing for gaps at the ends of one or both sequences without penalty. This is particularly useful when aligning a shorter sequence (such as a sequencing read or a gene) to a longer reference (such as a chromosome or genome).

The table below summarizes the key differences:

| Alignment Type | Matrix Initialization | Recursion Used        | Traceback Start           | Traceback End        | Typical Use Case                        |
|----------------|----------------------|-----------------------|--------------------------|----------------------|------------------------------------------|
| Global         | Gap penalties        | Global (Needleman-Wunsch) | Bottom-right cell         | Top-left cell         | Full-length comparison of similar sequences |
| Local          | Zeros                | Local (Smith-Waterman)    | Max cell anywhere         | First zero           | Motif/domain search, conserved regions      |
| Semi-global    | Zeros                | Global (Needleman-Wunsch) | Max in last row/column    | First zero in row/col | Short-to-long sequence alignment            |

Understanding these differences helps in selecting the appropriate algorithm for a given biological problem.
