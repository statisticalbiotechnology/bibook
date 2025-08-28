# Representations of pairwise alignments

Before delving into the intricacies of the alignment algorithms, it's crucial to understand the matrix representation of alignments. This approach is not just a computational convenience but a fundamental framework that encapsulates the complexity of aligning biological sequences. Here's an overview of why matrices are used in sequence alignment and the significance of the traceback process.

## Understanding Matrix Representation in Sequence Alignments

### The Essence of Matrices in Alignments

Aligning sequences involves comparing each character of one sequence with every character of another to identify regions of similarity or difference. Given the combinatorial explosion of possible alignments with even moderate-length sequences, a brute-force approach to comparison is computationally infeasible. This is where matrices come into play, offering a structured and efficient way to represent all possible alignments between two sequence.

### How Matrices Facilitate Alignment

A matrix for sequence alignment is essentially a grid. Each cell in this grid represents a comparison point between a character from one sequence and another from the second sequence. Each cell of the matrix correspond to the characters in the two sequences being aligned. And a path from the first left cell to the bottom right cell of the matrix represent a global alignment.


```{figure} ./img/movement.png
---
name: fig-movement
scale: 50 %
alt: One step of an alignment
---
A step in the path representing an alignment.

If your alignment is in position $a_i,b_j$ of to sequences $\mathbf{a}=(a_1,\ldots,a_N)$ and $\mathbf{b}=(b_1,\ldots,b_N)$ then the next step in the alignment is either a <span style="color:green">match/mismatch</span> between $a_{i+1},b_{j+1}$ (diagonal movement), a <span style="color:magenta">delete</span> $a_{i+1},-$, (a vertical movement), or an <span style="color:orange">insert</span> $-,b_{j+1}$, (a horizontal movement).

```

### Alignments as paths between the cells

The path we follow tells us how the sequences should be aligned. A move diagonally signifies a match or mismatch, indicating that the characters from both sequences should be aligned with each other. Moving up or to the left indicates a gap in one of the sequences.


### Examles of alignments

An example of some matrix representations of an alignment is found in {numref}`fig-matrix-representation`.

```{figure} ./img/matrix_representation.png
---
name: fig-matrix-representation
scale: 50 %
alt: Matrix representation of three alignments
---
Matrix representation of three full length alignments of `ACGTACT` and `ACTACGT`, as well as one partial alignment of the same sequences.

The path of Alignment 1 is given in red, Alignment 2 in green, and Alignment 3 in blue. We also included a partial (local) alignment, Alignment 4 in orange.

| Alignment 1 (red) | Alignment 2 (green) | Alignment 3 (blue) | Alignment 4 (orange)
| :---------- | :---------- | :---------- | :---------- |
| `---ACGTACT`  | `ACGTAC-T`    | `ACGTACT----` | `ACG` |
| `ACTACGT---`  | `AC-TACGT`    | `----ACTACGT` | `ACG` |

```

The matrix representation of alignments and the traceback process are not merely computational techniques; they are integral to understanding the relationships between biological sequences. By breaking down the alignment process into a series of quantifiable steps, matrices allow us to systematically explore the vast space of alignment possibilities, leading us to the optimal alignment that best reflects the evolutionary or functional relationship between the sequences involved.

### Excercises

````{exercise} Matrix Representations of Alignments
-------
:label: ex-matrixrep
-------

* Which alignment does the following path represent?

```{figure} ./img/matrix_representation_exer.png
---
name: fig-matrix-representation-exer
scale: 50 %
alt: Exercise on Matrix representation of alignments
---

```


```{dropdown} **Reveal Answer**

GATTA-GA   
G-CTACTA

```
````

## Why Matrix Representations Are Useful

Matrices are more than just a convenient way to visualize alignments — they form the mathematical foundation that allows us to apply efficient algorithms. The key property that makes them powerful is that they provide a **structured way to exploit the principle of optimal substructure**.

### Optimal Substructure in Alignments

When aligning two sequences, the score of the *best alignment* up to a given position $(i,j)$ depends only on:
- the characters $a_i$ and $b_j$, and  
- the best scores of alignments leading into that position (from $(i-1,j-1)$, $(i-1,j)$, or $(i,j-1)$).  

This means that once we have computed the optimal score up to $(i,j)$, we no longer need to remember *how* we reached that point. The matrix cell itself summarizes all relevant information about the past.

In other words, **any optimal scoring pathway of the subsequent part of the alignment is independent of the specific choices made earlier**. This independence allows us to build solutions incrementally, storing partial results in the matrix and avoiding redundant computations.
