---
kernelspec:
  name: python3
  display_name: Python 3
---

# Needleman-Wunsch Algorithm: Global Alignment

## Background

Developed by Saul B. Needleman and Christian D. Wunsch in 1970 {cite}`needleman1970general`, the Needleman-Wunsch algorithm was one of the first computational approaches to sequence alignment. Its introduction marked a significant advancement in bioinformatics, enabling the systematic and automated comparison of biological sequences.

## Principle

The algorithm is based on dynamic programming, a method that breaks down complex problems into simpler, smaller subproblems, solving each just once and storing their solutions. In the context of sequence alignment, it constructs an optimal global alignment by comparing every character of one sequence with every character of another, considering the costs of matches, mismatches, and gaps. The logic is that we build up an alignment one position at a time, using previously computed optimal alignments. We begin in the top left corner of the matrix representation of the alignment, where we still have not aligned any characters, and hence a score of zero. We then proceed to fill in the matrix, one cell at a time, until we reach the bottom right corner, which contains the score of the optimal alignment. 

## The definitions of the problem, and the solution

Given two sequences $a_1,\ldots,a_N$ and $b_1,\ldots,b_M$, a scoring function $d(x,y)$, find a *global* alignment that gives an optimal (maximal) score.

The solution can be found by studying the dynamic programming matrix, $S$, of size $(N+1,M+1)$, by using the recursions defined in Equations {eq}`nw-init` and {eq}`nw-recursion`.

```{math}
:label: nw-init

S_{0,0} = 0, \\
S_{i,0} = d(x,-) \cdot i \ \textrm{for all}\ i, \\
S_{0,j} =  d(-,y) \cdot j\ \textrm{for all}\ j 
```

```{math}
:label: nw-recursion

S_{i,j}=\max
\begin{cases}
   S_{i-1,j-1} & + d(a_i,b_j)\\
   S_{i-1,j} & + d(a_i,-)\\
   S_{i,j-1} & + d(-,b_j).
\end{cases}
```

We will walk through these steps more carefully below.

## Process

1. **Initialization:** It starts by creating a scoring matrix where one sequence is aligned along the top and the other along the side. The first row and column are filled with gap penalties, increasing progressively to set up the basis for the algorithm according to Equation {eq}`nw-init`.

1. **Matrix Filling:** Each cell in the matrix is then filled based on the scores of adjacent cells (top, left, and diagonal), plus the score for matching or mismatching the corresponding characters, or introducing a gap. The choice of score at each cell reflects the highest score achievable from the possible alignments up to that point. Here we use the recursion given by Equation {eq}`nw-recursion`.

1. **Traceback:** Once the matrix is filled, the optimal alignment is determined by tracing back from the bottom-right corner to the top-left, following the path that resulted in the highest score. This path represents the optimal global alignment of the two sequences.

2. **Alignment Output:** The traceback path is used to construct the aligned sequences, introducing gaps as necessary, to maximize the alignment score based on the predefined scoring system. The final score of the alignment is given by the element, $S_{N,M}$.

## Applications and Importance

The Needleman-Wunsch algorithm is fundamental when the goal is to align entire sequences, providing a comprehensive view of their similarity. It's particularly valuable in evolutionary biology for comparing homologous sequences across different species, helping to infer phylogenetic relationships and evolutionary events. Moreover, it lays the foundation for understanding the principles of dynamic programming in bioinformatics, influencing the development of other alignment algorithms and tools.

Despite its computational intensity, especially for long sequences, the Needleman-Wunsch algorithm remains a crucial method for global sequence alignment, embodying the essence of comparing biological sequences in a mathematically rigorous and systematic way.

## Examples

### Example 1: Short sequences

Here is an example alignment of the sequences `GAC` and `ACG` using Needleman-Wunsch, when we have a scoring function,
$d(x,y)= \begin{cases}1 & \textrm{if} x=y\\ -1 & \textrm{otherwise } \end{cases}$.

We start by filling in the borders of the matrix using the Equation {eq}`nw-init`, as shown in {numref}`fig-nw-short-1`.

```{code-cell} python
:tags: [hide-input]
:label: fig-nw-short-1
:caption: Initialization of the gap penalty matrix for the Needleman-Wunsch algorithm.

from alignviz import Alignment

# d(x,y) = 1 for a match and -1 otherwise, with a gap penalty of -1.
short = Alignment("GAC", "ACG", match=1, mismatch=-1, gap=-1)

# frames() gives the initialised borders, then one figure per cell of the
# recursion in the order they are filled in, and finally the traceback.
frames = short.frames()

short.figure(upto=short.n_init)
```

We then recursively fill in the other elements of the matrix in a row wise manner using Equation {eq}`nw-recursion`, as shown in {ref}`fig-nw-recursion`.

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-step-1

frames[1]
```

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-step-2

frames[2]
```

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-step-3

frames[3]
```

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-step-4

frames[4]
```

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-step-5

frames[5]
```

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-step-6

frames[6]
```

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-step-7

frames[7]
```

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-step-8

frames[8]
```

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-step-9

frames[9]
```

(fig-nw-recursion)=
**Filling in the NW matrix**

::::{grid} 3

:::{grid-item}
![](#cell-nw-step-1)
:::

:::{grid-item}
![](#cell-nw-step-2)
:::

:::{grid-item}
![](#cell-nw-step-3)
:::

:::{grid-item}
![](#cell-nw-step-4)
:::

:::{grid-item}
![](#cell-nw-step-5)
:::

:::{grid-item}
![](#cell-nw-step-6)
:::

:::{grid-item}
![](#cell-nw-step-7)
:::

:::{grid-item}
![](#cell-nw-step-8)
:::

:::{grid-item}
![](#cell-nw-step-9)
:::
::::

We fill in the elements recursively, in a row-wise manner. Each cell's value is evaluated using Equation {eq}`nw-recursion`. The cell currently being computed is shaded, together with the (up to) three cells the recursion reads from, and the recursion itself is spelled out under each image. We store trackers of which step we used to reach a certain cell, indicated by red arrows. Note that for some cells there are multiple optimal steps, i.e. paths that have the same score.

Given the filled in matrix, we can now track the optimal path from the bottom right element of the matrix, following the arrows back to the top-left element, as shown in {numref}`fig-nw-bt`.

```{code-cell} python
:tags: [hide-input]
:label: fig-nw-bt
:caption: We trace the alignment backwards from the bottom-right corner to the top-left corner of the matrix, and mark the found optimal path with blue arrows.

short.figure(path=short.traceback())
```

### Example 2: Longer sequences

Here we align the sequences $a=$TGCATTA $b=$GCATTAC when $\displaystyle d(x,y)= \begin{cases}3 & \textrm{if} x=y\\-2 & \textrm{if} x= \textrm{- or } y=\textrm{-}\\  -1 & \textrm{otherwise } \end{cases}$. 

The resulting matrix is found in {numref}`fig-nw-long`.

```{code-cell} python
:tags: [hide-input]
:label: fig-nw-long
:caption: We follow the alignment backwards from the bottom-right corner to the top-left corner of the matrix, and mark the found optimal path with blue arrows.

# match=3, mismatch=-1, gap=-2 is the default scoring of alignviz
long = Alignment("TGCATTA", "GCATTAC")
long.figure(path=long.traceback())
```

### Exercise

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-exe1

exe1 = Alignment("GATTA", "GCTAC")
exe1.figure(path=exe1.traceback())
```

```{code-cell} python
:tags: [remove-cell]
:label: cell-nw-exe2

exe2 = Alignment("GCAGCTA", "GCTA")
exe2.figure(path=exe2.traceback())
```

::::{exercise} Needleman-Wunsch Alignment 1
:label: ex-nwexe1

Calculate the Needleman Wunsch Alignment of the following two sequences:

```
GATTA  
GCTAC
```

Use the following scoring scheme:
- Match: +3
- Mismatch: -1
- Gap penalty: -2

:::{dropdown} **Reveal Answer**
![](#cell-nw-exe1)
:::
::::

::::{exercise} Needleman-Wunsch Alignment 2
:label: ex-nwexe2

Calculate the Needleman-Wunsch Alignment of th following two sequences:

```
GCAGCTA   
GCTA
```

Use the following scoring scheme:
- Match: +3
- Mismatch: -1
- Gap penalty: -2

:::{dropdown} **Reveal Answer**
![](#cell-nw-exe2)
:::
::::

## Big-O Notation

Big-O notation is used in computational science for describing how the running time or memory usage of an algorithm scales with a given factor. E.g. if we expect the running time to scale as $g(x)$ we write that the algorithm has complexity $\mathcal{O}(g(x))$. A more formal definition can be found at [wikipedia](https://en.wikipedia.org/wiki/Big_O_notation).

In the case of Needleman-Wunsch we see that the number of calculations needed are proportional to the size of the dynamic programming matrix, which equals the product of the lengths of the sequences, M x N. This results in a time complexity of $ \mathcal{O}(MN) $, indicating that the time to complete the task scales proportionally with the product of the lengths of the two sequences.

In the same way memory usage also scales with $ \mathcal{O}(MN)$, as the scoring matrix used to store intermediate results requires memory proportional to its size.

Big-O notation serves as a quick and effective tool for comparing different algorithms. For example, it allows us to see at a glance how the Needleman-Wunsch algorithm compares to other sequence alignment algorithms in terms of efficiency.

A useful comparison is the complexity of our initial proposition, to enumerate and calculate the scores for all possible alignments of two sequences. This can be done by calculating the number of alignments with $k$ matches/mis-matches between the two sequences which is ${M \choose k}{N \choose k}$. If we assume that $N>M$ and sum this for all possible values of $k$, we get $\sum_{k=0}^M{M \choose k}{N \choose k}=\sum_{k=0}^M{M \choose M-k}{N \choose k}={N+M \choose M}=\frac{(M+N)!}{M!*N!}$ number of different alignments. This can be [shown](https://math.stackexchange.com/a/4134185) to follow $\mathcal{O}((\frac{e(N+M)}{M})^M)$ {cite}`lange2002mathematical, eddy2004dynamic`.

## Try it yourself

The figures in this chapter are all drawn by the small module `alignviz.py`, which
also drives an interactive app at
[alignviz.serve.scilifelab.se](https://alignviz.serve.scilifelab.se). There you can
type in two sequences of your own, change the scores for matches, mismatches and
gaps, and watch the dynamic programming matrix and its traceback redraw as you type.

It is well worth a few minutes of experimentation. Some things to try:

- Make the gap penalty harsher and harsher, and see at what point the alignment
  stops using gaps altogether.
- Look for sequence pairs where a cell has more than one incoming arrow. Those are
  the cases where several different alignments share the same optimal score.
- Align a short sequence against a much longer one, and see what a global alignment
  does with the overhanging ends. The Smith-Waterman and semi-global algorithms of
  the coming chapters are the answer to exactly that problem, and the app can draw
  them as well.

```{bibliography}
:filter: docname in docnames
```
