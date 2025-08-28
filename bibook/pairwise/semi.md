
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

