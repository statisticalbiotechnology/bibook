# Code: Global pairwise alignments

We will now make use of our definitions of a Needleman-Wunsch alignment to see how the algorithm transforms to actual code. The sections below will walk you through how this is done. 

````{toggle} Service functions for formatting and printing alignments
To facilitate the the reasoning in the subsequent cells, we first we define a couple of service functions that we will need later, for formating and printing alignments. It is not important that you understand what these functions do, for now.

```python
import numpy as np

# Print 2 sequences on top of each other
def print_alignment(seqA,seqB):
    print(seqA)
    print(seqB)

# Print the dynamic programming score matrix (S) together with the sequences.
# This helps visualize how the alignment scores are computed.
def print_dynamic(seqA,seqB,S,trace):
    seqA,seqB = "-" + seqA, "-" + seqB  # Add leading '-' for alignment with matrix indices
    m,n = len(seqA),len(seqB)
    print('{:^5}'.format(" "), end=""),  # Print header row
    for j in range(n):
        print('{:^5}'.format(seqB[j]), end="")
    print()
    for i in range(m):
        print ('{:^5}'.format(seqA[i]), end="")  # Print row label
        for j in range(n):
            print ('{:5.1f}'.format(S[i,j]), end="")  # Print score
        print()
    print()

# Format an alignment by inserting gaps ('-') in the sequences according to the trace matrix.
# Returns the aligned sequences as strings.
def format_alignment(seqA,seqB,S,trace):
    print("Best score: " + str(S[-1,-1]))  # Print the optimal alignment score
    outA,outB = "",""  # Output aligned sequences
    i,j = len(seqA),len(seqB)
    while i>0 or j>0:
        di,dj = trace[i,j]  # Get the direction from the trace matrix
        i += int(di)
        j += int(dj)
        if di == 0:
            outA = "-" + outA  # Insert gap in seqA
        else:
            outA = seqA[i] + outA  # Add character from seqA
        if dj == 0:
            outB = "-" + outB  # Insert gap in seqB
        else:
            outB = seqB[j] + outB  # Add character from seqB
    return outA,outB
```
````

## Scoring system for DNA sequences

We setup the scoring system we need for the alignment of DNA sequences. Here we use a score system where gaps score -2 and miss matches are scored -1 and matches get a score of 3. You should absolutely try to change the scoring function at some point to see how the values affect the alignments. 

```python
# Returns the penalty for introducing a gap in the alignment.
def gap_penalty():
    return -2.0

# Returns the score for aligning two letters (could be nucleotides or amino acids).
# Matches get a positive score, mismatches a negative score, and gaps use the gap penalty.
def match_score(letterA,letterB):
    if letterA == '-' or letterB == '-':
        return gap_penalty()  # Gap penalty if either is a gap
    elif letterA == letterB:
        return 3.0  # Match score
    else:
        return -1.0  # Mismatch penalty
```

## Global alignments by Needleman-Wunsch

Here follows an implementation of the [Needleman-Wunsch](http://www.sciencedirect.com/science/article/pii/0022283670900574) pairwise alignment method.  We want to align two sequences $a=a_1,\ldots,a_{M}$ and $b=b_1,\ldots,b_{N}$. 

As said in last chapter, the dynamic programming matrix $S$ is initiated as:
$S_{i0}=g \cdot i, \forall i,$
$S_{0j}=g \cdot j, \forall j$

This translates easily into two for loops filling in the upper row and the leftmost column.

### The trace matrix: keeping track of the optimal path

To reconstruct the optimal alignment after filling the score matrix $S$, we need to know which move (match/mismatch, insertion, or deletion) led to each cell's score. For this, we use a trace matrix. Each element in the trace matrix contains a pair of values $(di, dj)$ indicating the direction of the optimal step that led to the current cell:

- $(-1, -1)$: Diagonal move (match or mismatch, aligning $a_i$ with $b_j$)
- $(-1, 0)$: Upwards move (gap in sequence B, aligning $a_i$ with a gap)
- $(0, -1)$: Leftwards move (gap in sequence A, aligning $b_j$ with a gap)

By following these pointers from the bottom-right cell back to the top-left, we can reconstruct the optimal alignment by inserting gaps where needed. This process is called traceback.

Here we initiate pointers in the form of a trace matrix. Each element in the trace matrix contains the difference in index between the current cell and the optimal step that lead to the current cell. 

```python
# Initialize the dynamic programming matrices for global alignment.
# S: score matrix, trace: matrix to keep track of the optimal path.
def initiate_global_dp(m,n):
    S = np.zeros((m,n))       # Score matrix of size m x n, initialized to 0
    trace = np.zeros((m,n,2)) # Trace matrix to store direction of optimal move
    S[0,0] = 0.
    trace[0,0,:] = (0.,0.)
    # Initialize first column: gap penalties for seqA
    for i in range(1,m):
        S[i,0] = i * gap_penalty()
        trace[i,0,:] =(-1,0)  # Came from above (gap in seqB)
    # Initialize first row: gap penalties for seqB
    for j in range(1,n):
        S[0,j] = j * gap_penalty()
        trace[0,j,:] =(0,-1)  # Came from left (gap in seqA)
    return S,trace
```

Subsequently, the rest of $S$ is filled as:
$S_{ij}=\max\left\{
\begin{array}{ll}
S_{i-1,j-1} & +d(a_i,b_j)\\
S_{i-1,j} & +d(a_i,-)\\
S_{i,j-1} & +d(-,b_j)
\end{array}
\right.$

This recursion is easily transformed into for-loops. We are free to select the order the matrix is filled in so we select to fill in row per row, i.e. rows become the inner loop, the columns the outer loop.

Again we keep track of the move that lead to a certain position by filling in the `trace` matrix.

```python
# Perform global alignment using the Needleman-Wunsch algorithm.
# Returns the aligned sequences (with gaps) and optionally prints the score matrix.
def global_align(seqA,seqB,print_dynamic_matrix = False):
    # Add 1 to lengths for DP matrix (to account for initial gap row/column)
    m, n = len(seqA)+1, len(seqB)+1
    S,trace = initiate_global_dp(m,n)
    # Fill in the rest of the dynamic programming matrix
    for i in range(1,m):
        for j in range(1,n):
            # Calculate scores for match/mismatch, deletion, and insertion
            # Note: i-1 and j-1 because Python uses 0-based indexing
            match = S[i-1,j-1] + match_score(seqA[i-1],seqB[j-1]) 
            delete = S[i-1,j] + match_score(seqA[i-1],'-') 
            insert = S[i,j-1] + match_score('-',seqB[j-1]) 
            # Choose the maximum score and record the move in trace
            S[i,j] = max(match,delete,insert)
            if match >= max(insert,delete):
                trace[i,j,:] = (-1,-1)  # Diagonal move: match/mismatch
            elif delete >= insert:
                trace[i,j,:] = (-1,0)   # Up: gap in seqB (deletion)
            else:
                trace[i,j,:] = (0,-1)   # Left: gap in seqA (insertion)
    if print_dynamic_matrix:
        print_dynamic(seqA,seqB,S,trace)  # Optional: print the score matrix
    return format_alignment(seqA,seqB,S,trace)  # Traceback to get alignment
```

Now everything is set. We can try the code for any of our favorite sequences. One can toggle the printout of the dynamic programming matrix by a boolean flag as a third argument.

```python
# Example: align two short DNA sequences and print the result
seqA,seqB = global_align("GATTA","GCTAC",True)
print_alignment(seqA,seqB)
```

I add a couple of extra alignments, check them manually as an excercise before the exam. Also try to run a couple of them yourself.

```python
# Another example: align two longer sequences
seqA,seqB = global_align("TGCATTA","GCATTAC",True)
print_alignment(seqA,seqB)
```

```python
# Try aligning even longer sequences for practice
seqA,seqB = global_align("CTATCTCGCTATCCA","CTACGCTATTTCA",True)
print_alignment(seqA,seqB)
```
