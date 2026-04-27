# Introduction to Sequence Alignments

A **sequence alignment** arranges two or more DNA, RNA, or protein sequences so that corresponding positions appear in the same column, making regions of similarity and difference immediately visible. Alignments are the central operation of bioinformatics: they underpin gene annotation, protein function prediction, evolutionary analysis, and sequence database search.

## Illustration of a Sequence Alignment

```
GATTA-
GCT-AC
```

Each column represents one position in the alignment. Where both rows carry the same character the column is a **match** (G–G at position 1); different characters is a **mismatch** (A–C at position 2); a dash paired with a character is a **gap**, representing an insertion or deletion (indel) relative to the other sequence.

## What is needed to systematically construct a pairwise alignment?

To systematically construct a pairwise alignment we need three things.

- A scoring function, i.e. how is each position in the alignment scored.
- An alignment type, i.e. what type of alignment is to be constructed. 
- An alignment algorithm, i.e. how do we optimise the score of an alignment to form the desired alignment type

### Scoring function

A scoring function, $d(x,y)$, giving the score of a column of any characters $x$ and $y$.  The characters could be DNA or RNA bases of amino acids, for now, let’s think about them as DNA bases. A typical scoring function could be:

$$  
d(x,y)= 
  \begin{cases} 
  p & \text{if } x=y\\
  g & \text{if } x=- \text{ or } y=-\\
  n & \text{otherwise }
  \end{cases}
$$
  
Here, $p$, is called a match score, $n$, a mismatch score, and $g$ a gap penalty.

We score an alignment by summing up the scoring function over each position in an alignment. This allows us to evaluate if one alignment is better than another.

````{admonition} Example
:class: Note

Consider the alignment,
```
ACT  
A-T  
```
, using a match score of 2 and a gap penalty of -1. The scoring of this alignment is as follows:

- The first position has `A` aligned with `A`, which is a match. Assuming a match score of 2, this position contributes +2 to the total score.

- The second position has `C` aligned with a gap (`-`). With a gap penalty of -1, this position subtracts 1 from the total score.

- The third position has `T` aligned with `T`, another match. This contributes +2 to the total score.

Thus, the total score of the alignment `ACT` with `A-T` is calculated as, +2 -1 +2 = 3
````

### Exercises

````{exercise} Alignment Score Calculation for Given Alignment
:label: ex-alignscore

Given the following alignment:

```
Sequence 1: GATTATC
Sequence 2: GCTTG-C
```

Use the following scoring scheme:
- Match: +2
- Mismatch: 0
- Gap penalty: -2

**Questions:**
- What is the total alignment score for this alignment?

```{dropdown} **Reveal Answer**
Solution to [](#ex-alignscore)
To calculate the alignment score:

- Matches: G-G, T-T, T-T, C-C → 4 matches × +2 = +8
- Mismatches: A-C, A-G → 2 mismatches × 0 = 0
- Gaps: - → 1 gap × -2 = -2

Total alignment score = 8 - 0 - 2 = **6**
```
````

````{exercise} Another Alignment Score Calculation for Given Alignment
:label: ex-alignscore2

Given the following alignment:

```
Sequence 1: GATTACA
Sequence 2: GAT-GCA
```

Use the following scoring scheme:
- Match: +1
- Mismatch: -1
- Gap penalty: -3

**Question:**
- What is the total alignment score for this alignment with the specified gap penalty?

```{dropdown} **Reveal Answer**

To calculate the alignment score:

- Matches: G-G, A-A, T-T, C-C, A-A → 5 matches × +1 = +5
- Mismatches: A-G → 1 mismatch × -1 = -1
- Gaps: - → 1 gap × -3 = -3

Total alignment score = 5 - 1 - 3 = **1**
```
````

## Alignment Type

The alignment type specifies which part of each sequence must be covered by the optimal alignment.

**Global alignment** spans the full length of both sequences and penalises any unaligned ends. Use it when comparing closely related sequences of similar length — for example, orthologous genes from two species. The classic algorithm is [Needleman–Wunsch](needleman).

**Local alignment** finds the highest-scoring contiguous region of similarity anywhere within the two sequences, ignoring flanking regions entirely. Use it when only part of each sequence is expected to match — for example, detecting a conserved domain in otherwise divergent proteins. The classic algorithm is [Smith–Waterman](waterman).

**Semi-global alignment** requires one sequence to be fully covered while allowing unpenalised overhangs at either end of the other. Use it when a shorter sequence (e.g. a gene) should align entirely within a longer one (e.g. a genome). See [variants of Needleman–Wunsch](semi).

## Alignment Algorithm

Given the scoring function and alignment type, we have a definition of what we want to achieve, i.e., there is a definition of optimality. Now we come to the question of *how* we obtain such optimality.

### Exhaustive searches

The simplest imaginable strategy is to enumerate every possible alignment of two sequences, score each one, and return the best. This is called an **exhaustive search**. For short sequences it works: you can write out all alignments by hand. At each position in the alignment, the algorithm must decide whether to match/mismatch the next character from each sequence, insert a gap in sequence 1, or insert a gap in sequence 2 — three choices that can be made independently at every step. The total number of alignments therefore grows exponentially with sequence length.

To make this concrete, consider two sequences of length $n$. A rough lower bound on the number of distinct global alignments is $\binom{2n}{n}$ — the number of ways to interleave $n$ characters from each sequence before any gaps are added. For $n = 10$ this is already $\binom{20}{10} = 184{,}756$; for $n = 100$ it exceeds $10^{58}$. Biological sequences are routinely hundreds to thousands of residues long, placing an exhaustive search entirely beyond computational reach.

The exhaustive approach fails not because any individual alignment is hard to score, but because the *number of candidates* grows faster than any practical computer can handle. This is the fundamental problem that dynamic programming solves.

### Dynamic programming

Dynamic programming avoids the combinatorial explosion by recognising that the optimal alignment of two full sequences must contain within it the optimal alignment of every pair of their prefixes — a property called **optimal substructure**. Each sub-problem therefore only needs to be solved once.

In practice this means:

- **Fill a matrix.** Cell $(i, j)$ stores the best alignment score for the first $i$ characters of sequence 1 against the first $j$ characters of sequence 2. Each cell depends only on the three cells directly above, to the left, and diagonally above-left, so the whole table is filled in $O(mn)$ time.
- **Traceback.** Follow the recorded choices back from the end cell (global) or the highest-scoring cell (local) to read off the alignment.

This reduces an exponentially large search to a calculation that scales with the product of the sequence lengths, and it is the principle behind every alignment algorithm in the following chapters.

In the next chapters we will describe a set of such dynamic programming algorithms in detail, here is table summarizing them.

```{list-table} Sequence alignment types
:header-rows: 1

* - Alignment Type
  - Definition
  - Typical Use Case
  - Algorithm
* - **Global**
  - Aligns full length of both sequences
  - Closely related sequences of similar length
  - [Needleman–Wunsch](needleman)
* - **Local**
  - Aligns subsequences with highest similarity
  - Motif/domain detection, different-length sequences
  - [Smith–Waterman](waterman)
* - **Semi-global**
  - Aligns one sequence entirely, allows overhangs in the other
  - Gene-to-genome, substring matching
  - [Variants of Needleman–Wunsch](semi)
``` 

