# Introduction to Sequence Alignments

Sequence alignments form the cornerstone of all bioinformatics. It is just a slight stretch to claim that bioinformatics is the science of sequence alignments. Sequence alignment is the comparison of nucleotide or protein sequences to identify regions of similarity that may reveal functional, structural, or evolutionary relationships. This foundational technique has propelled our understanding of molecular biology, enabling the prediction of gene functions, the elucidation of evolutionary trajectories, and the bearing principle behind modern sequencing.

Sequence alignments serve as a bridge between raw genetic information and meaningful biological insights. They facilitate the comparison of newly sequenced genes or proteins against well-characterized databases, providing immediate context and understanding. Through these comparisons, researchers can infer the biological role of unknown genes, identify conserved domains crucial for protein function, and predict the evolutionary lineage of species.

Furthermore, sequence alignments play a pivotal role in advancing personalized medicine. By comparing an individual's genetic makeup with reference genomes, alignments help to pinpoint mutations linked to hereditary diseases or drug responses. This, in turn, guides the development of tailored treatments and preventative measures.

Sequence alignments are indispensable tools that not only advance fundamental scientific knowledge but also catalyze medical breakthroughs. In this chapter, we explore the techniques, applications, and impact of sequence alignments, laying the groundwork for a deeper understanding of their significance in bioinformatics.

## What is a Sequence Alignment?

### Definition of an Alignment

Sequence alignment is a method of arranging sequences of DNA, RNA, or proteins to identify regions of similarity that may be a consequence of functional, structural, or evolutionary relationships between the sequences. By comparing these sequences, alignments can reveal patterns that are not evident when analyzing single sequences in isolation.

### Importance of Sequence Alignments

Sequence alignments are invaluable in the study of biological sciences for several reasons:

- **Localization of Equivalent Regions:** Alignments allow scientists to identify corresponding regions among two or more sequences. This helps reveal conserved regions that often correspond to essential functional or structural elements, such as enzyme active sites or protein binding regions.

- **Quantification of Sequence Similarity:** Through alignments, it is possible to quantitatively measure how similar two or more sequences are. This similarity can be indicative of how closely related the sequences are in an evolutionary context, providing insights into the phylogenetic relationships between different organisms.

Understanding sequence alignments is crucial for students entering the field of biotechnology, as it lays the foundation for more advanced studies in modern sequencing and proteomics and in that capacity being fundamental for genetics, molecular biology, and evolutionary studies. In fact, alignments underpin nearly every aspect of modern biology.

### Illustration of a Sequence Alignment

To illustrate the concept of sequence alignment, consider the following example involving two DNA sequences:

```
GATTA-
GCT-AC
```

In this example, the alignment demonstrates a few key concepts:

- **Match:** A match occurs when the same nucleotide or amino acid is present in the same position in both sequences. For instance, the first nucleotide 'G' is a match.

- **Mismatch:** A mismatch is when different nucleotides (or amino acids) are aligned, suggesting a substitution event. For example, the nucleotide 'A' in sequence 1 and 'C' in sequence 2 at the second position are a mismatch.

- **Insertion and Deletion (Indels):** These occur when a nucleotide (or amino acid) is present in one sequence but not the other, resulting in gaps in the alignment. An insertion in sequence 2 (represented by the gap in sequence 1 after 'A') and a deletion in sequence 1 (indicated by the gap in sequence 2 after 'T') illustrate this concept.

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

If we want to find an optimal alignment of the full length sequences, we are searching a *global* alignment type. If we search the highest scoring stretch of an alignment, you should use a *local* alignment type. You can also use a semi-global alignment, searching for an optimal alignment, with the exception for any overshooting sequence terminals.

To expand on the alignment types, let's delve deeper into global, local, and semi-global alignments, their definitions, differences, and use cases:

### Global Alignment

#### Definition

Global alignment aligns two sequences in their entirety, from start to finish. It's particularly useful when the sequences are of similar length and you're interested in comparing them across their whole length.

#### Use Case

Global alignment is ideal for comparing closely related sequences to identify small differences, such as mutations in genes from individuals of the same species or closely related species. The Needleman-Wunsch algorithm is a classic method for performing global alignments.

### Local Alignment

#### Definition

Local alignment identifies the highest-scoring alignment for any substring of the compared sequences. Unlike global alignment, it focuses on regions of high similarity without forcing the alignment outside of these regions.

#### Use Case

Local alignments are best when the sequences under comparison are of different lengths or when only parts of the sequences are expected to be conserved. This approach is particularly useful in identifying functional domains within proteins or conserved motifs in DNA sequences. The Smith-Waterman algorithm is commonly used for local alignments.

### Semi-Global Alignment

#### Definition

Semi-global alignment, also known as glocal alignment, seeks an optimal alignment for one sequence within the other, allowing for overhangs at either end of the shorter sequence without penalty. This method combines aspects of both global and local alignments, ensuring alignment across the entire length of one sequence while allowing extensions in the other.

#### Use Case

Semi-global alignments are particularly useful when aligning sequences where one is a substring or a longer version of the other, such as in gene annotation. This method is beneficial when aligning a gene to a full genomic sequence, where the gene (shorter sequence) is expected to align entirely, but the genomic sequence (longer sequence) may have extensions (overhangs) that are not part of the gene.

Each alignment approach serves different purposes and is suited to particular scenarios in biological research. By understanding the nuances and applications of each method, researchers can select the most appropriate alignment strategy for their specific needs, whether they are comparing whole genomes, identifying conserved motifs, or annotating genes within longer sequences.

## Alignment Algorithm

Given the scoring function and alignment type, we have a definition of what we want to achieve, i.e., there is a definition of optimality. Now we come to the question of *how* we obtain such optimality.

### Exhaustive searches

Let's first consider why this is not as straightforward as you might first think. Given that computers are fast, why don't we just walk through all possible alignments that can be formed by two sequences and then select whatever alignment is optimal?

Combinatorial explosion refers to the exponential growth in the number of possibilities that need to be considered as the size of the problem increases. In the context of sequence alignment, this means assessing every possible way in which two sequences can be aligned by considering every possible insertion, deletion, and substitution.

For two sequences, the number of possible alignments grows exponentially with the length of these sequences. This is due to the fact that each position in the first sequence can be matched with any position in the second sequence or aligned with a gap, and vice versa. If you consider aligning a sequence of just 10 amino acids against another of the same length, the number of potential alignments, considering all possible gaps and matches, is astronomically high. Aligning sequences of 100 characters each could potentially result in evaluating more possibilities than there are atoms in the observable universe. The situation becomes untenably complex as the sequence lengths increase to typical biological lengths of hundreds or thousands of residues.

In conclusion, while exhaustive search methods might conceptually offer a way to ensure the discovery of the optimal alignment, their practical application is limited by the sheer number of combinations they generate. This makes them an impractical choice in the field of bioinformatics, where efficiency and speed are often as critical as accuracy.

### Dynamic programming

Dynamic programming is a powerful computational approach used in bioinformatics to efficiently compute sequence alignments without facing the combinatorial explosions typical of exhaustive searches. This method operates under the principle that the optimal alignment of two sequences can be derived by solving smaller, simpler subproblems. 

Each step in a dynamic programming approach to sequence alignment treats the calculation of an optimal alignment for a given pair of subsequences as independent of previous steps. This independence is key to the method's efficiency:

- **Decomposition**: The problem of aligning two sequences is broken down into aligning smaller subsequences. Each subproblem corresponds to finding the best alignment between prefixes of the original sequences.

- **Table Filling**: A matrix is used where each cell represents an optimal score for aligning the subsequences up to that point. The score in each cell is calculated based solely on its immediate predecessors (top, left, top-left diagonal), independent of how those scores were derived.

- **Pruning**: By filling the matrix systematically, dynamic programming avoids recalculating the alignment for every possible combination of sequence segments. Instead, it builds upon the optimal solutions of smaller subproblems. This significantly reduces the number of combinations that need to be considered.

- **Optimality**: Once the matrix is filled, the optimal alignment can be traced back from the cell representing the entire sequence (for global alignments) or from the highest scoring cell (for local alignments). This ensures that only the most relevant paths through the problem space are considered.

Dynamic programming thus allows bioinformatics algorithms to compute sequence alignments rapidly and efficiently, avoiding the computational infeasibility associated with an exhaustive search. This method ensures that each step is self-contained, using only the necessary data from directly related subproblems, thereby dramatically reducing the computational complexity and resource requirements.

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

