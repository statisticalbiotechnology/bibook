# Multiple Sequence Alignments

## Intoduction 

While pairwise alignments enable comparisons between two sequences, many biological analyses require a more comprehensive approach: the alignment of multiple sequences simultaneously. Multiple sequence alignments (MSAs) offer a holistic view of sequence similarities and differences, illuminating functional and evolutionary relationships across diverse genomes and proteomes.

### Exponential Time Complexity

In principle, one can use dynamic programming to form MSAs. We can archive such systems using the same ideas as for the pairwise alignments, but just expanding the dynamic programming matrix into a tensor of the same dimension as the number of sequences.  However, a major drawback is its exponential time complexity.  Aligning $T$ sequences each of length $N$ can result in a time complexity of $O(N^T)$, making it impractical for datasets with more than a few sequences.

### Iterative Approaches as a Practical Alternative

Given the inefficiencies of dynamic programming for MSAs, bioinformaticians have turned to heuristics that mimic the behaviour of full dynamic programming. Particularly iterative approaches, such as progressive alignment. These methods simplify the alignment process by breaking it down into manageable steps. Progressive alignment, for instance, first generates a guide tree based on pairwise similarities, then iteratively aligns sequences according to the tree's order. This reduces the computational burden and allows for the processing of large datasets in a reasonable time frame.

## Progressive Alignments

While dynamic programming offers an optimal solution for pairwise alignments, its inefficiency for MSAs necessitates a balance between accuracy and computational efficiency. Iterative and heuristic methods like progressive alignment and consistency-based algorithms strike this balance, providing near-optimal solutions in significantly less time.

A key method for constructing MSAs is progressive alignment, a stepwise approach that builds upon pairwise alignments to create a broader, more nuanced comparison. The progressive alignment method begins by generating a distance matrix, which quantifies the pairwise similarity between all sequences in the dataset. From this matrix, a guide tree is constructed, representing the evolutionary relationships inferred from the sequence similarities.

Progressive alignment then follows this guide tree, iteratively aligning sequences in order of their proximity on the tree. This process begins with the closest pair of sequences, which are aligned directly. Subsequent sequences are added progressively, with each new sequence being aligned to the growing consensus. This stepwise process creates an alignment that reflects both sequence similarities and evolutionary relationships.

The progressive alignment method offers several advantages. It efficiently handles large datasets, making it well-suited for genomics and proteomics research. By utilizing an evolutionary framework, it provides deeper insights into the conserved and divergent features of sequences, aiding in the identification of functionally important regions and evolutionary patterns.

In this chapter, we delve into the mechanics of progressive alignment, exploring its implementation and applications in bioinformatics. By understanding this approach, we gain valuable insights into the biological relationships and evolutionary trajectories that shape our world.

### Progressive Alignment Steps

The process builds alignments in a stepwise manner, starting with the most similar sequences and gradually incorporating others. Here's a detailed look at how to create a progressive alignment. Lets assume you start with $T$ sequences, $C=\{A_1, \ldots, A_T\}$.

1. Compute Pairwise Alignments: Initially, calculate all pairwise alignments between each pair of sequences in $C$ using a suitable algorithm such as Needleman-Wunsch or Smith-Waterman. This provides a matrix of scores or distances that indicate how similar or different each pair of sequences is.

2. Use the pairwise distances to identify the two sequences in $C$, $A_i$ and  $A_j$ that are closest (those with the highest similarity score or lowest distance), and form an alignment of those two sequences.  Remove $A_i$ and $A_j$ from $C$. We call the alignment of those two sequences $B$.

3. Add $B$ to $C$. Set the distance from $B$ to each of the other sequences in $C$ to the average of their distance from $A_i$ and $A_j$.

4. If there are more sequences than one in $C$, then go to step 2.

### Progressive Alignment Steps

::::{grid} 2 2 4 4
:class-container: text-center
:gutter: 1

:::{grid-item-card}
Input sequences
```none
IAMAPEPTIDE  
IAMPEPTIDE
IAMPEPPED
IAMAPEPPERD
```
These are the sequences we want to align.
:::

:::{grid-item-card}
Seq1 and Seq2
```none
IAMAPEPTIDE
IAM-PEPTIDE
```
We form a pairwise alignment of those two most similar sequences.
We then pool the resulting sequences with the other two sequences.
:::

:::{grid-item-card}
Seq3 and Seq4.
```none
IAM-PEPPE-D
IAMAPEPPERD
```
These sequences were the most similar among remaining ones.
:::

:::{grid-item-card}
(Seq1,Seq2) and (Seq3,Seq4).
```none
IAMAPEPTID-E
IAM-PEPTID-E
IAM-PEPP-E-D
IAMAPEPP-ERD
```
There are two generalized sequences left. We align those two, and we are then done.
:::
::::

### Applications of Progressive Alignments

Modern bioinformatics tools, such as MUSCLE {cite}`edgar2004muscle`, MAFFT {cite}`katoh2002mafft`, T-Coffee {cite}`notredame2000t`, and Kalign {cite}`lassmann2005kalign` incorporate iterative refinement as a core feature, allowing researchers to achieve high-quality alignments with minimal manual intervention. These tools use sophisticated algorithms to iteratively refine alignments, often integrating additional information such as structural data or evolutionary constraints. This integration enhances the overall robustness of the alignment process, making these tools highly effective for a wide range of applications in genomics, proteomics, and evolutionary biology.

## Iterative Refinement of MSAs

While progressive alignment provides a structured approach to generating multiple sequence alignments (MSAs), it is often just the starting point for more refined alignments. The initial alignment produced by progressive methods can be improved through iterative refinement, a crucial step for enhancing alignment accuracy and resolving inconsistencies that may arise during the initial process.

### The Need for Iterative Refinement

Progressive alignment can introduce biases, particularly when the initial pairwise alignments are not optimal or when sequences with varying evolutionary distances are involved. These biases can propagate through the alignment process, leading to suboptimal alignments, especially in regions where the evolutionary signals are weak or ambiguous. Iterative refinement addresses these issues by revisiting and adjusting the alignment, ensuring that it more accurately reflects the true evolutionary relationships and functional similarities among the sequences.

### Iterative Refinement Techniques

There are several approaches to iterative refinement, each with its own strengths and applications:

* Realignment of Subsets: One common method involves selecting subsets of sequences or regions of the alignment and realigning them independently. This approach is particularly useful for correcting local misalignments that occur due to the progressive alignment process. By focusing on smaller sections of the alignment, the algorithm can correct errors without the computational burden of realigning the entire sequence set.
* Profile-based Realignment: In this technique, the MSA is treated as a series of profiles, where each profile represents a group of already aligned sequences. The iterative refinement process involves realigning these profiles against each other or introducing new sequences into the existing profiles. This method is effective in maintaining the overall structure of the alignment while correcting specific errors.
* Consistency-based Refinement: This approach leverages consistency information, which ensures that the aligned pairs in the final MSA are consistent with the pairwise alignments. By re-evaluating the alignment through the lens of consistency, this method can reduce the occurrence of alignment errors that arise from conflicts between pairwise comparisons and the global alignment.
* Iterative Re-scoring and Re-alignment: Another approach involves re-scoring the alignment based on a specific objective function (such as a scoring matrix or evolutionary model) and then realigning the sequences to maximize this score. This iterative process continues until the alignment stabilizes and no further improvements can be made. This method ensures that the alignment converges to a local optimum that best fits the given scoring criteria.

### Applications and Advantages of Iterative Refinement

Iterative refinement is widely used in bioinformatics tools such as MUSCLE {cite}`edgar2004muscle` and MAFFT {cite}`katoh2002mafft`, which implement various refinement strategies to enhance MSA quality. The iterative process ensures that the final alignment is more accurate and reliable, making it particularly valuable for downstream analyses, such as phylogenetic tree construction, structural modeling, and functional annotation.

Moreover, iterative refinement is adaptable, allowing bioinformaticians to tailor the process based on the specific characteristics of the sequences being aligned. This flexibility makes it an indispensable tool in the alignment of sequences with complex evolutionary histories, such as those involving duplications, deletions, or recombination events.

## References

```{bibliography}
:filter: docname in docnames
```
