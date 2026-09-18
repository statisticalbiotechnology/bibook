---
authors:
  - name: Anders Andersson
---

# Maximum Parsimony

Character-based methods infer phylogenetic relationships directly from the states of individual characters. In a sequence alignment, each column is treated as a character, and the nucleotide or amino acid observed in each sequence is its character state. Because the columns are evaluated separately, these methods retain information about which taxa share a state at each position rather than reducing every pair of sequences to a single distance.

Maximum Parsimony is a character-based method that favors the tree requiring the fewest changes in character states.

## The parsimony criterion

For a proposed topology, ancestral character states are assigned so that each alignment column requires as few changes as possible. The scores of all columns are added. A **most-parsimonious tree** has the smallest total number of changes.

This applies **Occam's razor**: prefer the hypothesis that explains the observations with the fewest assumptions—in phylogenetics, the tree requiring the fewest inferred substitutions. Maximum Parsimony therefore favors the shortest reconstructed evolutionary history, although the actual history may have involved additional, unobserved changes.

## A four-taxon example

For four taxa there are three possible unrooted tree topologies:

```text
((A,B),(C,D))   ((A,C),(B,D))   ((A,D),(B,C))
```

Suppose one alignment column is:

| Taxon | A | B | C | D |
|---|---:|---:|---:|---:|
| State | G | G | T | T |

The first topology explains the column with one change between the two groups. Each alternative requires at least two. Constant columns contribute no changes, and some variable columns do not distinguish among topologies. **Parsimony-informative** columns contain at least two states, each represented in at least two taxa.

## Searching tree space

The number of possible trees grows extremely rapidly with the number of taxa. Exhaustive evaluation is therefore feasible only for small datasets. Practical programs use heuristic searches: begin from one or more trees, rearrange branches, and retain changes that improve the score. A heuristic search may miss the global optimum, and several trees may tie for the best score.

## Strengths and limitations

Parsimony is intuitive, uses site-specific information and does not require an explicit probabilistic substitution model. However, it treats observed changes as a proxy for actual changes. Multiple substitutions can therefore cause underestimation.

A serious failure mode is **long-branch attraction**: rapidly evolving lineages may independently acquire similar states and be grouped together even though they are not closest relatives. Better taxon sampling and model-based methods can reduce this problem.

```{exercise}
For the four taxa A,B,C,D with the following character states of two nucleotides: A=A,A, B=G,G, C=A,T and D=T,G, determine which of the three four-taxon topologies has the smallest parsimony score for these sites. Draw one optimal assignment of ancestral states.
```
