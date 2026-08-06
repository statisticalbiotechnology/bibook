---
authors:
  - name: Anders Andersson
---

# Maximum Parsimony

Maximum Parsimony is a character-based method: it evaluates the nucleotides or amino acids in individual alignment columns rather than reducing each sequence pair to one distance.

## The parsimony criterion

For a proposed topology, ancestral character states are assigned so that each alignment column requires as few changes as possible. The scores of all columns are added. A **most-parsimonious tree** has the smallest total number of changes.

This applies Occam's razor: prefer the history that explains the observations with the fewest evolutionary events. It does not claim that evolution always follows the shortest possible path.

## A four-taxon example

For four taxa there are three possible unrooted binary topologies:

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
For the character states A=A, B=G, C=A and D=G, determine which of the three four-taxon topologies has the smallest parsimony score for this site. Draw one optimal assignment of ancestral states.
```
