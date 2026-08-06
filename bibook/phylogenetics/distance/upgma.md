---
authors:
  - name: Anders Andersson
---

# UPGMA

**UPGMA** (Unweighted Pair Group Method with Arithmetic Mean) constructs a rooted tree by repeatedly joining the closest clusters. It is simple and useful for understanding distance-based inference, but its biological assumptions are restrictive.

## Agglomerative clustering

Initially, each sequence forms a cluster. UPGMA repeats four steps:

1. Find the two clusters with the smallest distance.
2. Join them at a new internal node.
3. Place that node at half their distance.
4. Replace the two clusters by their union and update its distances to all other clusters.

For clusters $A$ and $B$, UPGMA defines their distance as the arithmetic mean of all pairwise distances between their members:

$$
d(A,B)=\frac{1}{|A||B|}\sum_{i\in A}\sum_{j\in B}d(i,j).
$$

The word *unweighted* means that every original sequence contributes equally; clusters are therefore weighted by how many sequences they contain.

## Worked example

Consider this distance matrix:

| | S1 | S2 | S3 | S4 |
|---|---:|---:|---:|---:|
| S1 | 0.0 | 0.2 | 0.3 | 0.4 |
| S2 | 0.2 | 0.0 | 0.3 | 0.4 |
| S3 | 0.3 | 0.3 | 0.0 | 0.1 |
| S4 | 0.4 | 0.4 | 0.1 | 0.0 |

S3 and S4 join first at height $0.1/2=0.05$. Their average distances to S1 and S2 are both $(0.3+0.4)/2=0.35$.

S1 and S2 then join at height $0.2/2=0.10$. The distance between the two two-member clusters is

$$
\frac{0.3+0.4+0.3+0.4}{4}=0.35.
$$

They join at height $0.35/2=0.175$. The branches from the two internal nodes have lengths $0.175-0.10=0.075$ and $0.175-0.05=0.125$. The result is

```text
((S1:0.10,S2:0.10):0.075,(S3:0.05,S4:0.05):0.125);
```

## Ultrametric assumption

UPGMA produces an **ultrametric** tree: every leaf is equally distant from the root. Biologically, this corresponds to a strict molecular clock in which all sampled lineages have accumulated substitutions at the same rate.

If one lineage evolves much faster, UPGMA can cluster by similar observed distance rather than true ancestry and infer the wrong topology. Neighbor Joining was developed to avoid this strict-clock requirement.

## Strengths and limitations

UPGMA is fast, deterministic and easy to inspect. It is appropriate when distances are approximately ultrametric, and it is widely useful as a clustering method. For general phylogenetic inference, however, its equal-rate assumption must be justified.

```{exercise}
Repeat the first UPGMA step after changing $d(S3,S4)$ from 0.1 to 0.5. Which pair joins first, and at what height? Do not implement the algorithm in code; show the matrix reasoning.
```
