---
authors:
  - name: Anders Andersson
---

# Neighbor Joining

Neighbor Joining (NJ) is a distance-based method that constructs an unrooted tree without assuming that all lineages evolve at the same rate.

## Why not simply join the closest pair?

The pair of taxa with the smallest evolutionary distance are not necessarily each other’s closest relatives. Their distance may be small because both are relatively slowly evolving. NJ corrects for this by considering each taxon's total distance from all other taxa.

For $n$ current clusters, define

$$
r_i=\sum_k d(i,k)
$$

and the joining criterion

$$
Q(i,j)=(n-2)d(i,j)-r_i-r_j.
$$

The pair with the smallest $Q$ value is joined. This favors pairs that are close to each other relative to their average divergence from the remaining taxa.

## Algorithm

1. Begin with one cluster per sequence and a distance matrix.
2. Calculate $r_i$ and the $Q$ matrix.
3. Join the pair with the smallest $Q(i,j)$.
4. Estimate the two branch lengths to their new node $u$.
5. Update distances using

   $$d(u,k)=\frac{d(i,k)+d(j,k)-d(i,j)}{2}.$$

6. Repeat until the tree is complete.

For a selected pair $i,j$, the branch from $i$ to $u$ is

$$
\delta(i,u)=\frac{1}{2}d(i,j)+\frac{r_i-r_j}{2(n-2)},
$$

and $δ(j,u)=d(i,j)-δ(i,u)$.

## UPGMA and NJ compared

| Property | UPGMA | Neighbor Joining |
|---|---|---|
| Output | Rooted | Unrooted |
| Strict molecular clock | Required | Not required |
| Joining rule | Smallest cluster distance | Smallest corrected $Q$ value |
| Leaf-to-root distances | Equal | May differ |
| Main advantage | Simplicity | Handles rate variation among lineages |

NJ is fast and often provides a useful starting tree for more computationally intensive methods. It still compresses each sequence pair into one distance and thus discards information about individual alignment columns.

```{exercise}
Explain why an unrooted NJ tree does not identify the oldest taxon. What additional information could be used to place a root?
```
