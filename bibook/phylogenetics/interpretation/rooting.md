---
authors:
  - name: Anders Andersson
---

# Rooting Trees

Many phylogenetic inference methods produce an unrooted tree. Rooting gives the tree a direction through time and is therefore necessary for interpreting the order of evolutionary divergences.

## Rooted and unrooted trees

An unrooted tree describes splits among taxa but does not identify their common ancestor or the direction of evolution. Neighbor Joining, Maximum Parsimony and Maximum Likelihood commonly produce unrooted trees. Placing the root on different branches can imply different sequences of divergence without changing the unrooted relationships.

```{warning}
A taxon drawn at the left or bottom of an unrooted tree is not thereby ancestral or early-diverging. Those claims require a root.
```

## Outgroup rooting

An **outgroup** is known *a priori* to have diverged before the taxa of primary interest, the **ingroup**. The root is placed on the branch connecting the outgroup to the ingroup.

The outgroup must be outside the ingroup but sufficiently related for trustworthy alignment and phylogenetic signal. A very distant outgroup can introduce alignment ambiguity, substitution saturation and long-branch attraction. The outgroup can consist of one or several taxa.

Other rooting approaches include midpoint rooting and molecular-clock models. Midpoint rooting places the root halfway along the longest leaf-to-leaf path and assumes approximately clock-like (constant rate) evolution.

````{exercise}
Consider the following unrooted four-taxon topology:

```text
((A,B),(C,D))
```

This tree has five branches on which the root could be placed. Draw the five resulting rooted trees. For each tree, identify the first divergence from the root and describe which taxa form sister groups. Explain why the unrooted topology alone cannot determine which of these evolutionary histories is correct.
````
