---
authors:
  - name: Anders Andersson
---

# Rooting Trees and Evaluating Support

Tree inference and tree evaluation answer different questions. An optimality criterion selects a preferred tree, while rooting gives it a direction through time and support analysis assesses which parts are stable under perturbation of the data.

## Rooted and unrooted trees

An unrooted tree describes splits among taxa but does not identify their common ancestor or the direction of evolution. Neighbor Joining, Maximum Parsimony and Maximum Likelihood commonly produce unrooted trees. Placing the root on different branches can imply different sequences of divergence without changing the unrooted relationships.

```{warning}
A taxon drawn at the left or bottom of an unrooted tree is not thereby ancestral or early-diverging. Those claims require a root.
```

## Outgroup rooting

An **outgroup** is known independently to have diverged before the taxa of primary interest, the **ingroup**. The root is placed on the branch connecting the outgroup to the ingroup.

The outgroup must be outside the ingroup but sufficiently related for trustworthy alignment and phylogenetic signal. A very distant outgroup can introduce alignment ambiguity, saturation and long-branch attraction. Whenever possible, use several appropriate outgroups and examine whether the root is stable.

Other rooting approaches include midpoint rooting and molecular-clock models. Midpoint rooting places the root halfway along the longest leaf-to-leaf path and assumes approximately clock-like evolution.

## Bootstrap resampling

The phylogenetic bootstrap examines sensitivity to the sampled alignment columns:

1. Sample $L$ columns with replacement from an alignment of length $L$.
2. Infer a tree from this pseudo-replicate alignment.
3. Repeat the procedure many times, often hundreds or thousands.
4. For each clade in the reported tree, calculate the fraction of replicate trees containing the same split.

Because sampling is with replacement, a replicate contains some original columns several times and omits others.

## Interpreting support

A bootstrap value of 85% means that the corresponding split occurred in 85% of bootstrap replicate trees under the chosen analysis procedure. It is not an 85% probability that the clade is true. Support depends on the alignment, taxon sampling, model, inference method and number of replicates.

High support indicates consistency under column resampling, but systematic error can be consistently supported. Low support indicates that the data do not robustly distinguish that split; it should not automatically be interpreted as evidence for a particular alternative.

## Reporting a tree responsibly

State the data and alignment, inference method, substitution model, rooting method, support procedure and number of replicates. Preserve branch-length units and avoid presenting unsupported nodes as established evolutionary events.

```{exercise}
An internal branch has 40% bootstrap support. Explain what this value means, what it does not mean, and give two possible reasons for the low support.
```
