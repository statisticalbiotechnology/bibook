---
authors:
  - name: Anders Andersson
---

# Evaluating Support

Tree inference and tree evaluation answer different questions. An optimality criterion (e.g. maximum likelihood or maximum parsimony) selects a preferred tree, while support analysis assesses which parts of the tree are stable under perturbation of the data and therefor relaiable.

## Bootstrap resampling

In statistics, **bootstrapping** is a resampling method used to examine how much a result might vary if the data had been sampled differently. Starting with a dataset of $n$ observations, many new datasets of size $n$ are created by (randomly) sampling observations from the original dataset **with replacement**. The statistic of interest is recalculated for every bootstrap replicate. Variation among the resulting estimates provides information about the sampling variability of the original result.

Sampling with replacement means that an observation can be selected more than once in a replicate, while another observation may not be selected at all. The bootstrap is especially useful when repeatedly collecting new independent datasets would be impractical.

```{note} Etymology
The statistical term *bootstrap* alludes to the expression [“pull oneself up by one's bootstraps”](https://en.wiktionary.org/wiki/pull_oneself_up_by_one%27s_bootstraps). The expression originally described an impossible task: lifting oneself by pulling on the straps of one's own boots. Statistical bootstrapping performs a similarly surprising feat. Without independently collected replicate datasets, it uses the observed dataset itself to construct pseudo-replicates and estimate sampling variability—figuratively pulling additional information from its own bootstraps.
```

The phylogenetic bootstrap examines sensitivity to the sampled alignment columns:

1. Sample $L$ columns with replacement from an alignment of length $L$.
2. Infer a tree from this pseudo-replicate alignment.
3. Repeat the procedure many times, often hundreds or thousands.
4. For each clade in the original tree, calculate the fraction of replicate trees containing the same split.

Because sampling is with replacement, a replicate contains some original columns several times and omits others.

## Interpreting support

A bootstrap value of 85% means that the corresponding split occurred in 85% of bootstrap replicate trees under the chosen analysis procedure. It is not an 85% probability that the clade is true. Support depends on the alignment, taxon sampling, model, inference method and number of replicates.

High support indicates consistency under column resampling, but systematic error can be consistently supported. Low support indicates that the data do not robustly distinguish that split; it should not automatically be interpreted as evidence for a particular alternative.

```{note}
When evaluating support for a clade, each replicate tree is checked for whether the same set of taxa forms a clade. The branching order within that clade does not need to be the same and does not affect its support value.
```

## Judging and displaying bootstrap support

There is no universal cutoff above which a clade should be regarded as reliable. As a commonly used rule of thumb for the standard phylogenetic bootstrap, values of 70% or more are often treated as reasonable support, while values of 90–95% or more are considered strong support. These boundaries are conventions rather than guarantees. Their interpretation depends on the data, inference method and bootstrap procedure, and even a highly supported clade can result from systematic error.

Bootstrap support is commonly displayed by writing the support value next to the corresponding internal branch or node in the reported tree. The number refers to the clade descending from that branch—or, in an unrooted tree, to the split defined by the branch. Values may be shown as percentages, such as 85, or as proportions, such as 0.85, so a figure legend should state which scale and support method are being used. Weakly supported branches are sometimes collapsed into polytomies or their values are omitted from the figure.

```{exercise}
An internal branch has 40% bootstrap support. Explain what this value means, what it does not mean, and give two possible reasons for the low support.
```
