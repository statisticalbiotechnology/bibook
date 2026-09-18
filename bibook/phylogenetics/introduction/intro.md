---
authors:
  - name: Anders Andersson
---

# Introduction to Phylogenetics

Phylogenetics uses heritable differences among organisms or genes to infer their evolutionary history. A phylogenetic tree is therefore a hypothesis of the past based on present-day observations, not a direct record of the past.

## Molecular evolution and speciation

First a few words about evolution. A population is a group of organisms of the same species occupying a particular space at a particular time (although there is some variation to the definition). Mutations continually introduce new variants (mainly substitutions, but also deletions, insertions and structural variants) into a population. **Natural selection** changes the frequencies of variants that affect reproductive success, while **genetic drift** changes frequencies by chance, especially in small populations. A new variant may disappear or become **fixed**, meaning that it eventually occurs in every member of the population.

When gene flow between two populations is reduced, mutations that arise may become fixed independently in each population. The populations consequently accumulate different substitutions. Selection can promote divergence when the populations experience different environments, while drift can produce divergence even without an adaptive difference. If reproductive isolation develops (the populations cannot cross breed and produce reproductive offspring), the populations form separate species.

Descendant species inherit most of their DNA from their common ancestral population. Shared and differing characters in modern sequences therefore contain an incomplete signal of their shared history. Phylogenetic methods use this signal to infer where evolutionary lineages diverged.

```{important}
A mutation is a change in a DNA molecule. A substitution is a mutation that has become fixed when comparing evolutionary lineages (e.g. different species). Phylogenetic models usually describe substitutions, although the two terms are sometimes used loosely.
```

## Molecular clocks
The molecular clock hypothesis states that substitutions accumulate at an approximately constant rate. Thus sequence distance can provide information about divergence time when a rate is known. The rate can be calibrated using fossils, dated geological events, ancient DNA or samples collected at known times. Substitution rates however differ among genes, sites and lineages because of generation time, selection, population size and molecular constraints. A universal molecular clock does therefor not exist.


## Phylogenies as inferred evolutionary histories

A **phylogeny** describes the evolutionary history of a group of replicating entities, such as species, viruses or genes. We normally observe only sequences from contemporary (today present) taxa. Common ancestors, ancestral sequences and divergence events are inferred from those observations. Phylogenetic divergence is assumed to be bifurcating, with each ancestral lineage splitting into two descendant lineages, making a tree a natural representation of evolutionary history.

Many different evolutionary histories could potentially explain the same observed data. A tree-building method uses explicit criteria—such as minimum evolutionary change or maximum likelihood—to select or compare hypotheses. The resulting tree should not be treated as certainty: its interpretation depends on the data, evolutionary model and statistical support.

## Parts of a phylogenetic tree

- **Leaves** or **terminal nodes** represent the sampled taxa or genes.
- **Internal nodes** represent inferred divergence events.
- A **branch** connects two nodes and represents an evolving lineage.
- The **root**, when known, represents the common ancestor of all observations and gives the direction of time.
- Two descendants sharing an immediate ancestor are often called **sister groups**.

An internal node is not one of the sampled descendants. For example, if A is sister to B, neither A nor B is the ancestor of the other.

## Topology and branch lengths

The **topology** is the branching order. A and B are more closely related to each other than either is to C only when A and B share a more recent common ancestor.

Branch lengths may contain additional information. In a **phylogram**, they are usually proportional to evolutionary change, commonly measured in expected substitutions per nucleotide or amino acid site. In a **cladogram**, only the topology is informative. Horizontal or vertical spacing chosen merely to make a tree readable has no biological meaning.

## Equivalent drawings and rotations

Rotating the descendants around an internal node does not change a topology. Thus `((A,B),(C,D))` describes the same relationships as `((B,A),(D,C))`. Compare trees by their groups of descendants, not by the left-to-right order of their leaves.

## Newick format

Newick is a compact, computer-readable representation of trees. Parentheses enclose descendants sharing a node, commas separate branches, and a semicolon ends the tree:

```text
((A,B),(C,D));
```

Branch lengths can be included and follow a colon:

```text
((A:0.10,B:0.10):0.08,(C:0.05,D:0.05):0.13);
```

Internal node labels, including support values, may appear after a closing parenthesis.

## Clades

A **clade** contains an ancestor and all of its descendants and is therefore monophyletic. A group assembled from several lineages while excluding some of the descendants of their most recent common ancestor, or some of its descendants, is not a clade. Such groupings may be polyphyletic or paraphyletic.

## Exercises

1. Explain why rotating branches around an internal node does not create a new tree.
2. In `((A,B),(C,D));`, which pairs are sister taxa? Is A more closely related to C than to D?
3. Write a Newick tree in which A and C are sister taxa and B is more distant to these.
4. Why is an inferred internal node not normally interpreted as one of the sampled sequences?
