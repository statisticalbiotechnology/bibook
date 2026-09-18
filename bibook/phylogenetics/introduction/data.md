---
authors:
  - name: Anders Andersson
---

# From Sequences to Phylogenetic Data

A robust phylogenetic analysis is more than just inference of the tree. The usual workflow is to collect homologous sequences, align homologous characters, infer a tree and evaluate its reliability.

## Homology

Characters compared in a phylogenetic analysis must be related by descent. Sequences sharing ancestry are **homologous**; homology is a relationship, not a percentage of similarity.

- **Orthologs** diverged through speciation and often reflect the species history.
- **Paralogs** diverged through gene duplication and reflect the history of gene copies.
- **Xenologs** have a history involving horizontal gene transfer.

If a paralog is mistaken for an ortholog, the inferred gene tree may disagree with the species tree even when the tree inference itself is correct.

## Gene trees and species trees

A sequence alignment normally informs a **gene tree**. It need not equal the species tree. Duplication, gene loss, horizontal transfer, recombination and incomplete lineage sorting can all produce genuine differences. Multiple carefully selected genes are therefore often used to infer organismal relationships.

Horizontal gene transfer is particularly important in prokaryotes. Conserved housekeeping genes (encoding essential functions in the cell) are frequently selected because many are transferred less often than genes involved in specialized functions, although no marker is universally immune.

## Choosing sequence data

The marker must vary at a suitable rate for the evolutionary timescale.

- Fast-evolving nucleotide regions can resolve recent divergences but may saturate over long times.
- Protein sequences retain useful signal across deeper divergences because several nucleotide changes can encode the same amino acid.
- Structural RNAs contain conserved regions useful for comparing distant organisms.
- Several genes can be concatenated to increase the number of characters, provided their histories and evolutionary properties are sufficiently compatible.

More data do not automatically remove systematic bias. Taxon sampling, contamination, sequence quality and model suitability also matter.

## Alignment defines the characters

Each column in a multiple sequence alignment is treated as a set of homologous characters inherited from a common ancestral position. An incorrect alignment may therefore creates false substitutions (or conservations) and can support a false tree.

Before tree inference, inspect the alignment for uncertain regions, excessive gaps, non-homologous sequence ends, frame disruptions and potential recombination. Ambiguous regions may need to be realigned, masked or removed.

```{exercise}
You want to infer the species relationships among four bacteria. For one species you accidentally select a distantly related paralog of the target gene. Predict how this could affect the tree, and explain why adding a more sophisticated substitution model would not solve the problem.
```
