---
authors:
  - name: Anders Andersson
---

# Evolutionary Distances

Distance-based methods summarize each pair of aligned sequences by a single number. The simplest number is the observed fraction of differences, but evolutionary models are needed when multiple changes can occur at the same site.

## The p-distance

For two aligned sequences, let $L$ be the number of compared positions and $m$ the number of mismatches. Positions containing gaps are commonly excluded. The **p-distance** is

$$
p = \frac{m}{L}.
$$

It is easy to interpret: if 8 of 100 compared sites differ, $p=0.08$. However, it measures observed differences rather than the number of evolutionary changes.

## Multiple substitutions and saturation

A site may change several times. A later substitution can overwrite an earlier one or restore the original state. Consequently, two identical bases need not have remained unchanged, and one observed mismatch can represent several substitutions.

As divergence increases, new substitutions increasingly strike previously changed sites. The observed p-distance then grows more slowly than the true number of substitutions. Eventually it approaches a ceiling—**substitution saturation**—and contains little information about additional evolutionary time.

## Jukes-Cantor correction

The Jukes-Cantor model assumes equal base frequencies and equal rates for all nucleotide substitutions. Under these assumptions, the corrected distance is

$$
d_{JC} = -\frac{3}{4}\ln\left(1-\frac{4p}{3}\right).
$$

The result is an estimate of substitutions per site. For small $p$, $d_{JC}\approx p$; the correction becomes larger as divergence increases. It is undefined when $p\geq0.75$, reflecting saturation under the model.

## Kimura two-parameter correction

Transitions exchange two purines ($A \leftrightarrow G$) or two pyrimidines ($C \leftrightarrow T$). Transversions exchange a purine and a pyrimidine. Transitions are often more frequent, so the Kimura two-parameter model distinguishes them.

Let $P$ be the fraction of compared sites with transitions and $Q$ the fraction with transversions. Then

$$
d_{K2P} = -\frac{1}{2}\ln(1-2P-Q)-\frac{1}{4}\ln(1-2Q).
$$

Like every correction, it depends on its assumptions and fails for sufficiently saturated data.

## Molecular clocks

The **molecular clock** hypothesis states that substitutions accumulate at an approximately constant rate. If each of two lineages evolves at rate $r$ substitutions per site per unit time after splitting $t$ time units ago, their expected distance is approximately

$$
d = 2rt.
$$

Thus sequence distance can provide information about divergence time when a rate is known. The rate can be calibrated using fossils, dated geological events, ancient DNA or samples collected at known times.

A mutation rate is the rate at which new mutations arise; a substitution rate is the rate at which changes become fixed. They are related but not generally identical. Rates can also differ among genes, sites and lineages because of generation time, selection, population size and molecular constraints. A strict clock is therefore an assumption to test, not a universal law. Relaxed-clock models allow rates to vary.

## From distances to a matrix

Calculating a distance for every pair of sequences produces a symmetric matrix with zeroes on its diagonal. UPGMA and Neighbor Joining use such a matrix, but they interpret it under different assumptions.

## Exercises

1. Two sequences differ at 12 of 100 ungapped positions. Calculate their p-distance.
2. Explain why p-distance increasingly underestimates evolutionary change as sequences diverge.
3. Why might two genes from the same species pair imply different molecular-clock dates?
4. What biological assumptions distinguish Jukes-Cantor from Kimura's two-parameter model?
