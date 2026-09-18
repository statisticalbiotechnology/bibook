---
authors:
  - name: Anders Andersson
---

# Bayesian Inference

Bayesian phylogenetic inference uses the same character data and substitution models as Maximum Likelihood (ML). In both approaches, the likelihood

$$
P(D\mid T,\theta)
$$

describes the probability of observing sequence data $D$ given a tree $T$ and model parameters $\theta$. The likelihood of each alignment site can be calculated with Felsenstein's pruning algorithm, exactly as described in the Maximum Likelihood chapter.

The two approaches use this likelihood differently. ML searches for the single tree and parameter values that maximize it:

$$
(\hat T,\hat\theta)_{ML}
=\underset{T,\theta}{\operatorname{argmax}}\;P(D\mid T,\theta).
$$

Bayesian inference instead uses the likelihood to update probability distributions over trees and parameters. It asks not only which tree scores best, but how posterior probability is distributed among alternative trees, branch lengths and model-parameter values.

## Bayes' theorem

Bayesian phylogenetics applies Bayes' theorem:

$$
P(T,\theta\mid D)
=\frac{P(D\mid T,\theta)P(T,\theta)}{P(D)}.
$$

The terms have distinct meanings:

- The **likelihood**, $P(D\mid T,\theta)$, measures how well a tree and its parameters explain the observed alignment.
- The **prior**, $P(T,\theta)$, represents assumptions about trees and parameter values before considering the current alignment.
- The **posterior**, $P(T,\theta\mid D)$, represents their relative probabilities after the prior has been updated by the data.
- The **marginal likelihood**, $P(D)$, is a normalizing constant that makes the posterior probabilities sum or integrate to one.

When only relative posterior probabilities are needed, the same relationship is often written

$$
P(T,\theta\mid D)
\propto P(D\mid T,\theta)P(T,\theta).
$$

Thus, **posterior $\propto$ likelihood $\times$ prior**.

## What does the prior describe?

Priors can be placed on all unknown parts of the analysis. A topology prior may assign equal prior probability to every allowed tree or may be based on a model of how lineages diversify. Branch-length priors describe which amounts of evolutionary change are plausible before examining the alignment. Other priors can describe base frequencies, substitution-rate parameters or variation in evolutionary rate among sites.

A prior is a probability distribution, not a single fixed guess. A broad prior allows many values, whereas a concentrated prior favors a narrower range. With a highly informative alignment, the likelihood will often dominate a reasonable prior. With limited or ambiguous data, the prior can have a stronger influence on the posterior. Priors should therefore be reported and, when their influence is uncertain, analyses can be repeated with alternative plausible priors.

### A simplified numerical example

Suppose, only for illustration, that two trees are being considered. Their likelihoods are

$$
P(D\mid T_1)=0.006
\qquad\text{and}\qquad
P(D\mid T_2)=0.004.
$$

If the trees have equal prior probabilities, their unnormalized posterior weights are

$$
0.006\times0.5=0.003
\qquad\text{and}\qquad
0.004\times0.5=0.002.
$$

Normalizing these two weights gives

$$
P(T_1\mid D)=\frac{0.003}{0.003+0.002}=0.6
$$

and $P(T_2\mid D)=0.4$. In this deliberately restricted example, $T_1$ has the higher posterior probability, but uncertainty remains: the analysis does not simply discard $T_2$.

```{note}
Real analyses consider many topologies and continuous model parameters. The numbers above illustrate Bayesian updating only; posterior probabilities in a real analysis are normalized over the complete set of trees and parameter values included by the model.
```

## Accounting for uncertain parameters

For each topology, ML normally reports the likelihood after branch lengths and other parameters have been optimized. Bayesian inference instead accounts for uncertainty in these parameters by **marginalizing** over them. The posterior probability of a topology is

$$
P(T\mid D)=\int P(T,\theta\mid D)\,d\theta.
$$

Conceptually, this adds the posterior contributions from all possible branch lengths and model-parameter values for that topology rather than relying only on their single best-fitting combination. Similarly, the posterior distribution of one branch length incorporates uncertainty in the topology and in the other parameters.

This gives an important contrast:

| Maximum Likelihood | Bayesian inference |
|---|---|
| Finds parameter values and a topology that maximize the likelihood | Estimates a posterior distribution over parameters and topologies |
| Conditions on the selected substitution model | Conditions on the selected substitution model **and** prior distributions |
| Commonly reports a best tree with bootstrap support | Commonly reports a consensus tree with posterior clade probabilities |
| Represents parameter uncertainty using methods such as confidence intervals or bootstrap replicates | Represents parameter uncertainty directly through posterior distributions and credible intervals |

Neither method removes uncertainty or guarantees that the inferred tree is correct. Both depend on the alignment, taxon sampling and evolutionary model.

## Sampling the posterior with MCMC

Calculating the denominator $P(D)$ directly would require summing over every possible tree and integrating over all possible parameter values. This is infeasible for realistic datasets. Bayesian phylogenetic programs therefore usually employ Markov chain Monte Carlo (**MCMC**) to sample from the posterior distribution.

An MCMC analysis starts from a tree and a set of parameter values. It then repeatedly proposes a change, such as rearranging part of the topology, modifying a branch length or changing a substitution-model parameter. Proposed states with higher posterior probability are readily accepted, but some lower-probability proposals are accepted as well. This allows the chain to move through tree space rather than becoming permanently trapped at the first local optimum it encounters.

After many iterations, and provided that the chain has converged, the frequency with which states are sampled approximates their posterior probability. The samples are correlated because each proposal begins from the current state, so the number of MCMC samples is not the same as the number of independent observations.

## Burn-in and convergence

Early MCMC samples may still reflect the arbitrary starting tree rather than the posterior distribution. An initial **burn-in** is therefore discarded. The remaining samples should not be interpreted until there is evidence that the chain has converged and explored the relevant regions of tree and parameter space.

Common diagnostics include:

- comparing independent runs started from different trees,
- inspecting traces of likelihoods and parameter values,
- checking that independent runs estimate similar clade frequencies, and
- examining effective sample sizes, which estimate how much independent information is contained in the correlated samples.

A long run is not automatically a good run. A chain can remain for a long time in only one region of tree space, which may give precise-looking but misleading summaries.

## Summarizing sampled trees

An MCMC analysis may retain thousands of sampled trees rather than producing only one tree. These samples are commonly summarized using a consensus tree. For a particular clade, its **posterior probability** is estimated by the fraction of retained trees containing that clade. If a clade occurs in 9,500 of 10,000 retained trees, its estimated posterior probability is

$$
\frac{9500}{10000}=0.95.
$$

This means that the clade has posterior probability 0.95 **under the specified alignment, likelihood model and priors**. It is not a model-independent guarantee that the clade is biologically correct. If the alignment or evolutionary model is seriously misleading, the posterior can confidently support the wrong relationship.

Continuous parameters are summarized by their posterior distributions. For example, a branch length can be reported using its posterior mean or median together with a **credible interval**, a range containing a specified proportion of its posterior probability.

## Posterior probabilities and bootstrap support

Posterior clade probabilities and bootstrap support are not interchangeable. They are produced by different procedures and answer different questions:

- A posterior probability is the probability assigned to a clade after combining the likelihood with the priors under a Bayesian model.
- Bootstrap support is the frequency with which a clade is recovered when alignment columns are resampled and the tree is inferred repeatedly.

Posterior probabilities are often numerically higher than bootstrap proportions, so the same numerical threshold should not automatically be interpreted in the same way for both measures.

## Strengths and limitations

Bayesian inference provides a coherent way to combine prior information with sequence evidence and to represent uncertainty across trees and model parameters. Instead of treating all non-optimal trees as irrelevant, it measures how posterior probability is distributed among alternatives. Marginalization also allows uncertainty in branch lengths and other parameters to be carried into summaries of topologies and clades.

These benefits come with costs and assumptions. Results may be sensitive to prior choices when the data contain little information. They can be biased by poor alignments, inadequate substitution models or limited taxon sampling, just as ML results can. MCMC analyses can also be computationally demanding, and failure to detect poor convergence can produce unreliable posterior summaries.

```{exercise}
1. Explain why Maximum Likelihood and Bayesian phylogenetic inference can use the same likelihood calculation but produce different kinds of results.
2. What information does a prior contribute to a Bayesian analysis? When would you expect it to have the greatest influence on the posterior?
3. Why should conclusions not be based on an MCMC run before its convergence has been assessed?
```
