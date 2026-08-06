---
authors:
  - name: Anders Andersson
---

# Bayesian Phylogenetic Inference

Bayesian phylogenetics estimates a posterior distribution over trees and model parameters:

$$
P(T,\theta\mid D) \propto P(D\mid T,\theta)P(T,\theta).
$$

The likelihood $P(D\mid T,\theta)$ is combined with **prior distributions** on topologies, branch lengths and evolutionary-model parameters. The result describes their relative support after observing the data.

Because the number of trees is enormous, Markov chain Monte Carlo (**MCMC**) is used to sample trees and parameters approximately in proportion to their posterior probabilities. An initial **burn-in** is discarded, and the remaining samples are summarized as a consensus tree. The fraction of retained trees containing a clade is its posterior probability.

Posterior clade probabilities and bootstrap support are not interchangeable. A posterior probability follows from a specified model and priors; bootstrap support measures how consistently a clade is recovered after resampling alignment columns. Posterior probabilities are often numerically higher.

Bayesian results are trustworthy only when the chains have converged and explored the relevant tree space. Independent runs, trace inspection and effective sample sizes are therefore essential diagnostics.

```{exercise}
What information does a prior contribute to a Bayesian analysis, and why should conclusions not be based on an MCMC run before its convergence has been assessed?
```
