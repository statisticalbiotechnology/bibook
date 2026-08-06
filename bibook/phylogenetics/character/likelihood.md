---
authors:
  - name: Anders Andersson
---

# Maximum Likelihood

Maximum Likelihood (ML) asks which tree and model parameters make the observed alignment most probable.

## Likelihood is conditional on a model

For tree $T$, model parameters $\theta$ and sequence data $D$, the likelihood is

$$
L(T,\theta)=P(D\mid T,\theta).
$$

It is the probability of the data given the tree—not the probability that the tree is correct. Substitution-model parameters can describe base frequencies, transition and transversion rates, and rate variation among sites. Branch lengths describe expected substitutions per site.

ML selects the topology and parameter values that maximize $L$. In practice, programs maximize the log-likelihood because products of many small probabilities become sums:

$$
\log L_{total}=\sum_s \log L_s,
$$

where $L_s$ is the likelihood of alignment site $s$ under the usual assumption that sites evolve independently.

## Unknown ancestral states

Internal-node states are not observed. Rather than choosing a single reconstruction, likelihood calculations sum over every possible state. For a DNA node, the contribution is summed over A, C, G and T, weighted by equilibrium frequencies and transition probabilities along the branches.

This calculation can be performed efficiently with the pruning algorithm: likelihoods are computed from the leaves toward an internal reference point and reused instead of enumerating every complete ancestral assignment.

## Comparing trees

For each candidate topology, branch lengths and other model parameters are optimized. The optimized likelihoods are then compared. Calculating one tree's likelihood is tractable, but the number of possible topologies grows super-exponentially, so large analyses rely on heuristic searches through tree space.

## Strengths and limitations

ML can model multiple substitutions, unequal base frequencies and rate variation. It is statistically well founded and often more accurate than simple distance or parsimony approaches when the model is suitable. Its results can nevertheless be biased by poor alignment, inappropriate models, insufficient taxon sampling or failure of the heuristic search. It is also computationally demanding.

```{exercise}
Explain why a tree with a higher likelihood is not necessarily known to be the true tree. Identify at least two assumptions shared by all trees in an ML comparison that could be wrong.
```
