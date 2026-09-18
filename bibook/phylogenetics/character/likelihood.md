---
authors:
  - name: Anders Andersson
---

# Maximum Likelihood

Maximum Likelihood (ML) is a character-based method: it evaluates the pattern of nucleotides or amino acids in each alignment column. Unlike Maximum Parsimony, however, it does not simply count the smallest possible number of changes. It uses an explicit model of sequence evolution to calculate the probability of observing each character pattern on a proposed tree.

The central question is:

> With this tree and this model of sequence evolution, how probable would the observed sequences be?

The preferred tree is the one for which the observed data have the highest probability after the branch lengths and other model parameters have been optimized.

## Likelihood is conditional on a model

For tree $T$, model parameters $\theta$ and sequence data $D$, the likelihood is

$$
L(T,\theta)=P(D\mid T,\theta).
$$

It is important to read this expression in the correct direction. Likelihood is the probability of the data **given** the tree and model. It is not the probability that the tree is correct given the data. Thus, a likelihood of $0.01$ does not mean that the tree has a one-percent chance of being correct.

The model specifies how character states change along branches. Its parameters may describe equilibrium base frequencies, transition and transversion rates, and rate variation among sites. Branch lengths usually measure expected substitutions per site. Together, the model and branch lengths determine the probability of every possible change along every branch.

## Transition probabilities under the Jukes-Cantor model

The Jukes-Cantor model provides a simple example. It assumes that A, C, G and T have equal equilibrium frequencies and that all six types of nucleotide substitution occur at the same rate. Let $\alpha$ be the expected substitution rate per site per unit time and let $t$ be the elapsed time. The expected branch length is then $d=\alpha t$ substitutions per site.

If a site begins with nucleotide $i$, the probability that nucleotide $i$ is observed after time $t$ is

$$
P_{ii}(t)=\frac{1}{4}+\frac{3}{4}e^{-4\alpha t/3}
=\frac{1}{4}+\frac{3}{4}e^{-4d/3}.
$$

For a particular different nucleotide $j$, where $i\neq j$, the probability is

$$
P_{ij}(t)=\frac{1}{4}-\frac{1}{4}e^{-4\alpha t/3}
=\frac{1}{4}-\frac{1}{4}e^{-4d/3}.
$$

There are three possible nucleotides different from $i$, so the probability of observing **any** different nucleotide is

$$
P_{i\neq j}(t)=3P_{ij}(t)
=\frac{3}{4}\left(1-e^{-4d/3}\right).
$$

For example, consider a branch of length $d=0.1$ expected substitutions per site. Under Jukes-Cantor,

$$
P_{ii}\approx 0.9064
\qquad\text{and}\qquad
P_{ij}\approx 0.0312.
$$

Thus, if the ancestral state is A, the descendant is A with probability approximately $0.9064$, while each of C, G and T has probability approximately $0.0312$. Their probabilities sum to one. Notice that $P_{ii}$ is not simply the probability that no substitution occurred (which would be 0.9). A nucleotide can change and later return to its original state. Similarly, an observed difference can result from more than one substitution.

As $d$ becomes very large, both $P_{ii}$ and each $P_{ij}$ approach $1/4$. The descendant state then contains almost no information about the ancestral state: repeated substitutions have produced saturation.

```{note}
Sequence data generally identify the expected amount of evolutionary change along a branch, not substitution rate and elapsed time separately. In the simple case,

$$
d=\alpha t,
$$

where $\alpha$ is the substitution rate and $t$ is elapsed time. The same branch length $d$ could therefore result from a high rate over a short time or a low rate over a long time.

This ambiguity does not prevent ordinary phylogenetic tree reconstruction. The transition probabilities depend on their product $d$, so branch length can be estimated directly in expected substitutions per site. Estimating evolutionary rates or absolute divergence times separately requires additional assumptions or calibration information, such as a molecular clock, fossils or dated samples.
```

## Calculating the likelihood of one alignment site

### A tree with one branch

Begin with the simplest possible tree: two observed sequences connected by a single branch of length $d$. At one alignment position, one sequence has A and the other has T:

```text
A -------- T
     d
```

Under Jukes-Cantor, the probability of changing from A to the particular nucleotide T is

$$
P_{AT}(d)=\frac{1}{4}-\frac{1}{4}e^{-4d/3}.
$$

If we condition on the first nucleotide being A, the probability of observing T at the other end is simply

$$
P(T\mid A,d)=P_{AT}(d).
$$

For a phylogenetic site likelihood, however, the A at the first end is also part of the observed data. Its equilibrium probability must therefore be included, which under Jukes-Cantor is 1/4:

$$
\begin{aligned}
L(A,T\mid d)
&=\pi_A P_{AT}(d)\\
&=\frac{1}{4}\left(\frac{1}{4}-\frac{1}{4}e^{-4d/3}\right)\\
&=\frac{1}{16}\left(1-e^{-4d/3}\right).
\end{aligned}
$$

For example, if $d=0.1$, then

$$
P_{AT}(0.1)\approx 0.03121
$$

and

$$
L(A,T\mid d=0.1)
=\frac{1}{4}(0.03121)
\approx 0.007802.
$$

### A tree with three branches

Consider a slightly more complex yet simplified tree with one internal node connected to three observed sequences:

```text
                 A (Sequence 1)
                       |
                       |
                       o
                     /   \
                    /     \
      A (Sequence 2)       G (Sequence 3)
```

For this alignment column, the observed states at the tips are A, A and G. Assume that all three branches have length $d=0.1$ substitutions per site. Under the Jukes-Cantor model, we calculated above that $P_{ii}=0.9064$ for retaining the same nucleotide and $P_{ij}=0.0312$ for changing to a particular different nucleotide.

In real phylogenetic data, the nucleotide at the internal node is not observed. It could be A, C, G or T. Maximum Likelihood therefore calculates the probability of the observed tip states for every possible internal state and adds the results.

Under Jukes-Cantor, each possible internal nucleotide has equilibrium probability $1/4$. If the internal state is A, its contribution is

$$
\frac{1}{4}P_{AA}P_{AA}P_{AG}
=\frac{1}{4}(0.9064)^2(0.0312)
\approx 0.006409.
$$

If the internal state is G, two branches must change from G to A, while one retains G:

$$
\frac{1}{4}P_{GA}P_{GA}P_{GG}
=\frac{1}{4}(0.0312)^2(0.9064)
\approx 0.000221.
$$

If the internal state is C, all three branches require a particular change:

$$
\frac{1}{4}P_{CA}P_{CA}P_{CG}
=\frac{1}{4}(0.0312)^3
\approx 0.0000076.
$$

The contribution for an internal T is the same:

$$
\frac{1}{4}P_{TA}P_{TA}P_{TG}
\approx 0.0000076.
$$

The likelihood of the observed site is the sum of these four alternatives:

$$
\begin{aligned}
L_s
&=\sum_{i\in\{A,C,G,T\}}\pi_iP_{iA}P_{iA}P_{iG}\\
&\approx 0.006409+0.000221+0.0000076+0.0000076\\
&\approx 0.006645.
\end{aligned}
$$

The internal state A makes the largest contribution because it provides the most probable explanation of the observed pattern. Nevertheless, Maximum Likelihood does not simply declare A to be the ancestral state and discard the alternatives. It includes the contributions from all four possible states.

### Generalizing to a larger tree

For a larger tree, the same type of calculation is applied across all branches: the individual branch probabilities are multiplied to obtain the probability of one particular assignment of states to the internal nodes. Because these ancestral states are unknown, the calculation must be repeated for every possible combination of internal-node states, and the resulting probabilities must be summed. With four possible nucleotides at each of $n$ unknown internal nodes, there are $4^n$ possible assignments. This number grows rapidly, but **Felsenstein's pruning algorithm** avoids evaluating every complete assignment separately. Instead, the following calculation is repeated at each internal node.

For a node $v$ with two descendants connected by branches of lengths $t_1$ and $t_2$, define $L_v(i)$ as the likelihood of the observations below that node conditional on state $i$ at node $v$. If its descendants have conditional likelihoods $L_1$ and $L_2$, then

$$
L_v(i)=
\left[\sum_j P_{ij}(t_1)L_1(j)\right]
\left[\sum_k P_{ik}(t_2)L_2(k)\right].
$$

At a leaf, the conditional likelihood is one for the observed nucleotide and zero for the other three. For example, if a leaf contains G, then $L(G)=1$ and $L(A)=L(C)=L(T)=0$. Applying the equation repeatedly moves the calculation from the leaves toward an internal reference point. At that point, the four conditional likelihoods are weighted by the equilibrium base frequencies $\pi_i$:

$$
L_s=\sum_i \pi_i L_{root}(i).
$$

Under Jukes-Cantor, every $\pi_i=1/4$.

```{note}
The reference point used in the likelihood calculation does not necessarily represent the biological root. With the time-reversible substitution models commonly used for phylogenetic inference, the likelihood of an unrooted tree does not depend on where this computational root is placed.
```

## From one site to a complete alignment

The procedure above gives the likelihood $L_s$ of one alignment site. Under the usual assumption that sites evolve independently, the likelihood of the complete alignment is the product of the site likelihoods:

$$
L_{total}=\prod_s L_s.
$$

These probabilities rapidly become extremely small as more sites are included. Programs therefore use the log-likelihood, which converts the product into a sum:

$$
\log L_{total}=\sum_s \log L_s,
$$

where a larger value (one that is less negative) represents a better fit. For example, a log-likelihood of $-1200$ is higher than a log-likelihood of $-1250$.

## Comparing trees

For each candidate tree topology, the branch lengths (and other model parameters, if used) are adjusted to find the highest likelihood obtainable for that topology. The optimized likelihoods of the candidate trees are then compared. It would be unfair to compare a tree with well-optimized branch lengths against one with arbitrary branch lengths.

Calculating the likelihood of one tree is tractable, but the number of possible topologies grows extremely rapidly with the number of taxa. Large analyses therefore use heuristic searches. A program starts with one or more initial trees, for example obtained with Neighbor-Joining, rearranges its branches, optimizes model parameters and retains changes that improve the likelihood. As with a parsimony search, a heuristic likelihood search is not guaranteed to find the global optimum.

The resulting score is meaningful only relative to the specified data and model. Likelihood scores from different alignments should not normally be compared directly, and the tree with the highest likelihood among those examined is not thereby proven to be the true evolutionary tree.

## Strengths and limitations

The branch lengths used in the examples above are model parameters. In an actual ML analysis, they are generally not known beforehand: for each proposed topology, the program adjusts the branch lengths to find values that increase the likelihood. The substitution model can also contain additional parameters. More complex models may allow unequal nucleotide frequencies, different rates for different types of substitution, or variation in evolutionary rate among alignment sites. These parameters can likewise be estimated from the data by finding values that give a higher likelihood. ML analysis therefore involves more than comparing tree shapes; it compares trees together with their fitted branch lengths and evolutionary-model parameters.

This model-based framework makes its assumptions explicit and can account for multiple substitutions at the same site. It uses all site patterns rather than reducing an alignment to pairwise distances, and it sums over uncertain ancestral states rather than committing to one reconstruction.

These advantages do not make the method assumption-free. A poor alignment, an inappropriate substitution model, insufficient taxon sampling or a failed heuristic search can bias the result. Standard models also simplify biological evolution—for example, they commonly assume that sites are independent even though RNA structure, codons and protein interactions can make their evolution interdependent. Maximum Likelihood is also computationally more demanding than simple distance methods.

```{exercise}
1. For a Jukes-Cantor branch of length $d=0.2$, calculate the probability that an ancestral A is observed as A at the descendant node, the probability that it is observed specifically as G, and the probability that it is observed as any nucleotide other than A. Confirm that the four possible descendant-state probabilities sum to one.
2. Explain why a tree with a higher likelihood is not necessarily known to be the true tree. Identify at least two assumptions shared by all trees in an ML comparison that could be wrong.
```
