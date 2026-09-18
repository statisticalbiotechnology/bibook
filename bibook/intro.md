---
authors:
  - name: Lukas Käll
  - name: Anders Andersson
---

# Bioinformatics for Biotechnology Students


[![Jupyter Book Badge](https://raw.githubusercontent.com/jupyter-book/jupyter-book/next/docs/media/images/badge.svg)](https://jupyterbook.org)
[![DOI](https://img.shields.io/badge/DOI-Cite_Jupyter_Book-blue)](https://jupyterbook.org/stable/cite/)

[**Download complete book as PDF**](exports/bibook.pdf)

This book is written with KTH's course [CB2442, Bioinformatics](https://www.kth.se/student/kurser/kurs/CB2442) in mind. Some of the material has been generated with ChatGPT.

```{tableofcontents}
```

## Run the notebooks

The notebooks in this book can be run directly on the KTH JupyterHub — no
installation needed. Sign in with your hub account, and the material is copied
into your own workspace. Following a link again later updates the material
**without** discarding the changes you have made.

[**Open the whole book in JupyterLab**](https://193.10.159.40.nip.io/hub/user-redirect/git-pull?repo=https://github.com/statisticalbiotechnology/bibook&urlpath=lab/tree/bibook/bibook&branch=main)

### Algorithms for Pairwise Alignments

- [Smith–Waterman in code](https://193.10.159.40.nip.io/hub/user-redirect/git-pull?repo=https://github.com/statisticalbiotechnology/bibook&urlpath=lab/tree/bibook/bibook/alignment/pairwise/sw_code.ipynb&branch=main)
- [Semi-global in code](https://193.10.159.40.nip.io/hub/user-redirect/git-pull?repo=https://github.com/statisticalbiotechnology/bibook&urlpath=lab/tree/bibook/bibook/alignment/pairwise/sg_code.ipynb&branch=main)

### Protein Sequence Alignment

- [Substitution matrices](https://193.10.159.40.nip.io/hub/user-redirect/git-pull?repo=https://github.com/statisticalbiotechnology/bibook&urlpath=lab/tree/bibook/bibook/alignment/protein/matrix.ipynb&branch=main)
- [Protein alignment in code](https://193.10.159.40.nip.io/hub/user-redirect/git-pull?repo=https://github.com/statisticalbiotechnology/bibook&urlpath=lab/tree/bibook/bibook/alignment/protein/prot_code.ipynb&branch=main)

### Multiple Sequence Alignments

- [Multiple sequence alignment](https://193.10.159.40.nip.io/hub/user-redirect/git-pull?repo=https://github.com/statisticalbiotechnology/bibook&urlpath=lab/tree/bibook/bibook/alignment/msa/msa.ipynb&branch=main)
- [Sequence logos](https://193.10.159.40.nip.io/hub/user-redirect/git-pull?repo=https://github.com/statisticalbiotechnology/bibook&urlpath=lab/tree/bibook/bibook/alignment/msa/seqlogo.ipynb&branch=main)
- [Viterbi for profile HMMs](https://193.10.159.40.nip.io/hub/user-redirect/git-pull?repo=https://github.com/statisticalbiotechnology/bibook&urlpath=lab/tree/bibook/bibook/alignment/msa/viterbi.ipynb&branch=main)

### Sequence Retrieval

- [BLAST in code](https://193.10.159.40.nip.io/hub/user-redirect/git-pull?repo=https://github.com/statisticalbiotechnology/bibook&urlpath=lab/tree/bibook/bibook/alignment/retrieval/blast_code.ipynb&branch=main)

:::{note}
The first launch takes a few seconds while your workspace is prepared. The
notebooks that fetch TCGA data will additionally download a few hundred
megabytes the first time they are run.
:::

<p xmlns:cc="http://creativecommons.org/ns#" >This work by <span property="cc:attributionName">Lukas Käll</span> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>
