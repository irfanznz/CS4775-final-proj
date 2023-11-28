# CS4775-final-proj

## Introduction

This is a project where we try to reimplement a profile-HMM based multiple sequence alignment program.

## NEED TO DO:

### Model-related stuff

- Fix correctness of emission probabilities for match states
- Build a transition probabilities matrix from transition counts (probably refactor transition counts function too to include begin and end states according to HMM topology)
- Get background probabilities of amino acids for use as emission probabilities for insert states
- Convert probabilities into log form

### MSA-related stuff

- Implement Viterbi algorithm for finding the most likely path the sequence takes in the HMM
- This corresponds to the alignment of the sequence(s) to the model, so we're basically done here

### Benchmarking

- Use BAliBASE to benchmark our alignments (mostly for use in slides, pretty much plug and play BAliBASE has a scoring program already we can compare our alignments to that of MUSCLE or something)
