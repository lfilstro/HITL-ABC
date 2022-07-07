# Approximate Bayesian Computation with Domain Expert in the Loop

This repository contains the code used to run the experimental work of the article "Approximate Bayesian Computation with Domain Expert in the Loop", published at ICML 2022.

arXiv : https://arxiv.org/abs/2201.12090

## Instructions to run the experiments (Section 4)

### Experiment under model misspecification (Section 4.1)

Check the "misspecification" sub-folder. You will find there a notebook reproducing the experiment. Note that it uses some datasets, which are provided, as well as the R codes which were used to create them.

### Experiment in low-simulation regime (Section 4.2)

Check the "lowsim" sub-folder. The codes there use some datasets, which are provided, as well as the R code which was used to create them (gk_generateData.R)

Results of the HITL-ABC method are obtained by running
> python gk_hitl.py [nSim] [criterion]
where criterion can be "normal" or "random".

Results of Barnes' method is obtained by running
> python gk_Barnes.py [nSim]