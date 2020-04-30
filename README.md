# proDA-Paper

Currently available as a preprint:

Constantin Ahlmann-Eltze and Simon Anders: *proDA: Probabilistic Dropout Analysis for Identifying Differentially Abundant Proteins in Label-Free Mass Spectrometry*. [biorXiv 661496](http://www.biorxiv.org/content/10.1101/661496v1) (Jun 2019)


This repository contains the code to reproduce the figures for the paper describing the 
[proDA](https://github.com/const-ae/proDA) R package.

## Data

There are two datasets that are used for demonstration:

* Data on the phosphorylation dynamics from a study by Erik de Graaf et al.<sup>[1](#myfootnote1)</sup>
* Data studying the interaction landscape of Ubiquitin signalling by Xiaofei Zhang et al.<sup>[2](#myfootnote2)</sup>

Both can be found in the `data/` folder.

## Analysis

There are three additional folders that contain R markdown notebook that were used to generate the plots
for the paper:

* `approach_intuition` contains the code to give an overview of the ideas underlying `proDA`
    - [Mean-variance relation](https://htmlpreview.github.io/?https://github.com/const-ae/proDA-Paper/blob/master/approach_intuition/mean_variance_relation.nb.html)
    - [Location of missing values](https://htmlpreview.github.io/?https://github.com/const-ae/proDA-Paper/blob/master/approach_intuition/missing_value_location.nb.html)
    - [Probabilistic dropout model](https://htmlpreview.github.io/?https://github.com/const-ae/proDA-Paper/blob/master/approach_intuition/probabilistic_dropout_model.nb.html)
* `compare_performance` contains the code to run `DEP`, `QPROT`, `Perseus`, `DAPAR`, `EBRCT`,
   `Triqler` and `proDA` on the
  de Graaf data and make the validation and comparison plots
    - Null comparison on the de Graaf data set [notebook](http://htmlpreview.github.io/?https://github.com/const-ae/proDA-Paper/blob/master/compare_performance/null_comparison.nb.html)
    - de Graaf semi-synthetic dataset performance comparison [notebook](https://htmlpreview.github.io/?https://github.com/const-ae/proDA-Paper/blob/master/compare_performance/compare_performance.nb.html)
* `ubiquitination` contains the code that was used to analyze the Ubiquitination data
    - [Analysis notebook](https://htmlpreview.github.io/?https://github.com/const-ae/proDA-Paper/blob/master/ubiquitination/Ubiquitination_Analysis.nb.html)



## Sources

<a name="myfootname1">1.</a> de Graaf, E. L., Giansanti, P., Altelaar, A. F. M. & Heck, A. J. R. Single-step Enrichment by Ti4 + -IMAC and Label-free Quantitation Enables In-depth Monitoring of Phosphorylation Dynamics with High Reproducibility and Temporal Resolution . Mol. Cell. Proteomics 13, 2426–2434 (2014).

<a name="myfootname2">2.</a> Zhang, X. et al. An Interaction Landscape of Ubiquitin Signaling. Mol. Cell 65, 941–955.e8 (2017).
