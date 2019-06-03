# proDA-Paper

This repository contains the code to reproduce the figures for the paper describing the 
[proDA](https://github.com/const-ae/proDA) package.

There are two datasets that are used for demonstration:
* Data on the phosphorylation dynamics from a study by Erik de Graaf et al.<sup>[1](#myfootnote1)</sup>
* Data studying the interaction landscape of Ubiquitin signalling by Xiaofei Zhang et al.<sup>[2](#myfootnote2)</sup>

Both can be found in the `data/` folder.

There are three additional folders that contain R markdown notebook that were used to generate the plots
for the paper:
* `approach_intuition` contains the code to give an overview of the ideas underlying `proDA`.
* `compare_performance` contains the code to run `DEP`, `QPROT`, `Perseus`, and `proDA` on the
  de Graaf data and make the validation and comparison plots
* `ubiquitination` contains the code that was used to analyze the Ubiquitination data




<a name="myfootname1">1.</a> de Graaf, E. L., Giansanti, P., Altelaar, A. F. M. & Heck, A. J. R. Single-step Enrichment by Ti4 + -IMAC and Label-free Quantitation Enables In-depth Monitoring of Phosphorylation Dynamics with High Reproducibility and Temporal Resolution . Mol. Cell. Proteomics 13, 2426–2434 (2014).

<a name="myfootname2">2.</a> Zhang, X. et al. An Interaction Landscape of Ubiquitin Signaling. Mol. Cell 65, 941–955.e8 (2017).
