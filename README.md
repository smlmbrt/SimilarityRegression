# SimilarityRegression

This is the code repositiory for Similarity Regression (SR) a method to predict motif similarity using weighted alignments. Description of the directories:

* `ConstructSimilarityModels/`: jupyter notebooks, and R scripts used to train and select SR models. This directory contains a `README` that describes the notebooks in greater detail.
* `Examples/` : Contains example data and a jupyter [notebook](https://github.com/smlmbrt/SimilarityRegression/blob/master/Example/Example%20Analysis%20Notebook.ipynb) with code to read TF gene/protein information from Cis-BP, and parse it into formats that can be used to train SR models, or score sequences using existing SR models.
* `Scripts/` contains python scripts and R code for aligning sequences to a Pfam HMM.
* `similarityregression/`: python module containing code to align DBDs, and score alignments using SR models.
* `CisBP/`: scripts to calculate E-score overlaps from data present in CisBP flat files.

`python` Dependancies: `numpy`, `pandas`, `biopython`, `sklearn`

`R` Dependancies: `caret`, `glmnet`, `PRROC`, `aphid`, `seqinr`

## Citation
Samuel A. Lambert, Ally Yang, Alexander Sasse, Gwendolyn Cowley, Mark X. Caddick, Quaid D. Morris, Matthew T. Weirauch, and Timothy R. Hughes (2019). **[Similarity Regression predicts evolution of transcription factor sequence specificity](https://www.nature.com/articles/s41588-019-0411-1)**. *Nature Genetics*. **51**:981–989.

**Abstract**

Transcription factor (TF) binding specificities (motifs) are essential for the analysis of gene regulation. Accurate prediction of TF motifs is critical, because it is infeasible to assay all TFs in all sequenced eukaryotic genomes. There is ongoing controversy regarding the degree of motif diversification among related species that is, in part, because of uncertainty in motif prediction methods. Here we describe **Similarity Regression**, a significantly improved method for predicting motifs, which we use to update and expand the Cis-BP database. Similarity regression inherently quantifies TF motif evolution, and shows that previous claims of near-complete conservation of motifs between human and Drosophila are inflated, with nearly half of the motifs in each species absent from the other, largely due to extensive divergence in C2H2 zinc finger proteins. We conclude that diversification in DNA-binding motifs is pervasive, and present a new tool and updated resource to study TF diversity and gene regulation across eukaryotes.
