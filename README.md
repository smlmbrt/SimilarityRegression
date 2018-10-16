# SimilarityRegression

This is the code repositiory for Similarity Regression (SR) a method to predict motif similarity using weighted alignments. Description of the directories:

* `ConstructSimilarityModels/`: jupyter notebooks used to train and select SR models. This directory contains a `README` that describes how to train SR models.
* `similarityregression/`: python module containing code to align DBDs, and score alignments using SR models.
* `CisBP/`: scripts to calculate E-score overlaps from data present in CisBP flat files.

`python` Dependancies: `numpy`, `pandas`, `biopython`, `sklearn`

`R` Dependancies: `caret`, `glmnet`, `PRROC`

## Citation
Samuel A. Lambert, Ally Yang, Alexander Sasse, Gwendolyn Cowley, Mark X. Caddick, Quaid D. Morris, Matthew T. Weirauch, and Timothy R. Hughes. **Similarity Regression predicts evolution of transcription factor sequence specificity**. *(In preparation)*

**Abstract**

Transcription factor (TF) binding specificities (motifs) are essential to the analysis of noncoding DNA and gene regulation. Accurate prediction of the sequence specificities of TFs is critical, because the hundreds of sequenced eukaryotic genomes encompass hundreds of thousands of TFs, and assaying each is currently infeasible. There is ongoing controversy regarding the efficacy of motif prediction methods, as well as the degree of motif diversification among related species. Here, we describe **Similarity Regression (SR)**, a significantly improved method for predicting motifs. We have updated and expanded the CisBP database using SR, and validate its predictive capacity with new data from diverse eukaryotic TFs. SR inherently quantifies TF motif evolution, and we show that previous claims of near-complete conservation of motifs between human and Drosophila are grossly inflated, with nearly half the motifs in each species absent from the other. We conclude that diversification in DNA binding motifs is pervasive, and present a new tool and updated resource to study TF diversity and gene regulation across eukaryotes.