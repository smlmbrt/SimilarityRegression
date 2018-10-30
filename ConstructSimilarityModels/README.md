This README describes the notebooks required to parse data from CisBP into dataframes for training with R, and evaluating different SR models.

#### Creating SR Models:

* `Create DBD alignments and training dataframes.ipynb`: This notebook shows how the motif sequences are parsed from the CisBP flat files, and split into TF families for SR model training. 
* `SimilarityRegression_Train.R`: This Rscript takes a CisBP family ID as input and finds the relevant folder with SR training data in it. It then scales the data, and trains 4 SR models (8 if you use the smoothed data - this isn't included in the paper because the models aren't interpretable, and usually perform worse). The script generates a `Models/` folder in that contains the feature scaling (mean/sd), model coefficients, model predictions (on test and final models), PR/NPV statistics/score cutoffs, and PR curves (separated by regression and classification). There is also a helper script (`SimilarityRegression_HelperFunctions.R`) with extra functions used for model training. 

#### Selecting SR Models:

Once the SR models are trained the following notebooks are used to select the best SR model (on a variety of performance metrics):

* `Select Best SR Models | Parse SR Performance | Join SR Weights with DNAproDB contact frequencies .ipynb`: This notebook shows how the best SR model is selected and parsed into a scoring json file. It also shows the selection of the %ID thresholds for families that don't have enough data to train SR models. It also shows how these metrics are parsed into relative performances (SR vs. %ID), and joined with contact frequencies dervived from DNAproDB.
* `Evaluate Heldout Predictions.ipynb`: This notebook has the code for Figure 2D, and Figure S6.



