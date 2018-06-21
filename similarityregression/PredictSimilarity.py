import numpy as np
import json

logistic = lambda x: 1 / (1 + np.exp(-x))

#Script to read the information required to score Similarity regression(SR) or Alignment %ID (PctID_L) models
def ReadSRModel(filename):
    with open(filename) as SRModel:
        srmodel = json.load(SRModel)
    if 'SR.Weights' in srmodel:
        srmodel['SR.Weights'] = np.asarray(srmodel['SR.Weights'])
    return(srmodel)
    
def ScoreAlignmentResult(resultDict, scoreDict):
    #Score The Sequence
    if scoreDict['ModelType'] == 'SequenceIdentity':
        Score = resultDict[scoreDict['ModelName']]
    else:
        weights = scoreDict['Weights']
        if scoreDict['ByPos'] == 'PctID':
            ByPos = resultDict['ByPos_PctID']
        else:
            ByPos = resultDict['ByPos_AvgB62']
            
        Score = scoreDict['Intercept'] + np.dot(weights, ByPos)
        
    #Check if it should be transformed re: logistic
    if scoreDict['ModelType'] == 'Classifier':
        Score = logistic(Score)
        
    #Classify Score Based on Thresholds
    if Score >= scoreDict['Threshold_0.75']:
        Classification = 'HSim'
    elif Score >= scoreDict['Threshold_0.25']:
        Classification = 'Amb'
    elif Score < scoreDict['Threshold_0.25']:
        Classification = 'Dis'
    else:
        Classification = 'NA'
    return(Score, Classification)