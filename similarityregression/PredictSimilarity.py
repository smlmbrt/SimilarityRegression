import numpy as np
import json

logistic = lambda x: 1 / (1 + np.exp(-x))

#Script to read the information required to score Similarity regression(SR) or Alignment %ID (PctID_L) models
def ReadSRModel(filename):
    with open(filename) as SRModel:
        srmodel = json.load(SRModel)
    #Convert to NP arrarys
    if 'SR.Weights' in srmodel:
        srmodel['SR.Weights'] = np.asarray(srmodel['SR.Weights'])
        srmodel['SR.FeatureScales.mean'] = np.asarray(srmodel['SR.FeatureScales.mean'])
        #Convert 0's to NAs
        sd = np.asarray(srmodel['SR.FeatureScales.sd'])
        sd[sd == 0] = np.nan
        srmodel['SR.FeatureScales.sd'] = sd
    #Check for Amb/Dis threshold
    if np.isnan(srmodel['Threshold.Dis']):
        srmodel['Threshold.Dis'] = None
    return(srmodel)
    
def ScoreAlignmentResult(resultDict, scoreDict):
    #Score The Sequence
    if scoreDict['Model.Class'] == 'SequenceIdentity':
        Score = resultDict[scoreDict['ModelName']]
    else:
        SRweights = scoreDict['SR.Weights']
        #Get postional scores
        ByPos = np.array(resultDict['ByPos.' + scoreDict['SR.Features'].replace('_','.')])
        #Normalize to features (f)
        f = (ByPos - scoreDict['SR.FeatureScales.mean'])/scoreDict['SR.FeatureScales.sd']
        f[np.isnan(f)] = 0 #Cleanup NAs
        Score = scoreDict['SR.Intercept'] + np.dot(SRweights, f)
        if scoreDict['SR.LogisticTransform'] == True:
            logistic = lambda x: 1 / (1 + np.exp(-x))
            Score = logistic(Score)
        
    #Classify Score Based on Thresholds
    ##Check if HSim/Amb
    Classification = np.nan
    if Score >= scoreDict['Threshold.HSim']:
        Classification = 'HSim'
    else:
        Classification = 'Amb'
    ##Check if Amb/Dis        
    if scoreDict['Threshold.Dis'] != None:
        if Score < scoreDict['Threshold.Dis']:
            Classification = 'Dis'
    return(Score, Classification)