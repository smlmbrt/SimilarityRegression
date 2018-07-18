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
    
def ScoreAlignmentResult(resultDict, scoreDict, applyidenticalRule = True):
    #Score The Sequence
    if scoreDict['Model.Class'] == 'SequenceIdentity':
        Score = resultDict[scoreDict['Model.Name']]
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
    #Check if 100% identical (gets rid of proteins w/ truncations)
    if (applyidenticalRule == True) and (resultDict['PctID_L'] == 1):
        Classification = 'HSim'
    
    return(Score, Classification)

#Functions to iterate over Sequence Dictionaries seq : [(TF_ID, M_ID), ...] or [(TF_ID, P_ID), ...]
def SeqDictIterator_XoverY(mdict, pdict):
    for seq_m, ids_m in mdict.items():
        seq_m = seq_m.split(',')
        ids_m.sort()
        m = (ids_m, seq_m)
        for seq_p, ids_p in pdict.items():
            seq_p = seq_p.split(',')
            ids_p.sort()
            p = (ids_p, seq_p)
            yield(m, p)

def SeqDictIterator_XoverX(pdict):
    import itertools
    for seq_i, seq_j in itertools.combinations(pdict.keys(),2):
        ids_i = pdict[seq_i]
        ids_j = pdict[seq_j]
        i = (ids_i, seq_i.split(','))
        j = (ids_j, seq_j.split(','))
        yield(i, j)
        
#Function to map
from similarityregression import PairwiseAlignment as pwsaln
def AlignAndScore_DictPairs(t, OutputClasses = ['HSim', 'Dis']):
    OutputClasses = set(OutputClasses)
    #Unpack input
    i_ids = t[0][0]
    i_seq = t[0][1]
    
    j_ids = t[1][0]
    j_seq = t[1][1]
    
    #Align Sequences
    aln_result = pwsaln.AlignDBDArrays(('i', i_seq), ('j', j_seq))
    
    #Score alignment
    aln_score = srpred.ScoreAlignmentResult(resultDict=aln_result, scoreDict=SRModel)
    aln_Class = aln_score[1]
    
    #Parse 2 results list
    results = []
    if aln_Class in OutputClasses:
        for i in i_ids:
            for j in j_ids:
                if i < j:
                    results.append(list(i) + list(j) + list(aln_score))
                else:
                    results.append(list(j) + list(i) + list(aln_score))
    return(aln_Class, results)

def SeqDictIterator_ParseIdentical2Results(d):
    results = []
    for ids in d.values():
        for i, j in itertools.combinations(ids, 2):
            if i < j:
                results.append(list(i) + list(j) + [1, 'HSim'])
            else:
                results.append(list(j) + list(i) + [1, 'HSim'])
    return(results)