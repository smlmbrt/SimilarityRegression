def ReturnBlossum62Dict():
    from Bio.SubsMat.MatrixInfo import blosum62
    subMat = {}
    for p, v in blosum62.items():
        subMat[p] = v
        subMat[p[::-1]] = v
    return(subMat)

def Blossum62Score(pair, mat):
    try:
        return mat[pair]
    except:
        import itertools
        aas = set(itertools.chain(*mat.keys()))
        pair = list(pair)
        for i, aa in enumerate(pair):
            if aa not in aas:
                pair[i] = 'X' #Replace with UNKNOWN AA
        return mat[tuple(pair)]

def AlnmtPctID(aln_x, aln_y):
    match = 0
    total = 0
    for sx, sy in zip(aln_x, aln_y):
        if (sx != '-') and (sy != '-'):
            #print sx, sy
            total += 1
            if sx == sy:
                match += 1
    if total == 0:
        PctMatch = 0
    else:
        PctMatch = match/float(total)
    return(float(match), PctMatch)

def PercentIdentityVect(aln_l, aln_s, Norm = 'L', SmoothingWindow_3 = False, subMat = ReturnBlossum62Dict()):
    DBDAlignmentLength = len(aln_l[0])
    numsegments = float(len(aln_l))
    if Norm == 'S':
        numsegments = 0
        #This finds the number of segments that don't soley consist of gaps
        for seq in aln_s:
            seq = seq.replace('-','')
            if len(seq) > 0:
                numsegments += 1 
    import numpy as np
    AAPercIdentityVect = np.zeros(DBDAlignmentLength)
    AAAvgScoreVect = np.zeros(DBDAlignmentLength)
    #AAMaxScoreVect = [0]*DBDAlignmentLength
    
    if SmoothingWindow_3 == True:
        windowsize = 3.0
        #Implement P%ID Matching while looking #(flank) AA positions (non-gap) before and ahead of current sites
        for i in range(0, len(aln_l)):
            seg_L = aln_l[i]
            seg_S = aln_s[i]
            
            if (set(seg_L) == set(['-'])) or (set(seg_S) == set(['-'])):
                continue
                
            pos = 0
            LastNonGapPos = None
            
            for sx, sy in zip(seg_L, seg_S):
                if (sx != '-') or (sy != '-'):
                    window_ID = 0.0
                    window_Score = 0.0
                    
                    #Check LastNonGapPos
                    if LastNonGapPos != None:
                        L_LNG = seg_L[LastNonGapPos]
                        S_LNG = seg_S[LastNonGapPos]
                        if (L_LNG != '-') and (S_LNG != '-'):
                            window_Score += Blossum62Score((L_LNG,S_LNG), subMat)
                        if L_LNG == S_LNG:
                            window_ID += 1
                    LastNonGapPos = pos #Set for the next iteration
                    
                    #Check Current Position for matches
                    if (sx != '-') and (sy != '-'):
                        window_Score += Blossum62Score((sx,sy), subMat)
                    if sx == sy:
                        window_ID += 1
                            
                    #Check for next position for matches
                    nextpos = None
                    np = pos + 1
                    while (nextpos == None) and (np < len(seg_L)):
                        L_NP = seg_L[np]
                        S_NP = seg_S[np]
                        if (L_NP != '-') or (S_NP != '-'):
                            nextpos = np
                        np += 1
                    if nextpos != None:
                        L_NP = seg_L[nextpos]
                        S_NP = seg_S[nextpos]
                        if (L_NP != '-') and (S_NP != '-'):
                            window_Score += Blossum62Score((L_NP,S_NP), subMat)
                        if L_NP == S_NP:
                            window_ID += 1
                    
                    #Norm by SmoothingWindow size (hardcoded 3)
                    window_ID = window_ID/windowsize
                    window_Score = window_Score/windowsize
                    
                    AAPercIdentityVect[pos] += window_ID
                    AAAvgScoreVect[pos] += window_Score
                pos += 1
    else:
        for i in range(0, len(aln_l)):
            seg_L = aln_l[i]
            seg_S = aln_s[i]
            
            #Skip the %ID if one alignment is all gaps
            #This would add spurious hits with gap == gaps -> 1
            if (set(seg_L) == set(['-'])) or (set(seg_S) == set(['-'])):
                continue

            pos = 0
            for sx, sy in zip(seg_L, seg_S):
                #Check %ID
                if sx == sy:
                    AAPercIdentityVect[pos] += 1
                #Check Blossum63       
                if (sx != '-') and (sy != '-'):
                    posScore = Blossum62Score((sx,sy), subMat)
                    AAAvgScoreVect[pos] += posScore
                    #if posScore > AAMaxScoreVect[pos]:
                    #    AAMaxScoreVect[pos] = posScore
                pos += 1
                            
    #Normalize vectors by numsegments
    AAPercIdentityVect = AAPercIdentityVect/numsegments
    AAAvgScoreVect = AAAvgScoreVect/numsegments
    
    #Return the resulting P%ID vectors 
    return(AAPercIdentityVect.tolist(), AAAvgScoreVect.tolist()) # Deimplemented: AAMaxScoreVect

def ReturnLongerThenShorterArray(x, y):
    if len(x[1]) > len(y[1]):
        return(x, y)
    else:
        return(y, x)

def AlignDBDArrays(i_x, i_y, ByPosNorm = 'L'):
    # 1) Find the longer array
    L, S = ReturnLongerThenShorterArray(i_x, i_y)
    
    L_Info = L[0]
    L_Array = L[1]
    S_Info = S[0]
    S_Array = S[1]

    span = len(S_Array)
    DBDAlignmentLength = len(L_Array[0])
    AAlen_L = len(''.join(L_Array).replace('-',''))
    AAlen_S = len(''.join(S_Array).replace('-',''))

    # 2) Find best ungapped alignment without overhangs
    ## Store the alignment stats
    almnts_S = []
    IDs = []
    pctIDs_OverlappingSites = []
    
    ## Start the loop of the small DBD array over the larger array
    for start in range(0,len(L_Array)- span + 1):
        end = start + span
        s_pos = 0
        currentaln = []
        #Loop over the longer protein to add gaps if the shorter array is internal to the longer array
        for p in range(0,len(L_Array)):
            if start <= p < end:
                currentaln.append(S_Array[s_pos])
                s_pos += 1
            else:
                currentaln.append( DBDAlignmentLength*'-' )
        almnts_S.append(currentaln)
        numID, pctID_OverlappingSites = AlnmtPctID(''.join(L_Array), ''.join(currentaln))
        IDs.append(numID)
        pctIDs_OverlappingSites.append(pctID_OverlappingSites)

    ## Find the position of the best alignment
    i_BestAln = IDs.index(max(IDs))
    S_BestAln = almnts_S[i_BestAln]
    MultiAlnFlag = sum([x == max(IDs) for x in IDs]) > 1

    # 3) Prepare output dictionary for alnmt 
    o = {}
    o['ArrayLenDifference'] = len(L_Array) - span
    o['MultiAlnFlag'] = MultiAlnFlag
    o['L_Info'] = L_Info
    o['S_Info'] = S_Info
    o['PctID_O'] = pctIDs_OverlappingSites[i_BestAln]
    o['PctID_L'] = 1.0*IDs[i_BestAln]/AAlen_L
    o['PctID_S'] = 1.0*IDs[i_BestAln]/AAlen_S
    o['L_Array'] = L_Array
    o['S_BestAln'] = S_BestAln
    o['ByPosNorm'] = ByPosNorm
    o['i_BestAln'] = i_BestAln
    
    # 4) Measure the Positional %ID (P%ID)
    ## Normalizing by the Longest DBD
    o['ByPos.PctID'], o['ByPos.AvgB62'] = PercentIdentityVect(L_Array, S_BestAln, Norm = ByPosNorm, SmoothingWindow_3 = False)
    o['ByPos.PctID.Smooth3'], o['ByPos.AvgB62.Smooth3'] = PercentIdentityVect(L_Array, S_BestAln, Norm = ByPosNorm, SmoothingWindow_3 = True)
    
    return(o)

def CalculateGapFeatures(L, S_in, iOffset, norm = 'L'):
    #Create Gap alignment based on matchpos offset
    if len(L) == len(S_in):
        S = S_in
    else:
        S = [None]*len(L)
        for i, val in enumerate(S_in):
            S[iOffset + i] = val
    
    NumGaps = len(L[0].split('|'))
    
    ID_Identical = np.zeros(NumGaps)
    ID_HaveGap = np.zeros(NumGaps)  #Both have or both don't have
    ID_LenDiff = np.zeros(NumGaps) #Differences are negative
    #ID_SeqSim = np.zeros(NumGaps) 
    
    for l, s in zip(L, S):
        l = l.split('|')
        if s != None:
            s = s.split('|')
            for i, g in enumerate(zip(l, s)):
                g_l, g_s = g
                g_l = g_l.replace('-', '')
                g_s = g_s.replace('-', '')
                
                #ID_Identical
                if g_l == g_s:
                    ID_Identical[i] += 1
                
                #Calc ID_HaveGap
                HaveGap = 0
                if (len(g_l) == len(g_s)) or ((len(g_l) > 0) and (len(g_s) > 0)):
                    HaveGap = 1
                ID_HaveGap[i] += HaveGap
                
                #Calc len difference
                l_diff = abs(len(g_l) - len(g_s))*-1
                ID_LenDiff[i] += l_diff
                    
    #Calculate Normalized Gap Scores   
    if norm == 'L':
        divisor = len(L)
    elif norm == 'S':
        divisor = len(S_in)
    else:
        divisor = None
    
    ID_Identical = ID_Identical/divisor
    ID_HaveGap = ID_HaveGap/divisor
    ID_LenDiff = ID_LenDiff/divisor
    
    return(list(ID_Identical), list(ID_HaveGap), list(ID_LenDiff))