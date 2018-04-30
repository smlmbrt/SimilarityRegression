#Requirements
from Bio import AlignIO
import itertools
from operator import itemgetter

def FastaIter(fileloc):
    # 1) Read the fasta file
    aln = AlignIO.read(fileloc, 'fasta')
    for record in aln:
        yield record.id, str(record.seq)

def ParseStockholmWithMatches(fileloc):
    # 1) Read the stockholm
    aln = AlignIO.read(fileloc, 'stockholm')
    
    # 2) find matchpos
    RF = ''
    with open(fileloc, 'r') as infile:
        for line in infile:
            if line.startswith('#=GC RF'):
                line = line.strip().split()
                RF += line[2]
    matchpos = []
    for i, v in enumerate(RF):
        if v == 'x':
            matchpos.append(i)
    
    #Name all the alignment positions
    PFamPos = list(RF)
    PFamPos[:RF.index('x')] = ['N-term']*RF.index('x')
    PFamPos[(RF.rfind('x') + 1):] = ['C-term']*(len(RF)- RF.rfind('x'))
    matchcount = 0
    for i, v in enumerate(PFamPos):
        if v == 'x':
            matchcount += 1
            PFamPos[i] = matchcount
    gPos = [i for i, v in enumerate(PFamPos) if v == '.']
    gPos = [map(itemgetter(1), g) for k, g in itertools.groupby(enumerate(gPos), lambda (i,x):i-x)]
    for l in gPos:
        name = 'GAP:'+ str(PFamPos[l[0]-1]) + '-'  + str(PFamPos[l[-1]+1])
        for i in l:
            PFamPos[i] = name
    return(aln, matchpos, RF, PFamPos)

        
def RFGapIntervals(RF, matchval = 'x'):
    gaps = []
    CurrentI = []
    for i, v in enumerate(RF):
        if (v == matchval) and (len(CurrentI) == 1):
            CurrentI.append(i)
            gaps.append(tuple(CurrentI))
            CurrentI = []
        if (v != matchval) and (len(CurrentI) == 0):
            CurrentI = [i]
        if (i == (len(RF)-1)) and (len(CurrentI) == 1):
            CurrentI.append(i + 1)
            gaps.append(tuple(CurrentI))
    return(gaps)