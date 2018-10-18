import sys
import os
import glob

#Functions
def ReadMatches(filename):
    matches = {}
    maxP = 57
    with open(filename, 'r') as r:
        d = r.read()
        d = d.split('>')[1:]
        for x in d:
            seq, m = x.strip().split('\n')
            s_pfam, s_seq, path = m.split('|')
            s_pfam = int(s_pfam)
            s_seq = int(s_seq)
            path = map(int, path.split(','))
            #PackDict 
            matches[seq] = (s_pfam, s_seq, path)
    return(matches)

def ParseMatches(matchdict, maxP):
    matchpos = {}
    for seq in matchdict.keys():
        s_pfam, s_seq, path = matchdict[seq]
        path_pfam = [None]*len(path)
        path_seq = [None]*len(path)

        #Find the first non-deletion (0) position in the path for sequence
        foundNonZero = False
        i = 0
        while foundNonZero == False:
            if path[i] != 0:
                path_seq[i] = s_seq
                foundNonZero = True
            i += 1
        ##Fill in the rest of the path 2 sequence relationship
        lastP = s_seq
        for i in range(path_seq.index(s_seq) + 1,len(path_seq)):
            #Advance sequence position if it's not a deletion
            if path[i] != 0:
                path_seq[i] = lastP + 1
                lastP += 1
        #Find the first non-insert (2) position in the path for Pfam
        foundNonInsert = False
        i = 0
        while foundNonInsert == False:
            if path[i] != 2:
                path_pfam[i] = s_pfam
                foundNonInsert = True
            i += 1
        #Fill in rest of PFam relationships
        lastP = s_pfam
        for i in range(path_pfam.index(s_pfam) + 1,len(path_seq)):
            #Advance Pfam position if it's not an insertion
            if path[i] != 2:
                path_pfam[i] = lastP + 1
                lastP += 1

        #Get Matches
        matches = ['-']*maxP
        for p_pfam, p_seq in zip(path_pfam, path_seq):
            #print p_pfam, p_seq
            if p_pfam != None and p_seq != None:
                matches[p_pfam - 1] = seq[p_seq - 1]
        matchpos[seq] = ''.join(matches)
    return(matchpos)

def HMM2MaxP(filename):
	with open(filename, 'r') as hmm:
		for line in hmm:
			line = line.strip().split()
			if line[0] == 'LENG':
				hmm_len = int(line[1])
	return(hmm_len)

#ARGS
HMM = sys.argv[1]
SEQS = sys.argv[2]
ALNTYPE = sys.argv[3]

#PrepOutput
loc_Output = 'DBDMatchPos_aphid/'
if os.path.isdir(loc_Output) == False:
    os.mkdir(loc_Output)

loc_RootName = loc_Output + SEQS.split('/')[-1].replace('.fa','')
#print loc_RootName
    
#1) Run R/Aphid Matches
os.system('Rscript RunAPHID.R %s %s %s %s'%(HMM, SEQS, ALNTYPE, loc_RootName)) 

#2) Parse Results
hmm_len = HMM2MaxP(HMM)
for aphidresultfile in glob.glob(loc_RootName + '*Viterbi*'):
    newname = aphidresultfile.replace('Viterbi', 'matchpos') + '.fa'
    #print aphidresultfile, newname
    aphid_matchpos = ParseMatches(ReadMatches(aphidresultfile), hmm_len)
    with open(newname, 'w') as outf:
    	for seq, matchpos in aphid_matchpos.items():
    		outf.write('>' + seq + '\n' + matchpos + '\n')