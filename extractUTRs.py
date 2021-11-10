import sys
import pandas as pd
import itertools
import more_itertools as mit
from collections import Counter
import numpy as np
import collections
from collections import OrderedDict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--assembledGTF", type=str, help="assembled transcripts in GTF format")
parser.add_argument("-c", "--annotationsGTF", type=str, help="CDS annotations in GTF format")
parser.add_argument("-f", "--RATTLE_FASTQ", type=str, help="RATTLE output 'transcriptome.fq' in FASTQ format")
parser.add_argument("-a", "--assembler", type=str, help="program used to assemble transcripts <Stringtie2> or <RATTLE>")
parser.add_argument("-o", "--output", type=str, help="output directory (default: current folder)")
args = parser.parse_args()

if args.assembledGTF == None or args.annotationsGTF == None or args.assembler == None:
    print('extractUTRs.py  unable to start:   The fields -t , -c and -a must be specified...')
    sys.exit()

if args.assembler == "RATTLE" and args.RATTLE_FASTQ == None:
    print('extractUTRs.py  unable to start:   The field -f must be specified...')
    sys.exit()

GTF1 = args.assembledGTF
GTF2 = args.annotationsGTF
RATTLE_FASTQ = args.RATTLE_FASTQ

if (args.output):
    outputDir=args.output
    if (outputDir[-1] != "/"):
        outputDir= outputDir+"/"
else:
    outputDir="./"


f = open(GTF1,'r')
counter = 0
for line in f:
    if line.startswith('#'):
        counter += 1

if counter == 0:
    assemblyGTF = pd.read_csv(GTF1, sep="\t", header=None)
else:
    assemblyGTF = pd.read_csv(GTF1, skiprows = list(range(counter)) , sep="\t", header=None)

# Filter out twin transcripts with lower coverage
df = assemblyGTF.iloc[:,0:7].duplicated(subset= None,keep=False)
a  = [i for i, x in enumerate(df) if x] ; b = [list(group) for group in mit.consecutive_groups(sorted(a))]

c = []
for i in range(0,len(b)):
    duplic = assemblyGTF.iloc[b[i],0:7]
    cut = None
    zerocut = False
    for j in range(0,len(b[i])-1):
        if any(duplic.iloc[j] != duplic.iloc[j+1]):
            if j == 0:
                c.append([b[i][0]]) ; zerocut = True
                continue
            elif j != 0 and cut == None:
                c.append(b[i][:j+1]) ; cut = j
                continue
            elif j != 0 and cut != None:
                c.append(b[i][cut+1:j+1]) ; cut = j
    if zerocut == True and cut == None:
        b[i].remove(b[i][0])
    if cut != None:
        b[i] = b[i][cut+1:]

b = b+c ; b.sort()

# correction of eventual error in grouping redundant transcripts with same length   
toCorrect = []
counter = -1
for i in b:
    counter += 1
    if len(i) > 1:
        for x in range(len(i)):
            for y in range(x + 1, len(i)):
                if all(assemblyGTF.iloc[i[x],0:7] == assemblyGTF.iloc[i[y],0:7]) == False:
                    toCorrect += [[counter,[i[x],i[y]]]]      
                    

unlist_toCorrect = list(itertools.chain.from_iterable(toCorrect)) 
concatenate_unlist_toCorrect = [x for x in unlist_toCorrect if not isinstance(x,int)]  
concatenate_unlist_toCorrect = list(itertools.chain.from_iterable(concatenate_unlist_toCorrect)) 
duplicates = [item for item,count in collections.Counter(concatenate_unlist_toCorrect).items() if count > 1] 
split_b = []
toRemove = []
for i in toCorrect:  # find element to be split off from b lists
    for j in range(len(i)-1):
        if i[1][j] in duplicates and [j] in split_b:
            continue
        if i[1][j] in duplicates:
            split_b.append([i[1][j]]) ; toRemove += [i[0]] 
        if i[1][j] not in duplicates:
            split_b.append([i[1][j]]) ; toRemove += [i[0]] 

split_b = list(itertools.chain.from_iterable(split_b)) 

for i in range(len(toRemove)):
    if split_b[i] in b[toRemove[i]]:
        b[toRemove[i]].remove(split_b[i]) 
    else:
        continue

counter = 0
for f in b:
    if len(f) == 1:
        counter += 1


tmp = []
counter_tmp = -1
for i in range(0,len(b)):
    if len(b[i]) > 1:
        tmp += [b[i]]
        counter_tmp += 1
    if len(b[i]) == 1:    
        bIndeces = list(range(0,len(b))) ; bIndeces.remove(i)
        for j in bIndeces:
            if all(assemblyGTF.iloc[b[i][0],0:7] == assemblyGTF.iloc[b[j][0],0:7]):
                counter_tmp += 1
                tmp += [b[j]] ; tmp[counter_tmp] += [b[i][0]]
                continue
            if all(assemblyGTF.iloc[b[i][0],0:7] != assemblyGTF.iloc[b[j][0],0:7]) and j == bIndeces[-2]:
                tmp += [b[j]]
                counter_tmp += 1
                break

tmp.sort()
b = list(k for k,_ in itertools.groupby(tmp))


if args.assembler == "Stringtie2":
    for j in b:
        coverage = []
        for row in assemblyGTF.loc[j].iterrows():   
            coverage.append(row[1][8].split('cov "')[1].split('"')[0])   
        coverage = [ float(x) for x in coverage ]
        j.pop(coverage.index(max(coverage)))      
    toDelete = list(itertools.chain.from_iterable(b))


if args.assembler == "RATTLE":
    fastq = open(RATTLE_FASTQ,'r') ; fastq = fastq.readlines() ; fastqHeaders = [] 
    for line in fastq:
        if line.startswith('@'):
            fastqHeaders += [line]

    for j in b:
        clusters_names = []
        for dupCluster in assemblyGTF.loc[j].iterrows():   
            clusters_names.append(dupCluster[1][8].split('transcript_id "')[1].split('"')[0])  
        ass_reads = []
        for name in clusters_names:
            for read in fastq:
                if read.startswith('@'+name+' '):
                    ass_reads.append(int(read.split('=')[-1].split('\n')[0]))  
        j.pop(ass_reads.index(max(ass_reads)))
    toDelete = list(itertools.chain.from_iterable(b))


assemblyGTF_filtered = assemblyGTF.drop(assemblyGTF.index[toDelete])

################################################################################################################################################
#########    Remove CDS sequence to transcript assembly output    #############################################
################################################################################################################################################

CDSannotationsGTF = open(GTF2,'r') ; CDSannotationsGTF = CDSannotationsGTF.readlines() 

################################################################################################################################################
## Join multiple exons to one annotation containing the most extreme coordinates of the CDS
allGenes = []
for annotation in CDSannotationsGTF:
    allGenes += [annotation.split('"')[-2]]

multipleExonsGenes = []
dup = [k for k,v in Counter(allGenes).items() if v>1]
for i in dup:
    multipleExonsGenes.append([l for l,value in enumerate(allGenes) if value == i])

toDelete = []
for dupGroup in range(len(multipleExonsGenes)):
    coordinatesdupGroup = []
    for position in range(len(multipleExonsGenes[dupGroup])):
        coordinatesdupGroup.extend([int(CDSannotationsGTF[multipleExonsGenes[dupGroup][position]].split('\t')[3]),int(CDSannotationsGTF[multipleExonsGenes[dupGroup][position]].split('\t')[4])])
    newAnnotation = '\t'.join(CDSannotationsGTF[multipleExonsGenes[dupGroup][position]].split('\t')[:3])+'\t'+str(min(coordinatesdupGroup))+'\t'+str(max(coordinatesdupGroup))+'\t'+'\t'.join(CDSannotationsGTF[multipleExonsGenes[dupGroup][position]].split('\t')[5:])
    CDSannotationsGTF[multipleExonsGenes[dupGroup][0]] = newAnnotation       
    toDelete += [multipleExonsGenes[dupGroup][1:]]

toDelete = list(itertools.chain.from_iterable(toDelete))
for toDelete in sorted(toDelete, reverse=True):
    del CDSannotationsGTF[toDelete]
################################################################################################################################################

# Create dictionary of generated transcripts and candidate UTRs positions:
#   Keys -> Assembled transcripts IDs
#   Values -> dictionary :
#               Keys -> position of gene's candidate UTR, in CDS annotation GTF file of annotated transcript overlapping Stringtie2 transcript AND type of candidate UTR (5' or 3'). Example : 234.5UTR
#               Values -> coordinates of candidate UTR

if args.assembler == "Stringtie2":
    overlaps = {}
    for row in assemblyGTF_filtered.iterrows():
        overlaps[row[1][8].split('"')[3]] = {}

if args.assembler == "RATTLE":
    overlaps = {}
    for row in assemblyGTF_filtered.iterrows():
        overlaps[row[1][8].split('transcript_id "')[1].split('"')[0]] = {}


GTF_5UTR = [] ; GTF_3UTR = []

counter = -1
for line in CDSannotationsGTF:
    counter           += 1
    chromosome         = line.split('\t')[0]
    coordCDS           = set(range(int(line.split('\t')[3]),int(line.split('\t')[4]) + 1))
    orientation        = line.split('\t')[6]
    chromAssemblyCoord = []
    for row in assemblyGTF_filtered.loc[assemblyGTF_filtered[0] == chromosome].iterrows(): 
        chromAssemblyCoord.append(set(range(int(row[1][3]),int(row[1][4]) + 1)))
    for i in range(len(chromAssemblyCoord)):
        if list(chromAssemblyCoord[i] & coordCDS) == []:
            continue
        if list(chromAssemblyCoord[i] & coordCDS) != []: # if intersection in non-empty, then write it down:
            utrs = [list(group) for group in mit.consecutive_groups(sorted(chromAssemblyCoord[i].difference(coordCDS)))]
            if len(utrs) == 0:
                continue
            if args.assembler == "Stringtie2":
                transcript_id = assemblyGTF_filtered.loc[assemblyGTF_filtered[0] == chromosome].iloc[i][8].split('"')[3]
            if args.assembler == "RATTLE":
                transcript_id = assemblyGTF_filtered.loc[assemblyGTF_filtered[0] == chromosome].iloc[i][8].split('transcript_id "')[1].split('"')[0]
            if len(utrs) == 2:
                if orientation == '+':
                    GTF_5UTR.append(line.replace('CDS','5UTR').replace(str(min(coordCDS)),str(min(utrs[0]))).replace(str(max(coordCDS)),str(max(utrs[0]))))
                    overlaps[transcript_id][str(counter)+'.5UTR'] = [min(utrs[0]),max(utrs[0]),orientation]
                    GTF_3UTR.append(line.replace('CDS','3UTR').replace(str(min(coordCDS)),str(min(utrs[1]))).replace(str(max(coordCDS)),str(max(utrs[1]))))
                    overlaps[transcript_id][str(counter)+'.3UTR'] = [min(utrs[1]),max(utrs[1]),orientation]
                else:
                    GTF_3UTR.append(line.replace('CDS','3UTR').replace(str(min(coordCDS)),str(min(utrs[0]))).replace(str(max(coordCDS)),str(max(utrs[0]))))
                    overlaps[transcript_id][str(counter)+'.3UTR'] = [min(utrs[0]),max(utrs[0]),orientation]
                    GTF_5UTR.append(line.replace('CDS','5UTR').replace(str(min(coordCDS)),str(min(utrs[1]))).replace(str(max(coordCDS)),str(max(utrs[1]))))
                    overlaps[transcript_id][str(counter)+'.5UTR'] = [min(utrs[1]),max(utrs[1]),orientation]
            elif len(utrs) == 1:
                if max(utrs[0])+1 == min(coordCDS): 
                    if orientation == '+':
                        GTF_5UTR.append(line.replace('CDS','5UTR').replace(str(min(coordCDS)),str(min(utrs[0]))).replace(str(max(coordCDS)),str(max(utrs[0]))))
                        overlaps[transcript_id][str(counter)+'.5UTR'] = [min(utrs[0]),max(utrs[0]),orientation]  
                    else:
                        GTF_3UTR.append(line.replace('CDS','3UTR').replace(str(min(coordCDS)),str(min(utrs[0]))).replace(str(max(coordCDS)),str(max(utrs[0]))))
                        overlaps[transcript_id][str(counter)+'.3UTR'] = [min(utrs[0]),max(utrs[0]),orientation]  
                elif min(utrs[0])-1 == max(coordCDS):
                    if orientation == '+':
                        GTF_3UTR.append(line.replace('CDS','3UTR').replace(str(min(coordCDS)),str(min(utrs[0]))).replace(str(max(coordCDS)),str(max(utrs[0]))))
                        overlaps[transcript_id][str(counter)+'.3UTR'] = [min(utrs[0]),max(utrs[0]),orientation]  
                    else:
                        GTF_5UTR.append(line.replace('CDS','5UTR').replace(str(min(coordCDS)),str(min(utrs[0]))).replace(str(max(coordCDS)),str(max(utrs[0]))))
                        overlaps[transcript_id][str(counter)+'.5UTR'] = [min(utrs[0]),max(utrs[0]),orientation]  
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

# Create 2 lists splitting {overlaps} keys and values. 
counter     = -1
coordUTRs   = []
indexCDSUTR = [[] for i in range(len(overlaps))]
for i in overlaps.keys():
    counter += 1
    coordUTRs.append(list(overlaps[i].values()))
    for j in overlaps[i].keys(): 
        indexCDSUTR[counter].append(j)


# Split strands
splitStrands_indexCDSUTR = []
splitStrands_coordUTRs   = []
splitStrands_counters    = []
counter = -1
for i in coordUTRs:
    counter += 1
    if i == [] or len(i) == 1:
        continue
    else:
        coordsRangesList = []
        for l in list(set([int(q.split('.')[0]) for q in indexCDSUTR[counter]])):
            coordsRange = set(range(int(CDSannotationsGTF[l].split('\t')[3]),int(CDSannotationsGTF[l].split('\t')[4]) + 1)) ; coordsRangesList += [coordsRange]
        if list(set.intersection(*coordsRangesList)) != []: # if some CDSs overlap (intersection non-empty), split by filament
            senses = []
            for coord in i:
                senses += coord[2]
            if len(set(senses)) == 1: # if there aren't CDSs on both strands
                continue
            else:
                indices = [k for k, x in enumerate(senses) if x == '-']
                minus_toKeep_indexCDSUTR = []
                minus_toKeep_coordUTRs   = []
                minus_indeces            = []
                for idx in indices:
                    minus_toKeep_indexCDSUTR += [indexCDSUTR[counter][idx]]
                    minus_toKeep_coordUTRs   += [i[idx]]
                    splitStrands_counters    += [counter]
                    minus_indeces += [idx]
                splitStrands_indexCDSUTR += [minus_toKeep_indexCDSUTR]
                splitStrands_coordUTRs   += [minus_toKeep_coordUTRs]
                for j in sorted(minus_indeces, reverse=True):
                    coordUTRs[counter].remove(i[j]) ; indexCDSUTR[counter].remove(indexCDSUTR[counter][j]) 
        else:
            continue


splitStrands_counters = sorted(list(set(splitStrands_counters)),reverse = False)
revrange = reversed(range(len(splitStrands_counters)))
for n in revrange:
    coordUTRs.insert(splitStrands_counters[n],splitStrands_coordUTRs[n]) ; indexCDSUTR.insert(splitStrands_counters[n],splitStrands_indexCDSUTR[n])




# keep candidate UTRs at the extremities of long assembled Stringtie transcripts by removing them in the 2 lists previously created
for string in range(len(coordUTRs)):
    if len(coordUTRs[string]) == 1:
        coordUTRs[string].remove(coordUTRs[string][0])
        indexCDSUTR[string].remove(indexCDSUTR[string][0])
    elif len(coordUTRs[string]) > 1:
        CDSindices = [] # list of indices of CDS (in the CDSannotationsGTF list) overlapping an assembled string
        for candidateUTR in indexCDSUTR[string]:
            CDSindices.append(int(candidateUTR.split('.')[0]))
        CDSindices = list(OrderedDict.fromkeys(CDSindices)) 
        CDSranges  = []
        for CDSidx in CDSindices:
            CDSranges.append(set(range(int(CDSannotationsGTF[CDSidx].split('\t')[3]),int(CDSannotationsGTF[CDSidx].split('\t')[4]))))
        coord1CDS, max_index = max((x, (i,j)) for i,row in enumerate(CDSranges) for j,x in enumerate(row)) # coord1CDS -> second coordinate of the corresponding overlapping CDS
        coord0CDS, min_index = min((x, (i,j)) for i,row in enumerate(CDSranges) for j,x in enumerate(row)) # coord0CDS -> first coordinate of the corresponding overlapping CDS
        extremitiesCDS       = [min_index[0],max_index[0]]
        electedUTR           = []
        for candidateUTR in range(len(coordUTRs[string])):
            UTRtype     = indexCDSUTR[string][candidateUTR].split('.')[1] 
            orientation = CDSannotationsGTF[int(indexCDSUTR[string][candidateUTR].split('.')[0])].split('\t')[6]
            if int(indexCDSUTR[string][candidateUTR].split('.')[0]) == CDSindices[min_index[0]] and coordUTRs[string][candidateUTR][1] == coord0CDS - 1:
                if orientation == '+' and UTRtype == '5UTR':
                    electedUTR.append(candidateUTR)
                elif orientation == '-' and UTRtype == '3UTR':
                    electedUTR.append(candidateUTR)
            elif int(indexCDSUTR[string][candidateUTR].split('.')[0]) == CDSindices[max_index[0]] and coordUTRs[string][candidateUTR][0] == coord1CDS + 2:
                if orientation == '+' and UTRtype == '3UTR':
                    electedUTR.append(candidateUTR)
                elif orientation == '-' and UTRtype == '5UTR':
                    electedUTR.append(candidateUTR)
        for electedUTR in sorted(electedUTR, reverse=True):
            del coordUTRs[string][electedUTR] ; del indexCDSUTR[string][electedUTR]




# Delete unelected UTRs

toDelete_3UTR = [] ; toDelete_5UTR = []

pos1 = -1
for string in coordUTRs:
    pos1 += 1
    if string != []:
        pos2 = -1
        for trashedUTR in string:
            pos2      += 1
            gene_name  = CDSannotationsGTF[int(indexCDSUTR[pos1][pos2].split('.')[0])].split('"')[-2]
            UTRtype    = indexCDSUTR[pos1][pos2].split('.')[1]
            if UTRtype == '5UTR': 
                indices = [h for h, w in enumerate(GTF_5UTR) if gene_name in w] 
                for idx in indices:
                    if str(trashedUTR[0])+'\t'+ str(trashedUTR[1]) in GTF_5UTR[idx]: 
                        toDelete_5UTR.append(idx) 
            elif UTRtype == '3UTR': 
                indices = [h for h, w in enumerate(GTF_3UTR) if gene_name in w] 
                for idx in indices:
                    if str(trashedUTR[0])+'\t'+ str(trashedUTR[1]) in GTF_3UTR[idx]: 
                        toDelete_3UTR.append(idx) 

for toDelete_3UTR in sorted(set(toDelete_3UTR),reverse = True):
    del GTF_3UTR[toDelete_3UTR]

for toDelete_5UTR in sorted(set(toDelete_5UTR),reverse = True):
    del GTF_5UTR[toDelete_5UTR] 

GTF_5UTR = list(OrderedDict.fromkeys(GTF_5UTR)) ; GTF_3UTR = list(OrderedDict.fromkeys(GTF_3UTR)) 


## remove small candidates UTRs 
genes = []
for line in GTF_5UTR:
    gene_name = line.split('"')[-2] ; genes += [gene_name]

duplicates = [item for item,count in collections.Counter(genes).items() if count > 1]

toFilter = []
for i in duplicates:
    indices = [h for h, w in enumerate(GTF_5UTR) if i in w]
    lengthList = []
    for idx in indices:
        length = int(GTF_5UTR[idx].split('\t')[4]) - int(GTF_5UTR[idx].split('\t')[3]) + 1 ; lengthList += [length]
    indices.remove(indices[lengthList.index(max(lengthList))]) # keep the longest sequence
    for garbage in indices:
        toFilter.append(garbage)

for toFilter in sorted(toFilter, reverse=True):
    del GTF_5UTR[toFilter]


genes = []
for line in GTF_3UTR:
    gene_name = line.split('"')[-2] ; genes += [gene_name]

duplicates = [item for item,count in collections.Counter(genes).items() if count > 1]

toFilter = []
for i in duplicates:
    indices = [h for h, w in enumerate(GTF_3UTR) if i in w]
    lengthList = []
    for idx in indices:
        length = int(GTF_3UTR[idx].split('\t')[4]) - int(GTF_3UTR[idx].split('\t')[3]) + 1 ; lengthList += [length]
    indices.remove(indices[lengthList.index(max(lengthList))])
    for garbage in indices:
        toFilter.append(garbage)

for toFilter in sorted(toFilter, reverse=True):
    del GTF_3UTR[toFilter]

GTF_3UTR_filtered = []
for line in GTF_3UTR:
    if '"SPNCRNA' not in line and '"SPATRNA' not in line and '"SPBTRNA' not in line and '"SPCTRNA' not in line and '"SPRRNA' not in line and '"SPSNORNA' not in line and '"SPSNRNA' not in line and '"SPMITTRNA' not in line: # don't take into account non-coding RNA
        GTF_3UTR_filtered += [line]

GTF_5UTR_filtered = []
for line in GTF_5UTR:
    if '"SPNCRNA' not in line and '"SPATRNA' not in line and '"SPBTRNA' not in line and '"SPCTRNA' not in line and '"SPRRNA' not in line and '"SPSNORNA' not in line and '"SPSNRNA' not in line and '"SPMITTRNA' not in line: 
        GTF_5UTR_filtered += [line]

if args.assembler == "Stringtie2":
    final_GTF_5UTR = open(outputDir+'Stringtie2_5UTR.gtf', 'w') ; final_GTF_5UTR.writelines(GTF_5UTR_filtered) ; final_GTF_5UTR.close()
    final_GTF_3UTR = open(outputDir+'Stringtie2_3UTR.gtf', 'w') ; final_GTF_3UTR.writelines(GTF_3UTR_filtered) ; final_GTF_3UTR.close()

if args.assembler == "RATTLE":
    final_GTF_5UTR = open(outputDir+'RATTLE_5UTR.gtf', 'w') ; final_GTF_5UTR.writelines(GTF_5UTR_filtered) ; final_GTF_5UTR.close()
    final_GTF_3UTR = open(outputDir+'RATTLE_3UTR.gtf', 'w') ; final_GTF_3UTR.writelines(GTF_3UTR_filtered) ; final_GTF_3UTR.close()

exit(0)
