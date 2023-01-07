## nohup python3 -u similarty_proteinpairsimilarityCUT.py > similarty_proteinpairsimilarityCUT1214out.txt 2>&1 &

### run.sh
#!/bin/bash
#
#for (( i==0; i<10; i++ ));
#do
#CUTstep=1
#CUTNUM=$((CUTstep * i))
#outputfile="similarty_proteinpairsimilarityCUTout_$CUTNUM.txt"
#cmd="python3 -u similarty_proteinpairsimilarityCUT.py $CUTNUM > $outputfile & "
#echo "$cmd"
#sleep 1
#python3 -u similarty_proteinpairsimilarityCUT.py $CUTNUM > $outputfile 2>&1 &
#done
#

import sys
cutnum=int(sys.argv[1])



import re
from re import search
import os, sys
from os import walk
import traceback
import subprocess
import pandas as pd
import time
import math
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import pearsonr

def minmaxavg(numlist):
    remove_list=[0,100]
    #numlist=[float(i) for i in numlist if i not in remove_list]
    
    numlist=[i for i in numlist if i <100]
    numlist=[i for i in numlist if i >0]
    numlist.sort()
    try:
        lstavg=sum(numlist)/len(numlist)
        lstavgT5=sum(numlist[-5:])/5
        #lstavgL5=sum(numlist1[:5])/5
    except:
        lstavg=np.nan
        lstavgT5=np.nan
        #lstavgL5=np.nan
    try:
        lstmin=min(numlist)
        lstmax=max(numlist)
    except:
        lstmin=np.nan
        lstmax=np.nan
        
    return(lstmin,lstmax,lstavg,lstavgT5)

def seq_similarity (seq1_list,seq2_list):
    #aminoacid_table1={'X':'0','-':'0','C':'1','G':'1','P':'1','A':'2','V':'2','I':'2','L':'2','M':'2','F':'2','Y':'2','W':'2','S':'3','T':'3','N':'3','Q':'3','D':'4','E':'4','R':'5','H':'5','K':'5'}
    #aminoacid_table2={'X':'00000','-':'00000','C':'00001','G':'00001','P':'00001','A':'00010','V':'00010','I':'00010','L':'00010','M':'00010','F':'00010','Y':'00010','W':'00010','S':'00100','T':'00100','N':'00100','Q':'00100','D':'01000','E':'01000','R':'10000','H':'10000','K':'10000'}
    #aminoacid_table3={'X':'a','-':'a','C':'b','G':'b','P':'b','A':'c','V':'c','I':'c','L':'c','M':'c','F':'c','Y':'c','W':'c','S':'d','T':'d','N':'d','Q':'d','D':'e','E':'e','R':'f','H':'f','K':'f'}
    s=0
    t=0
    for i,aa1 in enumerate(seq1_list):
        aa1=seq1_list[i]
        aa2=seq2_list[i]
        if aa1+aa2 != 0:
            t=t+1
            if aa1==aa2:
                s=s+1
    try:
        ss=s/t*100
    except:
        ss=0
    return(ss)

def pair_seqCUT_similarty(workpath,ppalignfile,cutnum,seqsimCCnames):
    timestamp=time.strftime("%m%d",time.localtime(time.time()))
    alignmentDF=pd.DataFrame()
    os.chdir(workpath)
    protalign=pd.read_csv(ppalignfile,header=None ) #,names=protalign_names)
    protalign1=protalign.T.set_index(1).T
    prot_conv=protalign1.to_dict('list')
    
    protalign_output='similarty_proteinpairsimilarity'+timestamp+'_'+re.sub('.csv','',ppalignfile)+'_'+str(cutnum)+'CUT'+'.csv'

    alnum=0
    for i,ppp1 in enumerate (protalign1.columns):
        for j,ppp2 in enumerate (protalign1.columns):
            seq1=prot_conv[ppp1][1]
            seq2=prot_conv[ppp2][1]
            
            try:
                seq1_list = [int(re.sub('\s|\[|\]|\'','',x)) for x in seq1.split(',')]
                seq2_list = [int(re.sub('\s|\[|\]|\'','',x)) for x in seq2.split(',')]
            except:
                seq1_list=[]
                seq2_list=[]
                continue

            alignmentDF.at[alnum,'seq1']=ppp1
            alignmentDF.at[alnum,'seq2']=ppp2

            sscut=[]
            sscut.append(seq_similarity(seq1_list,seq2_list))
            
            try:
                seqcutSTEP=int(len(seq1_list)/(cutnum))
                
                passtail=0
                for startlist in range (0,int(len(seq1_list)),seqcutSTEP):
                    if passtail==1:
                        continue
                    stoplist_next=startlist+seqcutSTEP+seqcutSTEP
                    if stoplist_next < len(seq1_list): 
                        stoplist=startlist+seqcutSTEP
                        seq1_list_cut=seq1_list[startlist:stoplist]
                        seq2_list_cut=seq2_list[startlist:stoplist]
                    else:
                        stoplist=len(seq1_list)
                        seq1_list_cut=seq1_list[startlist:]
                        seq2_list_cut=seq2_list[startlist:]
                        passtail=1

                    sscut.append(seq_similarity(seq1_list_cut,seq2_list_cut))
            except:
                seqcutSTEP=int(len(seq1_list))
                print (' *** ERROR protein sequence cannot be cutted ***')
 
            for i, seqsim_name in enumerate(seqsimCCnames):
                try:
                    alignmentDF.at[alnum,seqsim_name]="{:.2f}".format(sscut[i])
                except:
                    pass

            alnum=alnum+1

        os.chdir(workpath)
        alignmentDF.to_csv(protalign_output,index=False,encoding='utf-8')

    os.chdir(workpath)
    alignmentDF.to_csv(protalign_output,index=False,encoding='utf-8')
    return(protalign_output)

def proteinname_seqCUT_similarty(workpath,protalign_output,cutnum,seqsimCCnames):
    os.chdir(workpath)
    alignmentDF=pd.read_csv(protalign_output)
    alignmentDF_group=alignmentDF.groupby(['seq1'])
    similartydf=pd.DataFrame()
    seqsimCCnames=[]
    for cc in alignmentDF.columns:
        if re.search('seq_similarity',cc):
            seqsimCCnames.append(cc)
    #cutnum=len(seqsimCCnames)
    proteinsimilarity_output=re.sub('.csv','',protalign_output)+'_SIMminmaxavg'+'.csv'
    ppnum=0
    for i,ppp in enumerate(alignmentDF_group.first().index):
        ppp_groupDF=alignmentDF_group.get_group(ppp)
        similartydf.at[ppnum,'protein']=ppp
        lstmin_list=[]
        lstmax_list=[]
        lstavg_list=[]
        lstavgT5_list=[]
        for cut,seq_simnum in enumerate(seqsimCCnames):
            similarityslst=ppp_groupDF[seq_simnum].to_list()
            lstmin,lstmax,lstavg,lstavgT5=minmaxavg(numlist=similarityslst)
            lstmin_name='SIMmin_'+str(cutnum)+'CUT'+str(cut)
            lstmax_name='SIMmax_'+str(cutnum)+'CUT'+str(cut)
            lstavg_name='SIMavg_'+str(cutnum)+'CUT'+str(cut)
            lstavgT5_name='SIMavgT5_'+str(cutnum)+'CUT'+str(cut)
            try:
                similartydf.at[ppnum,lstmin_name]=lstmin
                similartydf.at[ppnum,lstmax_name]=lstmax
                similartydf.at[ppnum,lstavg_name]=lstavg
                similartydf.at[ppnum,lstavgT5_name]=lstavgT5
            except:
                print ('** no result: ',ppnum,ppp,lstmin,lstmax,lstavg,lstavgT5)
                pass
                
        ppnum=ppnum+1
    #simcolumnlist1=similartydf.columns.to_list()
    #simcolumnlist1.remove('protein')
    os.chdir(workpath)
    similartydf.to_csv(proteinsimilarity_output,index=False,encoding='utf-8')
    return(proteinsimilarity_output)

######################################## Main ########################################



timestamp=time.strftime("%Y%m%d%H%M",time.localtime(time.time()))
starttime=time.time()
workpath='/home1/drhuangwc/AAALAB/GPCR5b/Interspecies_20220915_AgonistAntagonist_datacode2ALLcolumns'
os.chdir(workpath)

#protalign_names=pd.read_csv('proteinseq20220915_align_withMOR1var.csv').loc[0].tolist()

ppalignfile='proteinseq20220915_align_withMOR1var.csv'

#for cutnum in [1,2,3,5,7,8,10]:

i=0
print ('cut=',cutnum)
seqsimCCnames=[]
if cutnum<2:
    cutnum=1
    cutnum_range=1
else:
    cutnum_range=cutnum+1
for i in range (cutnum_range):
    listname='seq_similarityCUT'+str(cutnum)+'_'+str(i)
    seqsimCCnames.append(listname)
print (seqsimCCnames)

protalign_output=pair_seqCUT_similarty(workpath,ppalignfile,cutnum,seqsimCCnames)
proteinsimilarity_output=proteinname_seqCUT_similarty(workpath,protalign_output,cutnum,seqsimCCnames)


print ('done',cutnum,'\n output to : \n',protalign_output,' \n',proteinsimilarity_output)


if (time.time()-starttime)/3600 > 1 :
    print(time.strftime("%Y/%m/%d %H:%M:%S",time.localtime(time.time())), '\n Total = ', (time.time()-starttime)/3600, ' hrs \n' )
else:
    print(time.strftime("%Y/%m/%d %H:%M:%S",time.localtime(time.time())), '\n Total = ', (time.time()-starttime)/60, ' mins \n' )

   