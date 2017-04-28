#!/usr/bin/python

import pandas as pd
import numpy as np
from tqdm import tqdm
import sys


file_stringtie_gff    = sys.argv[1] #'/ref/analysis/Cre/braker/braker.try5_mario/guided/intron3000.merge.sorted.bam.gffguided.stringtie.gff.strandcor.gff'
file_transdecoder_gff = sys.argv[2] #'/ref/analysis/Cre/braker/braker.try5_mario/guided/transcripts.fasta.transdecoder.gff3'
#file_transdecoder_gff = 'test.gff'

df_stringtie_gff                  = pd.read_csv(file_stringtie_gff,sep='\t',header=None,comment='#')
df_transdecoder_gff               = pd.read_csv(file_transdecoder_gff,sep='\t',header=None)
df_stringtie_gff['transcript_id'] = df_stringtie_gff[8].apply(lambda x : x.split(';')[1].replace('transcript_id','').strip().strip('"'))

df_transdecoder_gff_ix = df_transdecoder_gff.set_index(0)
contiglist             = set(df_transdecoder_gff[0].values)
df_stringtie_gff_ix    = df_stringtie_gff.set_index('transcript_id')


'''
[1.2.3.4.5.6.7.8.9.21.22.23.24]
-> [[1,9],[21,24]]
'''
def get_cont(listin): 
    result = []
    comp   = [listin[0]]
    for n,val in enumerate(listin):
        try:
            if listin[n+1] - val > 1:
                comp.append(val)
                result.append(comp)
                comp = [listin[n+1]]
        except IndexError:
            comp.append(val)
            result.append(comp)
    return result

'''
GFF3 format
seqname - The name of the sequence. Typically a chromosome or a contig. Argo does not care what you put here. It will superimpose gff features on any sequence you like.
source - The program that generated this feature. Argo displays the value of this field in the inspector but does not do anything special with it.
feature - The name of this type of feature. The official GFF3 spec states that this should be a term from the SOFA ontology, but Argo does not do anything with this value except display it.
start - The starting position of the feature in the sequence. The first base is numbered 1.
end - The ending position of the feature (inclusive).
score - A score between 0 and 1000. If there is no score value, enter ".".
strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'. Argo does not do anything with this field except display its value.
GFF3: grouping attributes Attribute keys and values are separated by '=' signs. Values must be URI encoded.quoted. Attribute pairs are separated by semicolons. Certain, special attributes are used for grouping and identification (See below). This field is the one important difference between GFF flavors.'''
dic_gff = {'Chromosome' : [],
           'source'     : [],
           'feature'    : [],
           'start'      : [],
           'end'        : [],
           'score'      : [],
           'strand'     : [],
           'frame'      : [],
           'GFF3: grouping attributes' : []}


for target_contig in tqdm(contiglist):
    #arget_contig =  'STRG.10.2'
    dic    = {}
    tpos   = 1 # transcript pos
    sub_df = df_stringtie_gff_ix.loc[target_contig].reset_index()
    mask   = (sub_df[2] == 'exon')
    for c,a,b,d in sub_df[mask].sort_values(by=3)[[0,3,4,6]].values:
        for gpos in np.arange(a,b+1):
            dic[tpos] = c,d,gpos
            tpos += 1

    # sub_df_tranddecoder parse
    sub_df                 = df_transdecoder_gff_ix.loc[target_contig].reset_index()
    sub_df                 = sub_df[(sub_df[6] == '+')] # only predictions that follow stringtie
    if len(sub_df) == 0:
        continue
    selected_index         = (sub_df[sub_df[2]=='CDS'][4] - sub_df[sub_df[2]=='CDS'][3]).sort_values(ascending=False).index[0] # Select longest CDS from Transdecoder result. 
    selected_ORF           = sub_df.loc[selected_index][8].split(';')[1]
    #sub_df['Parent']       = sub_df[8].apply(lambda x : x.split(';')[1]==selected_ORF)
    sub_df['Parent'] = sub_df[8].apply(lambda x : (x.split(';')[1]==selected_ORF) or (x.split(';')[0]==selected_ORF.replace('Parent=','ID=')))
    sub_df                 = sub_df[sub_df['Parent']==True]
    sub_df['predictionID'] = sub_df[8].apply(lambda x : x.split(';')[0]) 
    # end 
    if len(sub_df) == 0 :
        continue
    ''' 
    # gene part
    mask = (sub_df[2] == 'gene')
    left,right = sub_df[mask][[3,4]].values[0]
    chromosome = dic[left][0]
    strand     = dic[left][1]
    cleft,cright = dic[left][2] , dic[right][2]
    genename     = sub_df[mask][0].values[0]
    info         = 'ID=%s;Name=%s'%(genename,genename)
    dic_gff['Chromosome'].append(chromosome)
    dic_gff['source'].append('ST/TD/KYJ')
    dic_gff['feature'].append('gene')
    dic_gff['start'].append(cleft)
    dic_gff['end'].append(cright)
    dic_gff['score'].append('.')
    dic_gff['strand'].append(strand)
    dic_gff['frame'].append('.')
    dic_gff['GFF3: grouping attributes'].append(info)
    # gene part end
    '''

    # mRNA part
    mask = (sub_df[2] == 'mRNA')
    for left,right,ID in sub_df[mask][[3,4,'predictionID']].values:
        chromosome = dic[left][0]
        strand     = dic[left][1]
        cleft,cright = dic[left][2] , dic[right][2]
        mRNAname     = ID.replace('ID=','')
        genename     = '.'.join(mRNAname.split('|')[0].split('.')[0:2])
        info         = 'ID=%s;Parent=%s;Name=%s'%(mRNAname,genename,mRNAname)
        dic_gff['Chromosome'].append(chromosome)
        dic_gff['source'].append('ST/TD/KYJ')
        dic_gff['feature'].append('mRNA')
        dic_gff['start'].append(cleft)
        dic_gff['end'].append(cright)
        dic_gff['score'].append('.')
        dic_gff['strand'].append(strand)
        dic_gff['frame'].append('.')
        dic_gff['GFF3: grouping attributes'].append(info)
    # mRNA part end

    # cds part
    mask = (sub_df[2] == 'CDS')
    for left,right,ID in sub_df[mask][[3,4,'predictionID']].values:
        chromosome = dic[left][0]
        strand     = dic[left][1]

        transcript_len = len(dic.values())
        if strand == '+':
            covs = [dic[x][2] for x in np.arange(left,right+1)]
        else:
            covs = [dic[transcript_len-x+1][2] for x in np.arange(left,right+1)]
        covs.sort()
        cont = get_cont(covs)
        if strand == '-':
            cont.sort(key=lambda x : x[0],reverse=True)
        frame     = 0
        accum_len = 0
        for n,(cleft,cright) in enumerate(cont):
            CDSname     = ID.replace('ID=','')
            mRNAname     = CDSname.replace('cds.','')
            info         = 'ID=%s.%d;Parent=%s;Name=%s'%(CDSname,n,mRNAname,mRNAname)
            dic_gff['Chromosome'].append(chromosome)
            dic_gff['source'].append('ST/TD/KYJ')
            dic_gff['feature'].append('CDS')
            dic_gff['start'].append(cleft)
            dic_gff['end'].append(cright)
            dic_gff['score'].append('.')
            dic_gff['strand'].append(strand)
            dic_gff['frame'].append(frame)
            dic_gff['GFF3: grouping attributes'].append(info)
            accum_len += (cright - cleft + 1)
            frame = (3 - accum_len%3)%3
    # cds part end

    # five_prime_UTR part
    mask = (sub_df[2] == 'five_prime_UTR')
    for left,right,ID in sub_df[mask][[3,4,'predictionID']].values:
        chromosome = dic[left][0]
        strand     = dic[left][1]

        transcript_len = len(dic.values())
        if strand == '+':
            covs = [dic[x][2] for x in np.arange(left,right+1)]
        else:
            covs = [dic[transcript_len-x+1][2] for x in np.arange(left,right+1)]
        covs.sort()
        cont = get_cont(covs)
        if strand == '-':
            cont.sort(key=lambda x : x[0],reverse=True)

        for n,(cleft,cright) in enumerate(cont):
            FPUname     = ID.replace('ID=','')
            mRNAname     = FPUname.replace('.utr5p1','')
            info         = 'ID=%s.%d;Parent=%s;Name=%s'%(FPUname,n,mRNAname,mRNAname)
            dic_gff['Chromosome'].append(chromosome)
            dic_gff['source'].append('ST/TD/KYJ')
            dic_gff['feature'].append('five_prime_UTR')
            dic_gff['start'].append(cleft)
            dic_gff['end'].append(cright)
            dic_gff['score'].append('.')
            dic_gff['strand'].append(strand)
            dic_gff['frame'].append('.')
            dic_gff['GFF3: grouping attributes'].append(info)
    # five_prime_UTR part end

    # three_prime_UTR part
    mask = (sub_df[2] == 'three_prime_UTR')
    for left,right,ID in sub_df[mask][[3,4,'predictionID']].values:
        chromosome = dic[left][0]
        strand     = dic[left][1]
 
        transcript_len = len(dic.values())
        if strand == '+':
            covs = [dic[x][2] for x in np.arange(left,right+1)]
        else:
            covs = [dic[transcript_len-x+1][2] for x in np.arange(left,right+1)]
        covs.sort()
        cont = get_cont(covs)
        if strand == '-':
            cont.sort(key=lambda x : x[0],reverse=True)
 
        for n,(cleft,cright) in enumerate(cont):
            TPUname     = ID.replace('ID=','')
            mRNAname     = TPUname.replace('.utr3p1','')
            info         = 'ID=%s.%d;Parent=%s;Name=%s'%(TPUname,n,mRNAname,mRNAname)
            dic_gff['Chromosome'].append(chromosome)
            dic_gff['source'].append('ST/TD/KYJ')
            dic_gff['feature'].append('three_prime_UTR')
            dic_gff['start'].append(cleft)
            dic_gff['end'].append(cright)
            dic_gff['score'].append('.')
            dic_gff['strand'].append(strand)
            dic_gff['frame'].append('.')
            dic_gff['GFF3: grouping attributes'].append(info)
    # three_prime_UTR part end


df_new_gff = pd.DataFrame(dic_gff,columns=['Chromosome','source','feature','start','end','score','strand','frame','GFF3: grouping attributes'])
df_new_gff['Name'] = df_new_gff['GFF3: grouping attributes'].apply(lambda x : x.split(';')[-1])
df_new_gff.sort_values(by=['Name','start']).to_csv(file_transdecoder_gff+'.genome.gff',sep='\t')

