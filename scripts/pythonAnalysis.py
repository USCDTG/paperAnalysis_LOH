##### IMPORTS #####
import os
import re
import io
import sys
import glob
import pysam
import pandas
import datetime
import argparse
import subprocess
import simplesam
from subprocess import Popen
from pysam import VariantFile
from posixpath import expanduser
from numpy.core.numeric import NaN

def runAC(bam,sample,chr,start,stop,df,output,finalOutput):
    # print("- Count Coverage")
    # print("- Reading BAM")
    samfile = pysam.AlignmentFile(bam,"rb")
    df = df.reset_index(drop=True)
    df_2 = pandas.DataFrame([samfile.count_coverage(contig = str(a), start = b, stop = c)] for a,b,c in zip(chr,start,stop))
    nucleotides = pandas.DataFrame([x[0][0],x[1][0],x[2][0],x[3][0]] for x in df_2[0].tolist())
    nucleotides.columns = ['A','C','G','T']
    # print("- Merging dataframes")
    final = pandas.merge(df,nucleotides,how='outer',left_index=True, right_index=True)
    final = final[final['REF'].str.len() <= 1]
    final = final[final['ALT'].str.len() <= 1]
    final.drop(final.index[final['ALT'] == '*'], inplace=True)
    # print("- Looping through values")
    cols_of_interest = ['A','T','C','G']
    save = final.loc[final[final[cols_of_interest]!=0].dropna(thresh=1).index]
    final=save.reset_index()
    index = 0
    refCOUNT = []
    for value in final['REF']:
        refCOUNT.append(final[value][index])
        index += 1
    index = 0
    altCOUNT = []
    for value2 in final['ALT']:
        altCOUNT.append(final[value2][index])
        index += 1
    final['REF_COUNT_CALC'] = refCOUNT
    final['ALT_COUNT_CALC'] = altCOUNT
    final['TOTAL'] = final['REF_COUNT_CALC'] + final['ALT_COUNT_CALC']
    final['cluster'] = sample
    final = final[final.TOTAL > 1]
    final.drop(columns='index',inplace=True)
    final.to_csv('%s/%s_alleleCounts.csv' % (finalOutput, sample), index=False)
    samfile.close()

inputVCF = sys.argv[1]
bam = sys.argv[2]
clusterName = sys.argv[3]
output = sys.argv[4]

with open(inputVCF,'r') as f:
    lines = [x for x in f if not x.startswith('##')]
f.close()
importedVCF = pandas.read_csv(io.StringIO(''.join(lines)), sep='\t', low_memory=False)
importedVCF['POSITION_oneBefore'] = importedVCF['POS'] - 1
importedVCF['ANN'] = importedVCF['INFO']
chr = importedVCF['#CHROM'].tolist()
start = importedVCF['POSITION_oneBefore'].tolist()
stop = importedVCF['POS'].tolist()

runAC(bam,clusterName,chr,start,stop,importedVCF,output,output)
