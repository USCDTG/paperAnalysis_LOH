import pysam
import pandas
import sys
import os
import glob
import argparse

parser = argparse.ArgumentParser(description = "SPLIT BAMS")
parser.add_argument('--bam',help='path to bam')
parser.add_argument('--clusters',help='path to clusters.csv')
parser.add_argument('--output',help='output dir name')
args = parser.parse_args()
bam = args.bam
clusters = args.clusters
output = args.output

def splitBAMsByCluster(clusters,bam,output):
    print(clusters)
    print(bam)
    print(output)
    df = pandas.read_csv(clusters)
    fh = pysam.Samfile(bam)
    cluster_map = dict(zip(df.Barcode, df.Cluster))
    os.makedirs(output)
    out_handles = {}
    for c in set(df.Cluster):
        o = pysam.Samfile(os.path.join(output, 'cluster' + str(c) + '.bam'), 'wb', template=fh)
        out_handles[c] = o
    for rec in fh:
        if rec.has_tag('CB') and rec.get_tag('CB') in cluster_map:
            cluster = cluster_map[rec.get_tag('CB')]
            out_handles[cluster].write(rec)
    for h in out_handles.values():
        h.close()
    clusterBamList = glob.glob(os.path.join(output,'cluster?.bam'))

splitBAMsByCluster(clusters,bam,output)
