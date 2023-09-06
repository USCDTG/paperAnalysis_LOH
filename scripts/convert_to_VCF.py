import pandas
import glob
import sys

def convertCSVtoVCF(sample,pathToUse):
    print('- Run as input directory of per-cluster BAMs')
    path = pathToUse
    listOfAlleleCounts = glob.glob(path + "/*.csv")
    print("Arranging all AlleleCount files")
    pathImport = []
    for filename in listOfAlleleCounts:
        df = pandas.read_csv(filename, index_col = None, low_memory=False)
        pathImport.append(df)
    ac = pandas.concat(pathImport, axis = 0, ignore_index=True)
    ac = ac.rename(columns={"cluster": "CLUSTER", "ID": "rsID"})
    ac['CLUSTER'] = ac['CLUSTER'].astype(str)
    ac['CLUSTER'] = ac['CLUSTER'].str.replace('cluster', '')
    ac['CLUSTER'] = ac['CLUSTER'].str.split('_').str[1]
    if ac['#CHROM'].str.contains('chr').any():
        ac['#CHROM'] = ac['#CHROM'].str.replace('chr', '')
    ac = ac.drop(['ANN'], axis = 1)
    ac['TOTAL'] = ac['ALT_COUNT_CALC'] + ac['REF_COUNT_CALC']
    ac['AF'] = ac['ALT_COUNT_CALC'] / ac['TOTAL']
    ac = ac.drop(['POSITION_oneBefore'], axis = 1)
    ac['QUAL'] = '.'
    ac['GT'] = 'germlineGT:0/1'
    ac[['#CHROM','POS','rsID','QUAL','REF','ALT','GT','REF_COUNT_CALC','ALT_COUNT_CALC','CLUSTER','TOTAL','AF']]
    ac = ac[ac.rsID != '.']
    print("Preparing for Pivot")
    ac = ac.rename(columns={"#CHROM": "CHR"})
    ac['FORMAT'] = 'GT:AD:DP'
    ac = ac.drop(['AF'], axis=1)
    ac['AD'] = ac['REF_COUNT_CALC'].astype(str) + ',' + ac['ALT_COUNT_CALC'].astype(str)
    ac['SAMPLE'] = './.:' + ac['AD'].astype(str) + ":" + ac['TOTAL'].astype(str)

    print("Pivot")
    toMerge = ac.pivot(columns='CLUSTER', values = 'SAMPLE', index = 'rsID')
    toMerge['ID'] = toMerge.index
    formatVCF = pandas.DataFrame(data = {'#CHROM': ac['CHR'], 'POS': ac['POS'], 'ID': ac['rsID'], 'REF': ac['REF'], 'ALT': ac['ALT'], 'QUAL': ac['QUAL'], 'FILTER':'.', 'INFO': ac['GT'], 'FORMAT': 'GT:AD:DP'}, columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
    testVCF = pandas.merge(formatVCF, toMerge, on='ID')
    testVCF = testVCF.drop_duplicates()
    testVCF = testVCF[testVCF['#CHROM'] != 'X']
    testVCF['#CHROM']  = testVCF['#CHROM'].astype(int)
    testVCF['POS']  = testVCF['POS'].astype(int)
    testVCF = testVCF.sort_values(['#CHROM','POS'])
    print("VCF Construction")
    header = """##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=germlineGT,Number=1,Type=String,Description="Genotype from germline exome variant calling">
##contig=<ID=1,length=248956422>
##contig=<ID=2,length=242193529>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##contig=<ID=chrM,length=16569> \n"""

    output_VCF = sample + "_tLOH.vcf"
    with open(output_VCF, 'w') as vcf:
        vcf.write(header)
    testVCF.to_csv(output_VCF, sep = "\t", mode = 'a', index = False)



sample = sys.argv[1]
path = sys.argv[2]
convertCSVtoVCF(sample,path)
