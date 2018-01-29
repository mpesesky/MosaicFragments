import pandas as pd
from Bio import SeqIO
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def get_count_table(df):
    countDF = pd.DataFrame(df.groupby('qseqid').size(), columns=['Occurrences'])
    countDF['qseqid'] = countDF.index
    return countDF

def import_data(datafile):
    return pd.read_table(datafile, sep="\t")

def import_gb(gbFile):
    records = list(SeqIO.parse(gbFile, 'genbank'))
    recordDict = {x.id: x for x in records}
    return recordDict

def identify_targets(dataframe, targetList):
    subframe = dataframe[['qseqid', 'sseqid', 'sstart', 'send']][dataframe['qseqid'] .isin(targetList)]
    return subframe

def convert_loci(df):
    df['sloci'] = df['sseqid'].map(lambda x: x.split("|")[1])
    return df

def surrounding_genes(records, targetFrame, distance, up):
    geneDict = {}
    for i, row in targetFrame.iterrows():
        targetRecord = records[row['sloci']]

        fragment = targetRecord[row['sstart']-distance:row['send'] + distance]
        if row['qseqid'] not in geneDict.keys():
            geneDict[row['qseqid']] = {}
        for feat in fragment.features:
            if feat.type == 'CDS':
                try:
                    annot = feat.qualifiers['product'][0]
                except KeyError:
                    try:
                        annot = feat.qualifiers['note'][0]
                    except KeyError:
                        annot = feat.type
                if annot in geneDict[row['qseqid']].keys():
                    geneDict[row['qseqid']][annot] += 1
                else:
                    geneDict[row['qseqid']][annot] = 1
    return geneDict

def top10(dictOfDicts, outfile):
    df = pd.DataFrame(dictOfDicts)
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    colNum = len(df.columns)
    df.rename(columns={df.columns[i]:alphabet[i] for i in range(0, colNum)}, inplace=True)

    for col in df.columns:
        df[col].nlargest(50).to_csv(outfile, mode='a', sep="\t", header=True)


#    topDF = df.nlargest(10, 'A')
#    colors = ['#8c510a', '#bf812d', '#dfc27d', '#f6e8c3', '#f5f5f5', '#c7eae5', '#80cdc1', '#35978f',
#              '#01665e', '#000000']
#    plt.ioff()
#    print(topDF)
#    topDF.transpose().plot.bar(stacked=True, color=colors)
#    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
#    ax.bar(topDF.iloc[0], color=colors[0])
#    bottomSpot = topDF.iloc[0]
#    for i in range(1,9):
#        ax.bar(topDF.iloc[i], color=colors[i], bottom=bottomSpot)
#        bottomSpot += topDF.iloc[i]
#    plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
#    plt.tight_layout()
#    plt.savefig(outfile)



parser = argparse.ArgumentParser(description="Tabulate gene annotations for genes surrounding common fragments")

parser.add_argument("FragmentTable", help="TSV output by blast_filtering.py")
parser.add_argument("GenBank", help="Genbank file of plasmid set")
parser.add_argument("-t", "--threshold", default=900, type=int, help="Minimum abundance to be included")
parser.add_argument("-n", "--distance", default=5000, type=int, help="Distance (in bp) on either side of fragment to include")
parser.add_argument("-o", "--out", default=None, help="Created a stacked bar chart with the given name")
parser.add_argument("-u", "--upstream", action='store_true', help="Look at upstream, rather than downstream genes")

args = parser.parse_args()

df = import_data(args.FragmentTable)
gb = import_gb(args.GenBank)

countDF = get_count_table(df)
abundantList = list(countDF['qseqid'][countDF['Occurrences'] > args.threshold])
print(abundantList)
targetFrame = identify_targets(df, abundantList)
targetFrame = convert_loci(targetFrame)
genes = surrounding_genes(gb, targetFrame, args.distance, args.upstream)
restructured = {}
for gene in genes.keys():
    restructured[gene] = {}
    for annot in genes[gene].keys():
        if ('transposase' in annot) or ('TnpA' in annot):
#            pass
            if 'transposase' in restructured[gene].keys():
                restructured[gene]['transposase'] += genes[gene][annot]
            else:
                restructured[gene]['transposase'] = genes[gene][annot]
#        elif 'hypothetical protein' in annot:
#            pass
        else:
            restructured[gene][annot] = genes[gene][annot]

top10(restructured, args.out)