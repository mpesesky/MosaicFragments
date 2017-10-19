import pandas as pd
import Fasta_one_line as fol
import argparse
import matplotlib.pyplot as plt
import numpy as np

def col_add(transposase, Abres):
    if (transposase == "Yes") and (Abres == "Yes"):
        return "Both"
    elif (transposase == "Yes") and (Abres == "No"):
        return "Transposase"
    elif (transposase == "No") and (Abres == "Yes"):
        return "Antibiotic resistance"
    else:
        return "Neither"

def add_genes(valueTable, transTable, abresTable):
    geneTable = pd.merge(transTable[['query', 'transposase']], abresTable[['query', 'antibiotic resistance']], on='query')
    geneTable['Associated Genes'] = geneTable.apply(lambda x: col_add(x.transposase, x['antibiotic resistance']), axis=1)

    return pd.merge(valueTable, geneTable, how='left', left_on='qseqid', right_on='query')

def import_genes(geneFile):
    df = pd.read_table(geneFile, sep="\t", index_col=0)
    df['query'] = df.apply(lambda x: "{}_{}_{}_numOccur:{}".format(x.locus, x.start, x.stop, x['count']), axis=1)
    return df

def full_len_filter(lenDict, valueTable):
    valueTable['qlen'] = valueTable['qseqid'].map(lenDict)
    abreviatedTable = valueTable[valueTable['length'].map(float) >= valueTable['qlen']]
    return abreviatedTable

def seq_lens(seqDict):
    lenDict = {}
    for header in seqDict.keys():
        lenDict[header.lstrip(">")] = len(seqDict[header])
    return lenDict

parser = argparse.ArgumentParser(description="Print fragment abundance")

parser.add_argument("BLASTFile", help="BLAST results file (outfmt=6) to base abundance on")
parser.add_argument("FastaFile", help="Fasta file of query sequences")
parser.add_argument("TransposaseFile", help="Transposase presence table file")
parser.add_argument("AbresFile", help="Antibiotic resistance gene presence table file")
parser.add_argument("-b", "--bins", default=20, type=int, help="Number of bins to divide abundances into")
parser.add_argument("-f", "--figure", default="Frag_abundance_fig.png", type=str, help="Name of output figure")
parser.add_argument("-M", "--binMax", metavar="M", type=float, default=0, help="Explicitly set upper end of bin range")
parser.add_argument("-m", "--binMin", metavar="m", type=float, default=0, help="Explicitly set lower end of bin range")
parser.add_argument("-s", "--binSize", metavar="S", type=float, default=10, help="Set bin size (must also set max)")

args = parser.parse_args()

if args.binMax == 0:
    binBoundaries = args.bins
else:
    binBoundaries = np.arange(args.binMin, args.binMax, args.binSize)

df = pd.read_table(args.BLASTFile, sep="\t")

fastaDict = fol.one_line_d(args.FastaFile)

trans = import_genes(args.TransposaseFile)

abres = import_genes(args.AbresFile)

lenDict = seq_lens(fastaDict)

shortDF = full_len_filter(lenDict, df)

countDF = pd.DataFrame(shortDF.groupby('qseqid').size(),columns=['Occurrences'])
countDF['qseqid'] = countDF.index
geneAddedDF = add_genes(countDF, transTable=trans, abresTable=abres)
#print(geneAddedDF[(geneAddedDF['Occurrences'].map(float) >= 100.0) & (geneAddedDF['Associated Genes'] != 'Transposase')])

fig, ax = plt.subplots(nrows=1, ncols=1)
xcut, xbins = pd.cut(geneAddedDF['Occurrences'], bins=binBoundaries, retbins=True, right=False, precision=2)
binSeries = pd.Series(xbins[1:]).map(float).round(2)
xcut = pd.cut(geneAddedDF['Occurrences'], bins=binBoundaries, right=False, labels=binSeries)
allBins = binSeries.rename('Occurrences').to_frame()
x = geneAddedDF.groupby(['Associated Genes', xcut]).size().unstack('Associated Genes', fill_value=0)
y = pd.merge(x, allBins, right_on='Occurrences', left_index=True, how='right').fillna(0)
del y['Occurrences']
print(y)
y = y.applymap(lambda n: n+1)
y.plot.bar(ax=ax, color=['c', 'm', 'k', 'g'])
ax.set_xticklabels(binSeries)

ax.set_ylabel("Number of mosaic fragments")
ax.set_yscale("log")
ax.set_xlabel("Number of occurrences")
#plt.ylim(ymax=10000)
plt.tight_layout()
plt.savefig(args.figure)
