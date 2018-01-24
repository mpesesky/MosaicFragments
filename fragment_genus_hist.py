import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

colNames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                "send", "evalue", "bitscore"]

def name_to_locus(header):
    ret = ""
    try:
        ret = header.split(".")[0].split("|")[1]
    except IndexError:
        print(header)
        exit()
    return ret

def add_genus(blastdf, orgDict):
    blastdf = blastdf[blastdf.qseqid != 'qseqid']
    blastdf['slocus'] = blastdf['sseqid'].map(name_to_locus)
    blastdf['genus'] = blastdf['slocus'].map(orgDict)
    return blastdf

def import_orgs(orgfile):
    orgDict = {}
    orghandle = open(orgfile)
    for line in orghandle:
        fields = line.split("\t")
        orgDict[fields[0]] = fields[1].split(" ")[0]
    orghandle.close()
    return orgDict

def plot_hist(outdf, binBoundaries, ax, outfile, top):
    genusGroups = outdf.groupby('qseqid').genus.nunique()
#    genusGroups += 1
    genusGroups.plot.hist(bins=binBoundaries, ax=ax, color='c')
    plt.gca().set_yscale("log")
    plt.ylim(ymin=0)
    plt.ylabel("Number of Fragments")
    plt.xlabel("Unique Genera")
    plt.tight_layout()
    plt.savefig(outfile)
    return genusGroups.nlargest(top)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plot hist of fragment diversity")

    parser.add_argument("BLAST", help="BLAST file input")
    parser.add_argument("Organisms", help="Organism to locus table file")
    parser.add_argument("-b", "--numbins", metavar="X", type=int, default=20, help="Split into X even bins (default)")
    parser.add_argument("-M", "--binMax", metavar="M", type=float, default=0, help="Explicitly set upper end of bin range")
    parser.add_argument("-m", "--binMin", metavar="m", type=float, default=0, help="Explicitly set lower end of bin range")
    parser.add_argument("-s", "--binSize", metavar="S", type=float, default=10, help="Set bin size (must also set max)")
    parser.add_argument("-o", "--output", type=str, default="FragHist.png", help="Output figure file name")
    parser.add_argument("-t", "--top", type=int, default=10, help="Number of top frags to return")
    parser.add_argument("-c", "--csv", default=None, help="Output csv for +Genus table")

    args = parser.parse_args()

    df = pd.read_table(args.BLAST, sep="\t", header=None, names=colNames)
    organisms = import_orgs(args.Organisms)

    if args.binMax > 0:
        bins = np.arange(args.binMin, args.binMax, args.binSize)
    else:
        bins = args.numbins

    df = add_genus(df, organisms)

    if args.csv is not None:
        df.to_csv(args.csv, sep="\t")



    fig, ax = plt.subplots(nrows=1, ncols=1)

    topX = plot_hist(df, bins, ax, args.output, args.top)


    print(topX)