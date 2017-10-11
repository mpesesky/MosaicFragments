import pandas as pd
from Bio import SeqIO
import argparse

def import_records(filename):
    return SeqIO.index(filename, 'genbank')

def inside(hi, low, start, stop):
    return ((start >= low) and (stop <= hi)) #or ((start <= low) and (stop >= hi))

def sinside(hi, low, i):
    return (i >= low) and (i <= hi)

def overlap(hi, low, start, stop):
    return ((sinside(hi, low, start) and not sinside(hi, low, stop)) or (sinside(hi, low, stop) and not sinside(hi, low, start)))

def keeper(start, end, feature):
    return (inside(end, start, feature.start, feature.end) or (overlap(end, start, feature.start, feature.end)))

def keeper2(start, end, feature):
    return not (((start < feature.start) and (end < feature.start)) or ((start > feature.end) and (end > feature.end)))

def touching_features(locus, start, end, recordDict, pad):
    seqRecord = recordDict[locus]
    featList = seqRecord.features
    features = []

    if pad > start:
        start = 0
    else:
        start = start - pad
    end = end + pad

    for feat in featList:
        if feat.type != 'CDS':
            continue
        if keeper2(start, end, feat.location):
            try:
                features.append(feat.qualifiers['protein_id'][0])
            except KeyError:
                features.append(feat.qualifiers['locus_tag'][0])

    return features

def contains_transposase(idList, transposaseIDs):
    for protid in idList:
        if protid in transposaseIDs:
            return "Yes"
    return "No"


def to_id(locus):
    return locus.split("|")[1]

def to_query(locus, start, stop, occur):
    return "{}_{}_{}_numOccur:{}".format(locus, start, stop, occur)

#def all_connections()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify mosaic fragments associated with transposases")

    parser.add_argument("Fragments", help="Mosaic fragment region file")
    parser.add_argument("GenBank", help="Plasmid Genbank file")
    parser.add_argument("HMM", help="Transposase HMM output")
    parser.add_argument("Outfile", help="Output table name")
    parser.add_argument("-n", "--colName", type=str, default='transposase', help="Name of gene column")
    parser.add_argument("-p", "--pad", type=int, default=0, help="Area around mosaic fragment to check")
    parser.add_argument("-b", "--blastFile", type=str, default=None, help="With -p, BLAST to find all linked plasmids")

    args = parser.parse_args()
    colNames = ['locus', 'start', 'stop', 'length', 'count']

    df = pd.read_table(args.Fragments, sep="\t", header=None, names=colNames)

    df['ID'] = df['locus'].map(to_id)

    gbDict = import_records(args.GenBank)

    hmmdf = pd.read_table(args.HMM, sep="\t")

    transposases = list(hmmdf['query name'].map(to_id))

    if args.pad == 0:
        df['prot_ids'] = df.apply(lambda x: touching_features(x.ID, x.start, x.stop, gbDict, args.pad), axis=1)
        df[args.colName] = df['prot_ids'].map(lambda x: contains_transposase(x, transposases))
#    else:


    df.to_csv(args.Outfile, sep="\t")
