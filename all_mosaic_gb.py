import pandas as pd
from Bio import SeqIO
import argparse


def blast_to_positions(blastTable):
    return blastTable[['qseqid', 'qstart', 'qend', 'sseqid', 'sstart', 'send']]


def accession_to_locus(acc):
    return acc.split("|")[1]


parser = argparse.ArgumentParser(description="Get all mosaic fragments from a genbank file")

parser.add_argument("Genbank", help="Genbank file to split")
parser.add_argument("MosaicBlast", help="BLAST output table indicating mosaic fragments")
parser.add_argument("Outfile", help="Desired name of gene database")

args = parser.parse_args()

blastdf = pd.read_table(args.MosaicBlast, sep="\t")
positionFrame = blast_to_positions(blastdf)
gb = SeqIO.index(args.Genbank, 'genbank')

positionFrame['qlocus'] = positionFrame['qseqid'].map(accession_to_locus)
positionFrame['slocus'] = positionFrame['sseqid'].map(accession_to_locus)
geneFrame = pd.DataFrame(columns=['qlocus', 'slocus', 'geneID', 'description'])

for i, row in positionFrame.iterrows():
    record = gb[row['qlocus']]
    mosaic = record[int(row['qstart']):int(row['qend'])]
    qList = []
    sList = []
    idList = []
    descList = []

    for feature in mosaic.features:
        if (feature.type == 'CDS') and ('pseudo' not in feature.qualifiers):
            qList.append(row['qlocus'])
            sList.append(row['slocus'])
            try:
                idList.append(feature.qualifiers['protein_id'][0])
            except KeyError:
                print(feature)
                exit()
            try:
                descList.append(feature.qualifiers['product'][0])
            except KeyError:
                print(feature)
                exit()
    geneDict = {'qlocus': qList, 'slocus':sList, 'geneID':idList, 'description':descList}
    tempFrame = pd.DataFrame(data=geneDict)
    geneFrame = geneFrame.append(tempFrame, sort=True)

geneFrame.drop_duplicates(inplace=True)
geneFrame.to_csv(args.Outfile, sep="\t")

